Module ModRamScb
! Contains subroutines and functions necessary to couple RAM and SCB

  use nrtype, ONLY: DP

  implicit none
  save

  INTEGER, allocatable :: INDEXPA(:,:,:,:) ! Index that maps the pitch angle on each point of a field line to the equatorial (RAM) one
  REAL(DP), allocatable :: flux3DEQ(:,:,:,:,:)
  INTEGER :: NRAD
  INTEGER :: INCFD, j, k, L
  REAL(DP), allocatable :: bfMirror(:), bfInterm(:)
  REAL(DP), allocatable :: chiMirror(:)
  REAL(DP), allocatable :: derivs2(:), dBdTheta(:), dThetadB(:), derivs4(:), &
                           distance(:), dyDummyDens(:)

contains

!==============================================================================
subroutine ramscb_allocate

  use ModRamGrids, ONLY: NS, NE, NPA
  use ModScbGrids, ONLY: nthe, npsi, nzeta, nXRaw

  implicit none

  ALLOCATE(INDEXPA(nthe,npsi,nzeta,NPA), Flux3DEq(NS,npsi,nzeta,NE,NPA), &
           bfMirror(NPA), bfInterm(NPA), chiMirror(NPA), derivs2(nthe), &
           dBdTheta(nthe), dThetadB(nthe), derivs4(nthe), distance(nthe), &
           dyDummyDens(nthe))
  
  NRAD = nXRaw
end subroutine ramscb_allocate


!==============================================================================
subroutine ramscb_deallocate

  implicit none

  DEALLOCATE(INDEXPA, Flux3DEq, bfMirror, bfInterm, chiMirror, derivs2, &
             dBdTheta, dThetadB, derivs4, distance, dyDummyDens)

end subroutine ramscb_deallocate

!==============================================================================
SUBROUTINE computehI
  !****************************************************************************
  ! Routine that outputs magnetic field quantities for input into RAM:
  ! h and I integrals, plus equatorial Bz
  ! values output at radOut; h, I values at muboun are output also, 
  ! as well as bounce-averaged HDENS
  !
  ! Added branch for DIPL (operational run with dipole field only)
  ! still need from here:
  !  - indexPA (once)
  !  - flux3D (every 300 s or the satellite dump time cadence)
  ! Author: S. Zaharia, 2011
  !
  ! Added flux volume per unit flux (integral of ds/B) output (Apr. 2013)
  !
  !  Copyright (c) 2016, Los Alamos National Security, LLC
  !  All rights reserved.
  !****************************************************************************

  use ModRamVariables, ONLY: Kp, F107, FLUX, FNHS, FNIS, BNES, HDNS, dBdt, &
                             dIdt, dIbndt, BOUNHS, BOUNIS, EIR, EIP, DPHI, PHI, &
                             RLZ, LZ, MU, DMU, WMU, MLT, PAbn, PA
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamStart, DT_hI, TimeRamElapsed, TOld
  use ModRamConst,     ONLY: RE, HMIN, ME
  use ModRamParams,    ONLY: electric, IsComponent, NameBoundMag
  use ModRamCouple,    ONLY: SwmfPot_II
  use ModRamGrids,     ONLY: nR, nT, nPa, nE, RadiusMax

  USE ModScbMain,      ONLY: rHour, HourDecimal, Day, HourDecShort, iElectric, &
                             prefixOut
  use ModScbGrids,     ONLY: nthe, npsi, nzeta, nXRaw, nAzimRAM
  use ModScbVariables, ONLY: bf, chiVal, nZetaMidnight, x, y, z, bnormal, PhiIono, &
                             radGrid, angleGrid, h_Cart, I_Cart, h_Cart_interp, &
                             i_Cart_interp, bZEq_Cart, flux_vol_cart, bZ, chi, &
                             hdens_Cart, radRaw, azimRaw, nThetaEquator, fluxVolume

  use ModRamFunctions, ONLY: RamFileName, COSD, FUNT, FUNI

  use ModScbSpline, ONLY: Spline_2D_point, spline, splint
  use ModScbInterp, ONLY: Interpolation_natgrid_2D

  USE nrmod,     ONLY: qtrap, qromb, locate
  use nrtype,    ONLY: DP, pi_d
  USE ModIOUnit, ONLY: UNITTMP_

  IMPLICIT NONE

  REAL(DP) :: BNESPrev(NR+1,NT), FNISPrev(NR+1,NT,NPA), BOUNISPrev(NR+1,NT,NPA)
  INTEGER :: I, IC, iCount, ierr, IP1, j1, k1, NWK, ncdfId, iFix, iDomain
  INTEGER, DIMENSION(3) :: dimlens 
  INTEGER, PARAMETER :: LIMIT = 10000, INF = 1
  REAL(DP) :: MUEQ
  REAL(DP), DIMENSION(2*nthe) :: WK
  REAL(DP) :: switch, start_time, end_time, yLowerLimit, chiHalf, bfHalf
  REAL(DP) :: r0(npsi, nzeta)
  REAL(DP) :: MUBOUN, RWU, DL1, CLC, dyDummy
  INTEGER, PARAMETER :: LENIW = LIMIT, LENW = 4*LENIW+1
  REAL(DP) :: bZEqDiff_Cart(nR,nT), Epot_Cart(0:nR,nT), beqdip_Cart(nR)
  INTEGER :: iwork(LENIW)  
  INTEGER :: NEVAL, LAST, myAlphaBegin, myAlphaEnd, nratio, alphaRangeDiff
  REAL(DP) :: work(LENW)
  INTEGER, PARAMETER :: NPTS2 = 3  ! # of break points plus 2
  REAL(DP) :: POINTS(NPTS2)
  REAL(DP) :: EPSABS, EPSREL ! If EPSABS < 0, only EPSREL is taken into account; decide on what REL error you need by 
  ! comparing values obtained for a dipole with Ejiri's values
  REAL(DP) :: ABSERR, length, valueIntegral, valueintegral2, valueintegraldens, valueintegraldens2, &
       valueintegraldens3, valueIntegral2Bis, valueintegral3
  INTEGER :: IC0(2)
  INTEGER :: alphaBeg(0:63), alphaEnd(0:63), param(0:63), my_array_type_hI(0:63) ! Assumes no more than 64 processors
  INTEGER :: sizes(3), subsizes(3), starts(3)
  INTEGER :: request(3, 0:63) 
  REAL(DP) :: VC(2)
  REAL(DP) :: h_value(npsi, nzeta, NPA), hdens_value(npsi, nzeta, NPA), I_value(npsi, nzeta, NPA), I_value_Yue(npsi,nzeta,NPA) 
  REAL(DP) :: dyDummy1(NPA), dyDummy2(NPA)

  CHARACTER(len=500) :: fileNamehi
  character(len=200) :: NameFileOut
  CHARACTER(len=2) :: DayChar
  CHARACTER(len=4) :: ST3_local, ST3_local2, ST3_local3, ST3_el
  CHARACTER(LEN=100) :: HEADER(4)
  INTEGER :: nm_local
  LOGICAL :: SKIP = .FALSE.
  REAL(DP) :: radOut = RadiusMax+0.25_dp ! 9.25_dp ! 6.75_dp ! 10.25_dp ! Outer domain radius

  INTEGER :: iSpecies
  REAL(DP) :: EKEV2(NE), MU2(NPA)
  REAL(DP), PARAMETER :: b0dip = 30570._dp
  REAL(DP) :: beqdip(npsi,nzeta+1)
  REAL(DP) :: DthI
!  REAL(DP), :: flux_volume(:,:)

  switch = 0.0_dp
  dimlens = (/1, 1, 1/)

  EPSABS = 0.0006_dp
  EPSREL = 1.E-6_dp  ! 1% error should be fine for most applications

  h_Cart = 0._dp
  I_Cart = 0._dp
  bZEq_Cart = 0._dp
  flux_vol_Cart = 0._dp
  bZEqDiff_Cart = 0._dp

  h_value = 0._dp
  hdens_value = 0._dp
  I_value = 0._dp
  I_value_Yue = 0._dp
  valueIntegral = 0._dp
  valueIntegral2 = 0._dp
  valueIntegraldens = 0._dp
  valueIntegral3 = 0._dp

  yLowerLimit = 0._dp

  DO j1 = 0, nR
     radRaw(j1) = 1.75_dp + (radOut-1.75_dp) * REAL(j1,DP)/REAL(nR,DP)
  END DO
  DO k1 = 1, nT
     azimRaw(k1) = 24._dp * REAL(k1-1,DP)/REAL(nT-1,DP)
  END DO

  INCFD = 1
  NWK = 2*(nthe-nThetaEquator+1)

  ! Always fill this matrix; it's used by RAM outputs
  ! Interpolate RAM flux on 3DEQ grid, for mapping
  DO iSpecies = 1, 4 ! ions & electrons
     DO L = 1, NPA
        DO K = 1, NE
           !Cubic spline interpolation
           CALL Spline_2D_point(radRaw(1:nR), azimRaw*pi_d/12.0, &
                                REAL(FLUX(iSpecies, 1:nR,1:nT,K,L),DP), &   
                                radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), &
                                flux3DEQ(iSpecies, 1:npsi,2:nzeta,K,L), iDomain) 
        END DO
     END DO
  END DO
  ! Periodicity
  flux3DEQ(:,:,1,:,:) = flux3DEQ(:,:,nzeta,:,:)

  Operational_or_research: SELECT CASE (NameBoundMag)
  CASE('DIPL') ! Dipole without SCB calculation

!     IF (nm_local > 1) RETURN ! indexPA only needs to be calculated once

     Alpha_loop_parallel_oper:  DO k = 2, nzeta 
        Psi_loop_oper: DO j = 1, npsi 
           ! If bf(theta) not strictly increasing; to weed out very small differences  
           ! Generally if needed it means have to increase number of theta points (or crowd them more near the equatorial plane)
           ! This might happen for DIPL once in a blue moon if ntheta is ridiculously small
           iCount = 0
           Monotonicity_oper: DO 
              i = nThetaEquator
              iFix = 0
              iCount = iCount+1
              fixMonotonicity_oper: DO WHILE (i < nthe)      
                 IF (bf(i+1,j,k) < bf(i,j,k)) THEN
                    bf(i+1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                    iFix = 1   
                 END IF
                 i = i+1
              END DO fixMonotonicity_oper
              IF (iFix == 0) EXIT Monotonicity_oper
              CYCLE Monotonicity_oper
           END DO Monotonicity_oper

           ! Define indexPA for 90 degree pitch angle (-1 for all distances along field line except equatorial plane)
           indexPA(:,:,:,1) = -1
           indexPA(nThetaEquator,:,:,1) = 1

           pitchAngleLoop_oper:        DO L = 2, NPA 
              ! Find mirror points thetaM for mu, where B = Beq/(1 - mu**2)
              IF (L /= NPA) bfMirror(L) = bf(nThetaEquator,j,k) / (1._dp - mu(L)*mu(L))
              IF (L == NPA) bfMirror(L) = bf(nthe,j,k)

              IF (bfMirror(L) > bf(nthe,j,k)) THEN ! No possibility of trapping, only passing particles
                 iCount = iCount + 1
                 bfMirror(L) = bf(nthe,j,k)
              END IF

              DO i = nThetaEquator, nthe
                 IF (1. - bf(nThetaEquator,j,k)/bf(i,j,k)*(1.-MU(L)**2) >= 0.) THEN
                    MUEQ = SQRT(1. - bf(nThetaEquator,j,k)/bf(i,j,k)*(1.-MU(L)**2))
                    indexPA(i,j,k,L) = locate(MU(1:NPA), MUEQ)
                    indexPA(nthe+1-i,j,k,L) = indexPA(i,j,k,L) ! Mirror across magnetic equator.
                 ELSE
                    indexPA(i,j,k,L) = -1
                    indexPA(nthe+1-i,j,k,L) = -1
                 END IF
              END DO

           END DO pitchAngleLoop_oper
        END DO Psi_loop_oper
     END DO Alpha_loop_parallel_oper

     RETURN 

  CASE default  ! Regular calculation of h, I, bounce-averaged charge xchange etc.
     Alpha_loop_parallel:  DO k = 2, nzeta 
        Psi_loop: DO j = 1, npsi ! nZetaMidnight, nZetaMidnight !C 2, nzeta

           r0(j, k) = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
           DO i = nThetaEquator, nthe
              distance(i) = SQRT(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2) ! Distance from center of Earth
           END DO
           !CALL spline(chiVal(nThetaEquator:nthe),distance(nThetaEquator:nthe), 1.E31_dp, 1.E31_dp, &
           !            dyDummyDens(nThetaEquator:nthe))
           CALL spline(chi(nThetaEquator:nthe,j,k),distance(nThetaEquator:nthe), 1.E31_dp, 1.E31_dp, &
                       dyDummyDens(nThetaEquator:nthe))

           ! Compute length of field line (j,k)
           length = 0._dp
           DO i = 1, nthe-1
              length = length + SQRT((x(i+1,j,k)-x(i,j,k))**2 + (y(i+1,j,k)-y(i,j,k))**2 + (z(i+1,j,k)-z(i,j,k))**2)
           END DO

           ! Boundary choices at the ends (look into DPCHIC etc. for explanation)
           IC0(1) = 0 
           IC0(2) = -5 

           ! If bf(theta) not strictly increasing; to weed out very small differences (which nevertheless will make DPCHIC return an error and stop) 
           ! Generally if needed it means have to increase number of theta points (or crowd them more near the equatorial plane)
           iCount = 0
           Monotonicity: DO 
              i = nThetaEquator
              iFix = 0
              iCount = iCount+1
              fixMonotonicity: DO WHILE (i < nthe)      
                 IF (bf(i+1,j,k) < bf(i,j,k)) THEN
                    bf(i+1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                    iFix = 1   
                 END IF
                 i = i+1
              END DO fixMonotonicity
              IF (iFix == 0) EXIT Monotonicity
              CYCLE Monotonicity
           END DO Monotonicity

           ! Initialize Hermite spline for chi(B), in order to find thetaMirror; actually, chiMirror in the general case (constTheta /= 0)
           !CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe),  &
           !     dThetadB(nThetaEquator:nthe), INCFD, WK(1:NWK), NWK, IERR)
           CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, bf(nThetaEquator:nthe,j,k), chi(nThetaEquator:nthe,j,k),  &
                dThetadB(nThetaEquator:nthe), INCFD, WK(1:NWK), NWK, IERR)
           
           IF (ierr < 0) THEN ! Means monotonicity of theta(B) not established
              PRINT*, 'FATAL error in computehI: ierr = ', ierr
              STOP
           END IF

           iCount = 0

           ! Define indexPA for 90 degree pitch angle (-1 for all distances along field line except equatorial plane)
           indexPA(:,:,:,1) = -1
           indexPA(nThetaEquator,:,:,1) = 1

           ! L = 1 is the 90 degree pitch angle (mu = 0),  L = NPA is the 0 degree pitch angle (mu = 1)
           pitchAngleLoop:        DO L = 2, NPA ! No need to do anything for L = 1 (equator => I = 0; h not zero, but obtained by extrapolation or just 
              ! set the same as h(L=2) - see below)

              ! Find mirror points thetaM for mu, where B = Beq/(1 - mu**2)
              IF (L /= NPA) bfMirror(L) = bf(nThetaEquator,j,k) / (1._dp - mu(L)*mu(L))
              IF (L == NPA) bfMirror(L) = bf(nthe,j,k)

              IF (bfMirror(L) > bf(nthe,j,k)) THEN ! No possibility of trapping, only passing particles
                 iCount = iCount + 1
                 bfMirror(L) = bf(nthe,j,k)
                 !chiMirror(L) = chiVal(nthe)
                 chiMirror(L) = chi(nthe,j,k)
              ELSE
                 ! Find chiMirror, assuming monotonicity of bf with chiVal
                 !CALL DPCHFE (nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
                 !     dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bfMirror(L), chiMirror(L), IERR)
                 CALL DPCHFE (nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chi(nThetaEquator:nthe,j,k), &
                      dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bfMirror(L), chiMirror(L), IERR)
              END IF

              DO i = nThetaEquator, nthe
                 IF (1. - bf(nThetaEquator,j,k)/bf(i,j,k)*(1.-MU(L)**2) >= 0.) THEN
                    MUEQ = SQRT(1. - bf(nThetaEquator,j,k)/bf(i,j,k)*(1.-MU(L)**2))
                    indexPA(i,j,k,L) = locate(MU(1:NPA), MUEQ)
                    indexPA(nthe+1-i,j,k,L) = indexPA(i,j,k,L) ! Mirror across magnetic equator.
                 ELSE
                    indexPA(i,j,k,L) = -1
                    indexPA(nthe+1-i,j,k,L) = -1
                 END IF
              END DO

              ! Compute integral for mu
              IC0(1) = -1   ! First derivative dB/dchi set to zero at equatorial plane
              VC(1) = 0._dp
              IC0(2) = -5 

              ! Initialize Hermite spline B(theta) (to be used in fInt and fScalarInt)
              !CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
              !     dBdTheta(nThetaEquator:nthe), INCFD, WK(1:NWK), NWK, IERR)
              CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, chi(nThetaEquator:nthe,j,k), bf(nThetaEquator:nthe,j,k), &
                   dBdTheta(nThetaEquator:nthe), INCFD, WK(1:NWK), NWK, IERR)

              IF (ierr < 0) THEN
                 PRINT*, 'computehI2: ierr = ', ierr
                 STOP
              END IF

              !valueIntegral = qromb(fInt, chiVal(nThetaEquator), chiMirror(L))  ! The integral for I is non-singular, qromb should do fine
              valueIntegral = qromb(fInt, chi(nThetaEquator,j,k), chiMirror(L))

              ! Split h integral or compute solely by DQAGI (which is pretty good for the whole interval)
              ! March 2006 - transformed h integral by y = 1/sqrt(bm-b) change of variable; spreads out the interval from yeq to infinity; not singular anymore 
              ! ---> better accuracy

              !chiHalf = 0.2_dp*(chiVal(nThetaEquator) + 4.*chiMirror(L))
              !CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
              !     dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chiHalf, bfHalf, IERR)          
              chiHalf = 0.2_dp*(chi(nThetaEquator,j,k) + 4.*chiMirror(L))
              CALL DPCHFE (nthe-nThetaEquator+1, chi(nThetaEquator:nthe,j,k),bf(nThetaEquator:nthe,j,k), &
                   dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chiHalf,bfHalf, IERR) 

              IF (bfMirror(L) - bfHalf < 0._dp) THEN
                 PRINT*, 'computehI: problem w/ mirror point; j, k, L = ', j, k, L
                 !   STOP
              END IF

              yLowerLimit = 1._dp/SQRT(bfMirror(L)-bfHalf)

              CALL DQAGI(fScalarIntInf, yLowerLimit, INF, EPSABS, EPSREL, valueIntegral2, &
                   ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! DQAGI is the transformed integral, employed in the vicinity of the mirror point

              CALL DQAGI(fScalarIntInfDens, yLowerLimit, INF, EPSABS, EPSREL, valueIntegraldens2, &
                   ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! DQAGI is the transformed integral (for Hdens), employed in the vicinity of the mirror point

              valueIntegral2 = 2._dp * valueIntegral2 ! Factor of 2 comes from the conversion of integral over chi to integral over y
              valueIntegralDens2 = 2._dp * valueIntegralDens2

              ! Previous way of doing the h integral - good but the DQAGI method is better
              ! DQAGS - from Slatec library - adaptive integration with end-point singularities - double precision - for the integral for h
              !C CALL DQAGS(fScalarInt, chiVal(nThetaEquator), chiMirror(L), EPSABS, EPSREL, valueIntegral3, &
              !C     ABSERR, NEVAL, IERR, LENIW, LENW, LAST, IWORK, WORK)  ! If getting error messages, it means the field line is not "smooth" enough - usually
              ! need a new calculation with a smaller blending coefficient; IERR = 4 however just means that the required accuracy cannot be attained, but the result
              ! however is the best that can be obtained 

              ! Replaced this with same integral splitting for Hdens (DQAGS was giving some Inf values)
              !         CALL DQAGS(fScalarIntDens, chiVal(nThetaEquator), chiMirror(L), EPSABS, EPSREL, valueIntegralDens, &
              !              ABSERR, NEVAL, IERR, LENIW, LENW, LAST, IWORK, WORK)   ! whole interval for density

              !CALL DQAG(fScalarInt, chiVal(nThetaEquator), chiHalf, EPSABS, EPSREL, 3, valueIntegral3, &
              !     ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! from eq. plane to chiHalf; regular integral (for h calculation)
              !CALL DQAG(fScalarIntDens, chiVal(nThetaEquator), chiHalf, EPSABS, EPSREL, 3, valueIntegralDens3, &
              !     ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! from eq. plane to chiHalf; regular integral of Hdens
              CALL DQAG(fScalarInt, chi(nThetaEquator,j,k), chiHalf, EPSABS, EPSREL, 3, valueIntegral3, &
                   ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! from eq. plane to chiHalf; regular integral (for h calculation)
              CALL DQAG(fScalarIntDens, chi(nThetaEquator,j,k), chiHalf, EPSABS, EPSREL, 3, valueIntegralDens3, &
                   ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! from eq. plane to chiHalf; regular integral of Hdens


              valueIntegral2 = valueIntegral2 + valueIntegral3
              valueIntegralDens2 = valueIntegralDens2 + valueIntegralDens3
              hdens_value(j,k,L) = 1.E10_dp * valueIntegralDens2/valueIntegral2  ! Re-normalize (see definition in hdens_rairden)

              valueIntegral2 = valueintegral2 * SQRT(bfMirror(L))  ! Normalization for calculation of h

              I_value(j, k, L) = 2._dp * valueIntegral * length / (pi_d * r0(j, k))
              h_value(j, k, L) = valueIntegral2 * length / (pi_d * r0(j, k))

           END DO pitchAngleLoop

           ! Find h for PA = 90 degrees - various extrapolation tried, but better keep value at NPA=2 (same thing done in RAM code)
           I_value(j, k, 1) = 0._dp
           h_value(j, k, 1) = h_value(j, k, 2) ! This extrapolation is also used in the RAM code
           hdens_value(j, k, 1) = hdens_value(j, k, 2)

        END DO Psi_loop
     END DO Alpha_loop_parallel

     ! Periodicity for indexPA ! Alpha_loop_parallel only from 2 to nzeta
     indexPA(:,:,1,:) = indexPA(:,:,nzeta,:)

     DO j = 1, npsi
        DO k = 2, nzeta 
           I_value_Yue(j,k, :) = I_value(j,k, :)*r0(j,k)
        END DO
     END DO
     h_value(:,1,:) = h_value(:,nzeta,:)
     I_value(:,1,:) = I_value(:,nzeta,:)
     I_value_Yue(:,1,:) = I_value_Yue(:,nzeta,:)

     ! Interpolation will be for B - Bdip
     beqdip(1:npsi,2:nzeta)=b0dip/(x(nThetaEquator,1:npsi,2:nzeta)**2+y(nThetaEquator,1:npsi,2:nzeta)**2)**1.5

        ! Interpolate data for output in POLAR coordinates (for RAM)

     IF (MINVAL(h_value) < 0._dp .OR. MINVAL(I_value)<0._dp) THEN
        PRINT*, 'hI: minval(h) = ', MINVAL(h_value)
        PRINT*, 'hI: minval(I) = ', MINVAL(I_value)
        !   STOP
     END IF

     CALL Interpolation_natgrid_2D(radRaw(1:nR), azimRaw, h_value(1:npsi,2:nzeta,1:NPA), &
                I_value(1:npsi,2:nzeta,1:NPA), bZ(nThetaEquator,1:npsi,2:nzeta)*bnormal-beqdip(1:npsi,2:nzeta), &
                fluxVolume(1:npsi,2:nzeta)/bnormal, hdens_value(1:npsi,2:nzeta,1:NPA), &
                h_Cart(1:nR,1:nT,1:NPA), I_Cart(1:nR,1:nT,1:NPA), &
                bZEqDiff_Cart(1:nR,1:nT), flux_vol_Cart(1:nR,1:nT), &
                hdens_Cart(1:nR,1:nT,1:NPA))

     ! Add Bdip 
     beqdip_Cart(1:nR) = b0dip/radRaw(1:nR)**3

     ! If 3D code domain does not extend to radOut: (ACHTUNG!!!)
     ! print*, 'computehI: max(x), radOut = ', real(maxval(abs(x)),sp), real(radOut,SP)
     IF (MAXVAL(ABS(x)) <= radOut) THEN
        h_Cart(nR,:,:) = h_Cart(nR-1,:,:)
        I_Cart(nR,:,:) = I_Cart(nR-1,:,:)
        bZEqDiff_Cart(nR,:) = bZEqDiff_Cart(nR-1,:)
        hdens_Cart(nR,:,:) = hdens_Cart(nR-1,:,:)
     END IF

     ! Add dipole field
     DO j = 1, nT
        bZEq_Cart(:,j)  = bZEqDiff_Cart(:,j) + beqdip_Cart
     END DO

     ! Make field quasi-dipole as R -> 2 RE
     DO i = 1, nR
        IF (radRaw(i) <= 3.0) bZEq_Cart(i,:) = EXP(-2.*(radRaw(i)-3.0)**2)*bZEq_Cart(i,:) + &
             (1.-EXP(-2.*(RadRaw(i)-3.0)**2))*beqdip_Cart(i)
     END DO

     IF (MINVAL(h_Cart) < 0._dp .OR. MINVAL(I_Cart)<0._dp) THEN
        PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart)
        PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart)
        !   STOP
     END IF

     ! Cubic spline interpolation with natural boundaries to get h and I at muboun 
     DO j = 1, nT
        DO i = 1, nR
           CALL spline(PA(1:NPA),h_Cart(i,j,1:NPA),1.E31_dp,1.E31_dp,dyDummy1(1:NPA))
           CALL spline(PA(1:NPA),I_Cart(i,j,1:NPA),1.E31_dp,1.E31_dp,dyDummy2(1:NPA))
           DO L = 1, NPA-1
              h_Cart_interp(i,j,L) = splint(PA(1:NPA),h_Cart(i,j,1:NPA),dyDummy1(1:NPA),PAbn(L))
              I_Cart_interp(i,j,L) = splint(PA(1:NPA),I_Cart(i,j,1:NPA),dyDummy2(1:NPA),PAbn(L))
           END DO
           ! Do not do the NPA inclusive in the not-a-knot interpolation above -> can lead to negative h,I(NPA)
           h_Cart_interp(i,j,NPA) = h_Cart_interp(i,j,NPA-1)
           I_Cart_interp(i,j,NPA) = I_Cart_interp(i,j,NPA-1)
        END DO
     END DO

     ! Update h and I values for RAM (note that the output of the hI files
     ! now uses RAM variables so the numbers will be different
     DthI = TimeRamElapsed-TOld
     DO I=2,NR+1
        DO J=1,NT
           BNESPrev(I,J)=BNES(I,J)
           DO L=1,NPA
              FNISPrev(I,J,L)=FNIS(I,J,L)
              BOUNISPrev(I,J,L)=BOUNIS(I,J,L)
              ! SCB Variables -> RAM Variables
              FNHS(I,J,L)   = h_Cart(I-1,J,L)
              FNIS(I,J,L)   = i_Cart(I-1,J,L)
              BNES(I,J)     = bZEq_Cart(I-1,J)
              HDNS(I,J,L)   = hdens_Cart(I-1,J,L)
              BOUNHS(I,J,L) = h_Cart_interp(I-1,J,L)
              BOUNIS(I,J,L) = i_Cart_interp(I-1,J,L)
              !
              dIdt(I,J,L)=(FNIS(I,J,L)-FNISPrev(I,J,L))/DThI
              dIbndt(I,J,L)=(BOUNIS(I,J,L)-BOUNISPrev(I,J,L))/DThI
           ENDDO
           IF (J.LT.8.or.J.GT.18) THEN
              DO L=15,2,-1
                 if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
                    FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
                 endif
              ENDDO
           ENDIF
           BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
           dBdt(I,J)=(BNES(I,J)-BNESPrev(I,J))/DThI
        ENDDO
     ENDDO
     DO J=1,NT ! use dipole B at I=1
        BNES(1,J)=0.32/LZ(1)**3/1.e4
        dBdt(1,J) = 0.
        EIR(1,J) = 0.
        EIP(1,J) = 0.
        DO L=1,NPA
           FNHS(1,J,L) = FUNT(MU(L))
           FNIS(1,J,L) = FUNI(MU(L))
           BOUNHS(1,J,L)=FUNT(COSD(PAbn(L)))
           BOUNIS(1,J,L)=FUNI(COSD(PAbn(L)))
           HDNS(1,J,L)=HDNS(2,J,L)
           dIdt(1,J,L)=0.
           dIbndt(1,J,L)=0.
        ENDDO
     ENDDO
     TOld = TimeRamElapsed

!     END IF Master_rank

15   FORMAT(A, F6.1)
21   FORMAT(F8.2, F10.1, 3X, I2, 5X, 8(3X, E12.4))

  END SELECT Operational_or_research

  RETURN

END SUBROUTINE computehI

!==============================================================================
FUNCTION fInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator, chi
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local(:) 
  REAL(DP) :: bf_local(SIZE(chi_local))
  REAL(DP) :: dyDummy(SIZE(chi_local))
  INTEGER :: i, ierr
  REAL(DP) :: fInt(SIZE(chi_local))
  LOGICAL :: SKIP = .FALSE.

  !CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
  !     dBdTheta(nThetaEquator:nthe), INCFD, SKIP, SIZE(chi_local,1), chi_local, bf_local, IERR)
  CALL DPCHFE (nthe-nThetaEquator+1, chi(nThetaEquator:nthe,j,k), bf(nThetaEquator:nthe,j,k), &
       dBdTheta(nThetaEquator:nthe), INCFD, SKIP, SIZE(chi_local,1), chi_local, bf_local, IERR)


  !C fInt = SQRT(MAX(1._dp/(1._dp - bf_local/bfMirror(L)), 0._dp)) ! For function h(mu0)
  fInt = SQRT(MAX((bfMirror(L) - bf_local)/bfMirror(L), 0._dp)) ! For function I(mu0)

  RETURN
END FUNCTION fInt

!==============================================================================
FUNCTION fIntDer(t_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: bf, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: t_local(:) 
  REAL(DP) :: bf_local(SIZE(t_local)), bfprime_local(SIZE(t_local))
  REAL(DP) :: dyDummy(SIZE(t_local))
  INTEGER :: i, ierr
  REAL(DP) :: fIntDer(SIZE(t_local))

  bf_local = bfMirror(L) - t_local*t_local

  ! Find Bprime for the current argument array
  DO i = 1, SIZE(t_local)
     bfprime_local(i) = splint(bf(nThetaEquator:nthe,j,k), dBdTheta(nThetaEquator:nthe), &
          derivs4(nThetaEquator:nthe), bf_local(i))
  END DO

  fIntDer = 1._dp / bfprime_local

  RETURN
END FUNCTION fIntDer

!==============================================================================
FUNCTION fScalarIntDens(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: x, y, z, bf, chiVal, bnormal, nThetaEquator, chi
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: dyDummy(nthe), radius
  REAL(DP) :: rad(nthe)
  REAL(DP), PARAMETER :: HUGE = 1.E12_dp
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntDens
  LOGICAL :: SKIP = .FALSE.

  !CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
  !     dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chi_local, bf_local, IERR)
  !radius = splint(chiVal(nThetaEquator:nthe), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)
  CALL DPCHFE (nthe-nThetaEquator+1, chi(nThetaEquator:nthe,j,k), bf(nThetaEquator:nthe,j,k), &
       dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chi_local, bf_local, IERR)
  radius = splint(chi(nThetaEquator:nthe,j,k), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)


  fScalarIntDens = hdens_rairden(radius) * SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) 

  RETURN

END FUNCTION fScalarIntDens

!==============================================================================
FUNCTION fScalarInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator, chi
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: dyDummy
  REAL(DP), PARAMETER :: HUGE = 1.E12_dp
  INTEGER :: i, ierr
  REAL(DP) :: fScalarInt
  LOGICAL :: SKIP = .FALSE.

  !CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
  !     dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chi_local, bf_local, IERR)
  CALL DPCHFE (nthe-nThetaEquator+1, chi(nThetaEquator:nthe,j,k), bf(nThetaEquator:nthe,j,k), &
       dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chi_local, bf_local, IERR)


  !C if (chi_local > chiVal(nthe)) STOP 'Problem in fScalarInt.'

  IF (IERR < 0) THEN  ! IERR = 0 normal return, IERR > 0 non-fatal error
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  fScalarInt = SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) ! For function h(mu0)
!C  if (bf_local < bfMirror(L)) then
!C     fScalarInt = SQRT(1._dp/(bfMirror(L) - bf_local)) ! For function h(mu0)
!C  else
!C     fScalarInt = 0._dp
!C  end if

  RETURN

END FUNCTION fScalarInt

!==============================================================================
FUNCTION fScalarIntInf(y_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator, chi
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: y_local
  REAL(DP) :: bf_local, chi_local, dChidB_local
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntInf
  LOGICAL :: SKIP = .FALSE.

  bf_local = bfMirror(L) - 1._dp/(y_local*y_local)

  !CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
  !     dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bf_local, chi_local, dChidB_local, IERR)
  CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chi(nThetaEquator:nthe,j,k), &
       dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bf_local, chi_local, dChidB_local, IERR)


  IF (IERR < 0) THEN  
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  IF(y_local*y_local>0._dp) THEN
     fScalarIntInf = dChidB_local/(y_local*y_local)
  ELSE
     PRINT*, 'Problem in fScalarIntInf; y_local = ', y_local
  END IF

  RETURN
END FUNCTION fScalarIntInf

!==============================================================================
FUNCTION fScalarIntInfDens(y_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  USE ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator, chi
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: y_local
  REAL(DP) :: bf_local, chi_local, dChidB_local, radius
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntInfDens
  LOGICAL :: SKIP = .FALSE.

  bf_local = bfMirror(L) - 1._dp/(y_local*y_local)

  !CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
  !     dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bf_local, chi_local, dChidB_local, IERR)
  CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chi(nThetaEquator:nthe,j,k), &
       dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bf_local, chi_local, dChidB_local, IERR)


  IF (IERR < 0) THEN  
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  !radius = splint(chiVal(nThetaEquator:nthe), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)
  radius = splint(chi(nThetaEquator:nthe,j,k), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)

  IF(y_local*y_local>0._dp) THEN
     fScalarIntInfDens = dChidB_local/(y_local*y_local) * hdens_rairden(radius)
  ELSE
     PRINT*, 'Problem in fScalarIntInf; y_local = ', y_local
  END IF

  RETURN
END FUNCTION fScalarIntInfDens

!==============================================================================
FUNCTION hdens_rairden(radius)
  USE nrtype, ONLY : DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: radius
  REAL(DP) :: rad_local
  REAL(DP) :: hdens_rairden
  REAL(DP), PARAMETER :: a0=13.326_dp, a1=-3.6908_dp, a2=1.1362_dp, &
       a3=-0.16984_dp, a4=0.009552_dp

  rad_local = MIN(radius, 6.4_dp)  ! Rairden function increases for R > 6.4
  hdens_rairden = a0+a1*rad_local+a2*rad_local**2 + a3*rad_local**3 + a4*rad_local**4
  hdens_rairden = 1.E-10_dp * 10._dp**hdens_rairden  ! Normalized by 10^10 to allow integrals to be not too large

  RETURN

END FUNCTION hdens_rairden

END MODULE ModRamScb
