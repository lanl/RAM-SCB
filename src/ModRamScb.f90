Module ModRamScb
! Contains subroutines and functions necessary to couple RAM and SCB

  use nrtype, ONLY: DP
  use ModRamGrids, ONLY: NS, NT, NR, NE, NPA, RadiusMax
  use ModScbGrids, ONLY: nthe, npsi, nzeta, nAzimRam, nXRaw

  implicit none

  INTEGER :: INDEXPA(nthe,npsi,nzeta,NPA) ! Index that maps the pitch angle on each point of a field line to the equatorial (RAM) one
  REAL(DP) :: flux3DEQ(NS, npsi,nzeta,NE,NPA)
  INTEGER, PARAMETER :: NRAD = nXRaw+1
  INTEGER :: INCFD, j, k, L
  REAL(DP) :: bfMirror(NPA), bfInterm(NPA)
  REAL(DP) :: chiMirror(NPA)
  REAL(DP), DIMENSION(nthe) :: derivs2, dBdTheta, dThetadB, derivs4, &
                               distance, dyDummyDens

contains

!==============================================================================
SUBROUTINE computehI(flux_volume)
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

  use ModRamVariables, ONLY: Kp, F107, FLUX
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamStart
  use ModRamConst,     ONLY: RE, HMIN, ME
  use ModRamParams,    ONLY: electric, IsComponent, NameBoundMag
  use ModRamCouple,    ONLY: SwmfPot_II

  USE ModScbMain,      ONLY: rHour, HourDecimal, Day, HourDecShort, iElectric, &
                             prefixOut
  use ModScbVariables, ONLY: bf, chiVal, nZetaMidnight, x, y, z, bnormal, PhiIono, &
                             radGrid, angleGrid, h_Cart, I_Cart, h_Cart_interp, &
                             i_Cart_interp, bZEq_Cart, flux_vol_cart, bZ, &
                             hdens_Cart, radRaw, azimRaw, nThetaEquator

  use ModRamFunctions, ONLY: RamFileName

  use ModScbSpline, ONLY: Spline_2D_point, spline, splint
  use ModScbInterp, ONLY: Interpolation_natgrid_2D

  USE nr,        ONLY: qtrap, qromb, locate
  use nrtype,    ONLY: DP, pi_d
  USE ModIOUnit, ONLY: UNITTMP_

  IMPLICIT NONE

  INTEGER :: I, IC, iCount, ierr, IP1, j1, k1, NWK, ncdfId, iFix, iDomain
  INTEGER, DIMENSION(3) :: dimlens 
  INTEGER, PARAMETER :: LIMIT = 10000, INF = 1
  REAL(DP) :: LZ(NRAD), RAD(NRAD), CONE(NRAD+4), PHI(NT), MLT(NT), PA(NPA), PAbn(NPA),  &
       DMU(NPA), WMU(NPA), MU(NPA)
  REAL(DP) :: MUEQ
  REAL(DP), DIMENSION(2*nthe) :: WK
  REAL(DP) :: switch, start_time, end_time, yLowerLimit, chiHalf, bfHalf
  REAL(DP) :: r0(npsi, nzeta)
  REAL(DP) :: MUBOUN, RWU, DL1, DR, IR1, DPHI, CLC, dyDummy
  INTEGER, PARAMETER :: LENIW = LIMIT, LENW = 4*LENIW+1
  INTEGER, PARAMETER :: nXRaw_local = nXRaw+1, nYRaw = nAzimRAM
  REAL(DP) :: bZEqDiff_Cart(nXRaw_local,nYRaw), Epot_Cart(0:nXRaw_local,nYRaw), beqdip_Cart(nXRaw_local)
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
  REAL(DP), PARAMETER :: radOut = RadiusMax+0.25_dp ! 9.25_dp ! 6.75_dp ! 10.25_dp ! Outer domain radius

  INTEGER :: iSpecies
  REAL(DP) :: EKEV2(NE), MU2(NPA)
  REAL(DP), PARAMETER :: b0dip = 30570._dp
  REAL(DP) :: beqdip(npsi,nzeta+1)

  REAL(DP), INTENT(IN) :: flux_volume(:,:)

  switch = 0.0_dp
  dimlens = (/1, 1, 1/)

  EPSABS = 0.01_dp
  EPSREL = 1.E-2_dp  ! 1% error should be fine for most applications

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

  DO j1 = 0, nXRaw_local
     radRaw(j1) = 1.75_dp + (radOut-1.75_dp) * REAL(j1,DP)/REAL(nXRaw_local,DP)
  END DO
  DO k1 = 1, nYRaw
     azimRaw(k1) = 24._dp * REAL(k1-1,DP)/REAL(nYRaw-1,DP)
  END DO

  !       NPA=no. of grids in equatorial pitch angle
  !       NPA=72  loss cone at 200 km
  DL1 = (radOut-2.25_dp)/REAL(NRAD-2, dp)     ! Grid size of L shell
  IR1=DL1/0.25_dp
  DR=DL1*RE               ! Grid size for RAD = RO
  DO I=1,NRAD
     LZ(I)=2.+(I-2)*DL1
     RAD(I)=RE*LZ(I)
  END DO

  DPHI=2._dp*PI_D/(NT-1)      ! Grid size for local time [rad]
  DO J=1,NT
     PHI(J)=(J-1)*DPHI     ! Magnetic local time in radian
     MLT(J)=PHI(J)*12./PI_D  ! Magnetic local time in hour
  END DO
  IP1=(MLT(2)-MLT(1))/0.5_dp

  !.......CONE - pitch angle loss cone in degree
  DO I = 1,NRAD
     CLC = (RE+HMIN)/RAD(I)
     CONE(I) = 180._dp/pi_d*ASIN(SQRT(CLC**3/SQRT(4.-3.*CLC)))
  END DO
  CONE(NRAD+1)=2.5_dp          ! to calculate PA grid near 0 deg
  CONE(NRAD+2)=1.5_dp
  CONE(NRAD+3)=1._dp
  CONE(NRAD+4)=0._dp


  !.......PA is equatorial pitch angle in deg - PA(1)=90, PA(NPA)=0.
  !       MU is cosine of equatorial PA
  PA(1) = 90._dp
  MU(1) = 0._dp
  PA(NPA) = 0._dp
  MU(NPA) = 1._dp
  RWU = 0.98_dp
  WMU(1) = (MU(NPA)-MU(1))/32._dp

  !                                    ! |_._|___.___|____.____|______.______|
  DO L = 1, 46                   !    MU    <  DMU   >    <     WMU     >
     WMU(L+1) = WMU(L)*RWU
     DMU(L) = 0.5_dp*(WMU(L)+WMU(L+1))
     MU(L+1) = MU(L)+DMU(L)
     PA(L+1) = 180._dp/pi_d*ACOS(MU(L+1))
  END DO

  PA(48) = 18.7_dp
  MU(48) = COS(pi_d/180._dp*PA(48))
  DMU(47) = (MU(48)-MU(47))
  IC = 2
  DO L = 48, NPA-1
     PA(L+1) = CONE(IC)
     IF(L == 49) THEN
        PA(50) = 16._dp
     ELSE
        IC = IC+1
     ENDIF
     MU(L+1) = COS(pi_d/180._dp*PA(L+1))
     DMU(L) = (MU(L+1)-MU(L))       ! Grid size in cos pitch angle
     WMU(L) = 2._dp*(DMU(L-1) - 0.5*WMU(L-1))
     IF(L > 55)WMU(L) = 0.5_dp*(DMU(L)+DMU(L-1))
  END DO
  DMU(NPA) = DMU(NPA-1)
  WMU(NPA) = DMU(NPA-1)
  DO L = 1, NPA-1
     MUBOUN = MU(L)+0.5_dp*WMU(L)
     PAbn(L) = 180._dp/pi_d*ACOS(MUBOUN)          ! PA at boundary of grid
  ENDDO
  PAbn(NPA)=0._dp

  INCFD = 1
  NWK = 2*(nthe-nThetaEquator+1)

  ! Always fill this matrix; it's used by RAM outputs
  ! Interpolate RAM flux on 3DEQ grid, for mapping
  DO iSpecies = 1, 4 ! ions & electrons
     DO L = 1, NPA
        DO K = 1, NE
           !Cubic spline interpolation
           CALL Spline_2D_point(radRaw(1:nXRaw), azimRaw*pi_d/12.0, &
                                REAL(FLUX(iSpecies, 1:nXRaw,1:nYRaw,K,L),DP), &   
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
           CALL spline(chiVal(nThetaEquator:nthe),distance(nThetaEquator:nthe), 1.E31_dp, 1.E31_dp, &
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
           CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe),  &
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
                 chiMirror(L) = chiVal(nthe)
              ELSE
                 ! Find chiMirror, assuming monotonicity of bf with chiVal
                 CALL DPCHFE (nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
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
              CALL DPCHIC (IC0, vc, switch, nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
                   dBdTheta(nThetaEquator:nthe), INCFD, WK(1:NWK), NWK, IERR)

              IF (ierr < 0) THEN
                 PRINT*, 'computehI2: ierr = ', ierr
                 STOP
              END IF

              valueIntegral = qromb(fInt, chiVal(nThetaEquator), chiMirror(L))  ! The integral for I is non-singular, qromb should do fine

              ! Split h integral or compute solely by DQAGI (which is pretty good for the whole interval)
              ! March 2006 - transformed h integral by y = 1/sqrt(bm-b) change of variable; spreads out the interval from yeq to infinity; not singular anymore 
              ! ---> better accuracy

              chiHalf = 0.2_dp*(chiVal(nThetaEquator) + 4.*chiMirror(L))
              CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
                   dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chiHalf, bfHalf, IERR)          

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

              CALL DQAG(fScalarInt, chiVal(nThetaEquator), chiHalf, EPSABS, EPSREL, 3, valueIntegral3, &
                   ABSERR, NEVAL, IERR, LIMIT, LENW, LAST, IWORK, WORK)   ! from eq. plane to chiHalf; regular integral (for h calculation)

              CALL DQAG(fScalarIntDens, chiVal(nThetaEquator), chiHalf, EPSABS, EPSREL, 3, valueIntegralDens3, &
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

     IF (electric=='IESC'.OR.electric=='WESC'.OR.electric=='W5SC') THEN ! Need potential mapping along SC B-field lines
        CALL Interpolation_natgrid_2D(radRaw(1:nXRaw_local), azimRaw, h_value(1:npsi,2:nzeta,1:NPA), &
                I_value(1:npsi,2:nzeta,1:NPA), &
                bZ(nThetaEquator,1:npsi,2:nzeta)*bnormal-beqdip(1:npsi,2:nzeta), &
                flux_volume(1:npsi,2:nzeta)/bnormal, hdens_value(1:npsi,2:nzeta,1:NPA), &
                h_Cart(1:nXRaw_local,1:nYRaw,1:NPA), I_Cart(1:nXRaw_local,1:nYRaw,1:NPA), &
                bZEqDiff_Cart(1:nXRaw_local,1:nYRaw), flux_vol_Cart(1:nXRaw_local,1:nYRaw), &
                hdens_Cart(1:nXRaw_local,1:nYRaw,1:NPA), PhiIono(1:npsi,2:nzeta), &
                Epot_Cart(1:nXRaw_local,1:nYRaw))
        Epot_Cart(0,:) = Epot_Cart(1,:) ! 3Dcode domain only extends to 2 RE; at 1.75 RE the potential is very small anyway
     ELSE
        CALL Interpolation_natgrid_2D(radRaw(1:nXRaw_local), azimRaw, h_value(1:npsi,2:nzeta,1:NPA), &
                I_value(1:npsi,2:nzeta,1:NPA), bZ(nThetaEquator,1:npsi,2:nzeta)*bnormal-beqdip(1:npsi,2:nzeta), &
                flux_volume(1:npsi,2:nzeta)/bnormal, hdens_value(1:npsi,2:nzeta,1:NPA), &
                h_Cart(1:nXRaw_local,1:nYRaw,1:NPA), I_Cart(1:nXRaw_local,1:nYRaw,1:NPA), &
                bZEqDiff_Cart(1:nXRaw_local,1:nYRaw), flux_vol_Cart(1:nXRaw_local,1:nYRaw), &
                hdens_Cart(1:nXRaw_local,1:nYRaw,1:NPA))
     END IF

     ! Add Bdip 
     beqdip_Cart(1:nXRaw_local) = b0dip/radRaw(1:nXRaw_local)**3

     ! If 3D code domain does not extend to radOut: (ACHTUNG!!!)
     ! print*, 'computehI: max(x), radOut = ', real(maxval(abs(x)),sp), real(radOut,SP)
     IF (MAXVAL(ABS(x)) <= radOut) THEN
        h_Cart(nXRaw_local,:,:) = h_Cart(nXRaw_local-1,:,:)
        I_Cart(nXRaw_local,:,:) = I_Cart(nXRaw_local-1,:,:)
        bZEqDiff_Cart(nXRaw_local,:) = bZEqDiff_Cart(nXRaw_local-1,:)
        hdens_Cart(nXRaw_local,:,:) = hdens_Cart(nXRaw_local-1,:,:)
        IF (iElectric == 1) Epot_Cart(nXRaw_local,:) = Epot_Cart(nXRaw_local-1,:)
     END IF

     ! Add dipole field
     DO j = 1, nYRaw
        bZEq_Cart(:,j)  = bZEqDiff_Cart(:,j) + beqdip_Cart
     END DO

     ! Make field quasi-dipole as R -> 2 RE
     DO i = 1, nXRaw_local
        IF (radRaw(i) <= 3.0) bZEq_Cart(i,:) = EXP(-2.*(radRaw(i)-3.0)**2)*bZEq_Cart(i,:) + &
             (1.-EXP(-2.*(RadRaw(i)-3.0)**2))*beqdip_Cart(i)
     END DO

     IF (MINVAL(h_Cart) < 0._dp .OR. MINVAL(I_Cart)<0._dp) THEN
        PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart)
        PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart)
        !   STOP
     END IF

!     Master_rank:  IF (rank == 0) THEN
        ! Cubic spline interpolation with natural boundaries to get h and I at muboun 
        DO j = 1, nYRaw
           DO i = 1, nXRaw_local
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

        ! Writes output file for RAM (hI...dat)
        ! Use correct naming convention.
     filenamehI=trim(prefixOut)//RamFileName('hI_output', 'dat', TimeRamNow)
     PRINT*, 'Writing to file ', TRIM(filenamehI); CALL flush(6)
     OPEN(UNITTMP_, file = TRIM(filenamehI), action = 'write', status = 'unknown')
!     WRITE(UNITTMP_, 15) 'T= ', rHour
     WRITE(UNITTMP_, 25) '   Lsh        MLT    NPA    h(Lsh,MLT,NPA)  hBoun(Lsh,MLT,NPA) '//&
             'I(Lsh,MLT,NPA) IBoun(Lsh,MLT,NPA) Bz(Lsh,MLT) HDENS(Lsh,MLT,NPA) fluxVolume(Lsh,MLT)'
     DO i = 1, nXRaw_local
        DO j = 1, nYRaw
           DO L = 1, NPA
              WRITE(UNITTMP_, *) radRaw(i), azimRaw(j), L, h_Cart(i,j,L), h_Cart_interp(i,j,L), &
                   I_Cart(i,j,L), I_Cart_interp(i,j,L), bZEq_Cart(i,j), hdens_Cart(i,j,L), flux_vol_Cart(i,j)
           END DO
        END DO
     END DO
     CLOSE(UNITTMP_)

        SWMF_electric_potential:  IF (electric=='IESC') THEN
           NameFileOut=trim(prefixOut)//trim(RamFileName('IE_SCB','in',TimeRamNow))
           PRINT*, 'Writing to file ', NameFileOut
           OPEN(UNITTMP_, file = NameFileOut, &
                action = 'write', status = 'unknown')
           ! Write header to file.  
           WRITE(UNITTMP_, '(A, i4.4)') ' DOM   UT      Kp   F107   Ap   Rs  PHIOFS, Year ', TimeRamNow%iYear
           WRITE(UNITTMP_, '(f5.0, 2x, f6.3, 2x, f4.2, 2x, f5.1, 2x, 2(f4.1,2x), f3.1)') &
                REAL(TimeRamNow%iDay), REAL(TimeRamNow%iHour) + &
                REAL(TimeRamNow%iMinute)/60.0, kp, f107, 0.0, 0.0, 0.0
           WRITE(UNITTMP_, '(A)') 'SWMF ionospheric potential mapped along SCB lines'
           WRITE(UNITTMP_, '(A)') 'L     PHI       Epot(kV)'
           DO i = 0, nXRaw_local
              DO j = 1, nYRaw
                 WRITE(UNITTMP_, 22) radRaw(i), azimRaw(j)*2.*pi_d/24., Epot_Cart(i,j)
              END DO
           END DO
           CLOSE(UNITTMP_)

           ! Save traced potential to ModRamCouple::SwmfPot_II
           IF(IsComponent) SwmfPot_II(1:nXRaw_local+1,1:nYRaw) = &
                Epot_Cart(0:nXRaw_local,  1:nYRaw)
        END IF SWMF_electric_potential

        Weimer_electric_potential_along_SCB: IF (electric=='WESC' .or. electric=='W5SC') THEN
           SwmfPot_II(1:nXRaw_local+1,1:nYRaw) = Epot_Cart(0:nXRaw_local, 1:nYRaw)
           print*, 'computehI: SwmfPot_II here SwmfPot_II mm', maxval(swmfpot_II), minval(swmfpot_ii)
        END IF Weimer_electric_potential_along_SCB
        
!     END IF Master_rank

15   FORMAT(A, F6.1)
21   FORMAT(F8.2, F10.1, 3X, I2, 5X, 8(3X, E12.4))
22   FORMAT(F4.2, 2X, F8.6, 2X, E11.4)
25   FORMAT(A)

  END SELECT Operational_or_research

  RETURN

END SUBROUTINE computehI

FUNCTION fInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local(:) 
  REAL(DP) :: bf_local(SIZE(chi_local))
  REAL(DP) :: dyDummy(SIZE(chi_local))
  INTEGER :: i, ierr
  REAL(DP) :: fInt(SIZE(chi_local))
  LOGICAL :: SKIP = .FALSE.

  CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
       dBdTheta(nThetaEquator:nthe), INCFD, SKIP, SIZE(chi_local,1), chi_local, bf_local, IERR)

  !C fInt = SQRT(MAX(1._dp/(1._dp - bf_local/bfMirror(L)), 0._dp)) ! For function h(mu0)
  fInt = SQRT(MAX((bfMirror(L) - bf_local)/bfMirror(L), 0._dp)) ! For function I(mu0)

  RETURN
END FUNCTION fInt


FUNCTION fIntDer(t_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
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


FUNCTION fScalarIntDens(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: x, y, z, bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: dyDummy(nthe), radius
  REAL(DP) :: rad(nthe)
  REAL(DP), PARAMETER :: HUGE = 1.E12_dp
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntDens
  LOGICAL :: SKIP = .FALSE.

  CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
       dBdTheta(nThetaEquator:nthe), INCFD, SKIP, 1, chi_local, bf_local, IERR)
  radius = splint(chiVal(nThetaEquator:nthe), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)
  fScalarIntDens = hdens_rairden(radius) * SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) 

  RETURN

END FUNCTION fScalarIntDens


FUNCTION fScalarInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: dyDummy
  REAL(DP), PARAMETER :: HUGE = 1.E12_dp
  INTEGER :: i, ierr
  REAL(DP) :: fScalarInt
  LOGICAL :: SKIP = .FALSE.

  CALL DPCHFE (nthe-nThetaEquator+1, chiVal(nThetaEquator:nthe), bf(nThetaEquator:nthe,j,k), &
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


FUNCTION fScalarIntInf(y_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: y_local
  REAL(DP) :: bf_local, chi_local, dChidB_local
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntInf
  LOGICAL :: SKIP = .FALSE.

  bf_local = bfMirror(L) - 1._dp/(y_local*y_local)

  CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
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


FUNCTION fScalarIntInfDens(y_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  USE ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: y_local
  REAL(DP) :: bf_local, chi_local, dChidB_local, radius
  INTEGER :: i, ierr
  REAL(DP) :: fScalarIntInfDens
  LOGICAL :: SKIP = .FALSE.

  bf_local = bfMirror(L) - 1._dp/(y_local*y_local)

  CALL DPCHFD(nthe-nThetaEquator+1,  bf(nThetaEquator:nthe,j,k), chiVal(nThetaEquator:nthe), &
       dThetadB(nThetaEquator:nthe), INCFD, SKIP, 1, bf_local, chi_local, dChidB_local, IERR)

  IF (IERR < 0) THEN  
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  radius = splint(chiVal(nThetaEquator:nthe), distance(nThetaEquator:nthe), dyDummyDens(nThetaEquator:nthe), chi_local)

  IF(y_local*y_local>0._dp) THEN
     fScalarIntInfDens = dChidB_local/(y_local*y_local) * hdens_rairden(radius)
  ELSE
     PRINT*, 'Problem in fScalarIntInf; y_local = ', y_local
  END IF

  RETURN
END FUNCTION fScalarIntInfDens


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
