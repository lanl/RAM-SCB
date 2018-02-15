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
  REAL(DP), allocatable :: chiMirror(:,:)
  REAL(DP), allocatable :: derivs2(:), dBdTheta(:), dThetadB(:), derivs4(:), &
                           distance(:), dyDummyDens(:)

  ! Variabels for testing stuff -ME
  INTEGER :: LMAX = 1000
  INTEGER :: IOPT, LOUT, nEquator
  REAL(DP) :: ER, DSMAX, RLIM, PARMOD(10), DIR, PS
  REAL(DP) :: x0, y0, z0, xe, ye, ze, xf, yf, zf, RIN
  REAL(DP), DIMENSION(1000) :: xx, yy, zz, cVal, dTemp, bTemp, dCdB, &
                               dBdC, dyDummyDensTemp, bxtemp, bytemp, bztemp

  integer :: time1, clock_rate = 1000, clock_max = 100000
  real(dp) :: starttime,stoptime

contains

!==============================================================================
subroutine ramscb_allocate

  use ModRamGrids, ONLY: NS, NE, NPA
  use ModScbGrids, ONLY: nthe, npsi, nzeta, nXRaw

  implicit none

  ALLOCATE(INDEXPA(nthe,npsi,nzeta,NPA), Flux3DEq(NS,npsi,nzeta,NE,NPA), &
           bfMirror(NPA), bfInterm(NPA), chiMirror(NPA,2), derivs2(nthe), &
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

  use ModScbIO,     ONLY: trace
  use ModScbSpline, ONLY: Spline_2D_point, spline, splint
  use ModScbInterp, ONLY: Interpolation_SCB_to_RAM

  USE nrmod,     ONLY: locate
  use nrtype,    ONLY: DP, pi_d
  USE ModIOUnit, ONLY: UNITTMP_

  IMPLICIT NONE

  INTEGER :: i, iFix, iCount, iDomain, iSpecies, k1, j1
  REAL(DP), PARAMETER :: b0dip = 30570._dp ! This doesn't really matter
  REAL(DP) :: DthI, radOut = RadiusMax+0.25_dp
  
  ! Variables for Hermite Splines
  INTEGER :: NWK, IERR
  REAL(DP) :: IC0(2), VC(2), switch
  LOGICAL :: SKIP = .false.

  ! Variables for QuadPack Integrals
  INTEGER, PARAMETER :: LIMIT = 10000, LENW = 4*LIMIT+1
  INTEGER :: IWORK(LIMIT), NEVAL, LAST
  REAL(DP) :: EPSABS, EPSREL, ABSERR, WORK(LENW), WK(LENW)

  ! Variables for SCB
  REAL(DP) :: dyDummy1(NPA), dyDummy2(NPA)
  REAL(DP) :: length, r0(npsi,nzeta), BeqDip(npsi,nzeta+1)
  REAL(DP) :: valueIntegralI, valueIntegralH, valueIntegralHDens
  REAL(DP), DIMENSION(npsi,nzeta,NPA) :: H_value, I_value, Hdens_value

  ! Variables for RAM
  REAL(DP) :: MUeq
  REAL(DP) :: BeqDip_Cart(nR), BzeqDiff_Cart(nR,nT)
  REAL(DP) :: BNESPrev(NR+1,NT), FNISPrev(NR+1,NT,NPA), BOUNISPrev(NR+1,NT,NPA)
  REAL(DP) :: scalingI(nT,nPA), scalingH(nT,nPA), scalingD(nT,nPA)

  switch = -1.0_dp

  EPSABS = 0.0001_dp
  EPSREL = 1.E-6_dp

  h_Cart = 0._dp
  I_Cart = 0._dp
  bZEq_Cart = 0._dp
  flux_vol_Cart = 0._dp
  bZEqDiff_Cart = 0._dp

  h_value = 0._dp
  hdens_value = 0._dp
  I_value = 0._dp
  valueIntegralI = 0._dp
  valueIntegralH = 0._dp
  valueIntegralHDens = 0._dp

  DO j1 = 0, nR
     radRaw(j1) = 1.75_dp + (radOut-1.75_dp) * REAL(j1,DP)/REAL(nR,DP)
  END DO
  DO k1 = 1, nT
     azimRaw(k1) = 24._dp * REAL(k1-1,DP)/REAL(nT-1,DP)
  END DO

  INCFD = 1
  NWK = 2*(nthe)+1

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

  Alpha_loop_parallel_oper:  DO k = 2, nzeta
     Psi_loop_oper: DO j = 1, npsi
!!!!! Is this really needed? -ME
        ! If bf(theta) not strictly increasing; to weed out very small differences  
        ! Generally if needed it means have to increase number of theta
        ! points (or crowd them more near the equatorial plane)

        !! Do nThetaEquator -> nthe
        iCount = 0
        Monotonicity_up: DO
           i = nThetaEquator
           iFix = 0
           iCount = iCount+1
           fixMonotonicity_up: DO WHILE (i < nthe)
              IF (bf(i+1,j,k) < bf(i,j,k)) THEN
                 bf(i+1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                 iFix = 1
              END IF
              i = i+1
           END DO fixMonotonicity_up
           IF (iFix == 0) EXIT Monotonicity_up
           CYCLE Monotonicity_up
        END DO Monotonicity_up

        !! Do nThetaEquator -> 1
        iCount = 0
        Monotonicity_down: DO
           i = nThetaEquator
           iFix = 0
           iCount = iCount+1
           fixMonotonicity_down: DO WHILE (i > 1)
              IF (bf(i-1,j,k) < bf(i,j,k)) THEN
                 bf(i-1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                 iFix = 1
              END IF
              i = i-1
           END DO fixMonotonicity_down
           IF (iFix == 0) EXIT Monotonicity_down
           CYCLE Monotonicity_down
        END DO Monotonicity_down
!!!!!

!!!!! Compute IndexPA to go along with Flux3D
        ! Define indexPA for 90 degree pitch angle (-1 for all distances
        ! along field line except equatorial plane)
        indexPA(:,:,:,1) = -1
        indexPA(nThetaEquator,:,:,1) = 1

        pitchAngleLoop_oper:        DO L = 2, NPA
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
!!!!!
        END DO pitchAngleLoop_oper
     END DO Psi_loop_oper
  END DO Alpha_loop_parallel_oper

  ! Periodicity for indexPA
  indexPA(:,:,1,:) = indexPA(:,:,nzeta,:)


  Operational_or_research: SELECT CASE (NameBoundMag)

  CASE('DIPL') ! Dipole without SCB calculation
     ! Nothing needed for a pure dipole
     RETURN 

  CASE default  ! Regular calculation of h, I, bounce-averaged charge xchange etc.
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

     RhoLoop: DO j = 2, npsi-1
        PhiLoop: DO k = 2,nzeta
           bfmirror = bf(nThetaEquator,j,k)/(1._dp - mu**2)

           r0(j, k) = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
           DO i = 1, nthe
              distance(i) = SQRT(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2) ! Distance from center of Earth
           END DO
           ! Compute length of field line (j,k)
           length = 0._dp
           DO i = 1, nthe-1
              length = length + SQRT((x(i+1,j,k)-x(i,j,k))**2 + (y(i+1,j,k)-y(i,j,k))**2 + (z(i+1,j,k)-z(i,j,k))**2)
           END DO

           ! Initialize Hermite Splines
           !! Boundary choices at the ends (look into DPCHIC etc. for explanation)
           IC0(1) = 0
           IC0(2) = -5
           VC(1) = 0
           VC(2) = 0
           CALL DPCHIC(IC0,vc,switch,nthe,chiVal(:),bf(:,j,k),dBdTheta(:),INCFD,WK(1:NWK),NWK,IERR)
           ! Need to split this call into two parts so bf is monotonically increasing
           !! nThetaEquator -> nthe
           CALL DPCHIC(IC0,vc,switch,nthe-nThetaEquator+1,bf(nThetaEquator:nthe,j,k), &
                       chiVal(nThetaEquator:nthe),dThetadB(nThetaEquator:nthe),INCFD,WK(1:NWK),NWK,IERR)
           !! nThetaEquator -> 1
           CALL DPCHIC(IC0,vc,switch,nThetaEquator,bf(nThetaEquator:1:-1,j,k), &
                       chiVal(nThetaEquator:1:-1),dThetadB(nThetaEquator:1:-1),INCFD,WK(1:NWK),NWK,IERR)
           CALL spline(chiVal(:),distance(:), 1.E31_dp, 1.E31_dp, dyDummyDens(:))

           PitchAngleLoop: DO L = 2,NPA
              ! Find ChiMirror points
              IF (bfMirror(L).gt.bf(nthe,j,k)) THEN
                 bfmirror(L) = bf(nthe,j,k)
                 chiMirror(L,1) = chiVal(nthe)
                 chiMirror(L,2) = chiVal(1)
              ELSE
                 !! Find upper mirror point
                 CALL DPCHFE(nthe-nThetaEquator+1,bf(nThetaEquator:nthe,j,k),chiVal(nThetaEquator:nthe), &
                             dThetadB(nThetaEquator:nthe),INCFD,SKIP,1,bfMirror(L),chiMirror(L,1),IERR)
                 !! Find lower mirror point
                 CALL DPCHFE(nThetaEquator,bf(nThetaEquator:1:-1,j,k),chiVal(nThetaEquator:1:-1), &
                             dThetadB(nThetaEquator:1:-1),INCFD,SKIP,1,bfMirror(L),chiMirror(L,2),IERR)
              ENDIF
              !! Fix any issues with the mirror points
              if (chiMirror(L,2).lt.0.) chiMirror(L,2) = 0._dp ! Occasionally get a negative number when near 0
              if (chiMirror(L,1).lt.chiMirror(L,2)) then ! Occasionally messes up for low pitch angles
                 chiMirror(L,1) = pi_d - chiMirror(L,2)
                 if (chiMirror(L,1).lt.0.) chiMirror(L,1) = chiMirror(L,2) + 0.001
              endif

              ! Evaluate I Integral, can use DQAG since I has no singularities
              CALL DQAG(fIInt,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,3,valueIntegralI, &
                        ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

              ! Evaluate H Integral, need to use DQAGS since it has end point singularities
              CALL DQAGS(fHInt,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,valueIntegralH, &
                         ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

              ! Evaluate HDens Integral
              CALL DQAGS(fHDensInt,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,valueIntegralHDens, &
                         ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

              ! Populate SCB arrays
              I_value(j,k,L) = (length/(pi_d*r0(j,k))) * valueIntegralI/SQRT(Bfmirror(L)) ! Don't understand where the length and pi come from -ME
              H_value(j,k,L) = (length/(pi_d*2*r0(j,k))) * valueIntegralH*SQRT(Bfmirror(L))
              HDens_value(j,k,L) = 1.E10_dp * valueIntegralHDens/valueIntegralH ! Re-normalize (see definition in hdens_rairden)

              !if (I_value(j,k,L)+1.eq.I_value(j,k,L)) I_Value(j,k,L) = I_Value(j,k-1,L)
              !if (H_value(j,k,L)+1.eq.H_value(j,k,L)) H_Value(j,k,L) = H_Value(j,k-1,L)
              !if ((HDens_value(j,k,L)+1.eq.HDens_value(j,k,L)).or.(hDens_value(j,k,L).eq.0)) then
              !   Hdens_Value(j,k,L) = Hdens_Value(j,k-1,L)
              !endif
           END DO PitchAngleLoop
        END DO PhiLoop
     END DO RhoLoop
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       write(*,*) 'SCB h and I: ',stoptime-starttime

     ! Handle pitch angle of 90 degrees
     I_value(:,:,1) = 0._dp
     H_value(:,:,2) = H_value(:,:,3) ! These two should be changed to using an actual calculation
     H_value(:,:,1) = H_value(:,:,2) ! for near 90 degree pitch angle particles -ME
     HDens_value(:,:,1) = HDens_value(:,:,2)

     ! Handle periodicity
     I_value(:,1,:) = I_value(:,nzeta,:)
     H_value(:,1,:) = H_value(:,nzeta,:)

     ! Now interpolate from SCB grid to RAM grid
     !! Interpolation will be for B - Bdip
     beqdip(1:npsi,2:nzeta)=b0dip/(x(nThetaEquator,1:npsi,2:nzeta)**2+y(nThetaEquator,1:npsi,2:nzeta)**2)**1.5

     !! Interpolate data for output in POLAR coordinates (for RAM)
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

     CALL Interpolation_SCB_to_RAM(radRaw(1:nR), azimRaw, &
                                   h_value(2:npsi-1,2:nzeta,1:NPA), &
                                   I_value(2:npsi-1,2:nzeta,1:NPA), &
                                   bZ(nThetaEquator,2:npsi-1,2:nzeta)*bnormal-beqdip(2:npsi-1,2:nzeta), &
                                   fluxVolume(2:npsi-1,2:nzeta)/bnormal, &
                                   hdens_value(2:npsi-1,2:nzeta,1:NPA), &
                                   h_Cart(1:nR,1:nT,1:NPA), &
                                   I_Cart(1:nR,1:nT,1:NPA), &
                                   bZEqDiff_Cart(1:nR,1:nT), &
                                   flux_vol_Cart(1:nR,1:nT), &
                                   hdens_Cart(1:nR,1:nT,1:NPA))

       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       write(*,*) 'Interpolation: ',stoptime-starttime

     !! Add dipole field back
     beqdip_Cart(1:nR) = b0dip/radRaw(1:nR)**3
     DO j = 1, nT
        bZEq_Cart(:,j)  = bZEqDiff_Cart(:,j) + beqdip_Cart
     END DO

     !! Check for points from the RAM grid that were outside the SCB grid
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

     DIR = -1.0
     DSMAX = 0.1
     ER = 0.001
     RLIM = 20.0
     PS = 0._dp
     IOPT = floor(kp)
     PARMOD = (/1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
     scalingI = 0._dp
     scalingH = 0._dp
     scalingD = 0._dp
     do i = 1,nR
        do j = 1,nT
           if (h_Cart(i,j,1).eq.-1) then
              ! Sample code for now
              !! Get equatorial point
              x0 = LZ(i) * COS(MLT(j)*2._dp*pi_d/24._dp - pi_d)
              y0 = LZ(i) * SIN(MLT(j)*2._dp*pi_d/24._dp - pi_d)
              !! Trace from equatorial point to pole then pole to other pole
              CALL trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD, &
                         xe,ye,ze,xx(:),yy(:),zz(:),LOUT,LMAX,bxtemp,bytemp,bztemp)
              CALL trace(xe,ye,ze,-DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD, &
                         xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bxtemp,bytemp,bztemp)
              !! Get necessary at each point along the trace
              length = 0._dp
              nEquator = 0
              cVal(1) = 0._dp
              do k = 1,LOUT
                 dTemp(k) = SQRT(xx(k)**2+yy(k)**2+zz(k)**2) ! Distance from center of earth
                 bTemp(k) = SQRT(bxtemp(k)**2+bytemp(k)**2+bztemp(k)**2)/bnormal
                 if (k.gt.1) then
                    cVal(k) = cVal(k-1) + SQRT((xx(k)-xx(k-1))**2 + (yy(k)-yy(k-1))**2 + (zz(k)-zz(k-1))**2)
                    length = length + SQRT((xx(k)-xx(k-1))**2 + (yy(k)-yy(k-1))**2 + (zz(k)-zz(k-1))**2) ! Length of field line
                    if ((btemp(k).gt.btemp(k-1)).and.(nEquator.eq.0)) nEquator = k-1
                 endif
              enddo
              bZEq_Cart(i,j) = bztemp(nEquator)
              cVal = pi_d*cVal/cVal(LOUT)
              r0(i,j) = SQRT(xx(nEquator)**2+yy(nEquator)**2+zz(nEquator)**2)
              !! Now compute stuff (copied from above)
              bfmirror = btemp(nEquator)/(1._dp - mu**2)
              ! Initialize Hermite Splines
              !! Boundary choices at the ends (look into DPCHIC etc. for explanation)
              IC0(1) = 0
              IC0(2) = -5
              VC(1) = 0
              VC(2) = 0
              NWK = 2*(LOUT)+1
              CALL DPCHIC(IC0,vc,switch,LOUT,cVal(1:LOUT),btemp(1:LOUT),dBdC(1:LOUT),INCFD,WK(1:NWK),NWK,IERR)
              ! Need to split this call into two parts so bf is monotonically increasing
              !! nThetaEquator -> nthe
              CALL DPCHIC(IC0,vc,switch,LOUT-nEquator+1,btemp(nEquator:LOUT), &
                          cVal(nEquator:LOUT),dCdB(nEquator:LOUT),INCFD,WK(1:NWK),NWK,IERR)
              !! nThetaEquator -> 1
              CALL DPCHIC(IC0,vc,switch,nEquator,btemp(nEquator:1:-1), &
                          cVal(nEquator:1:-1),dCdB(nEquator:1:-1),INCFD,WK(1:NWK),NWK,IERR)
              CALL spline(cVal(1:LOUT), dTemp(1:LOUT), 1.E31_dp, 1.E31_dp, dyDummyDensTemp(1:LOUT))

              DO L = 2,NPA
                 ! Find ChiMirror points
                 IF (bfMirror(L).gt.btemp(LOUT)) THEN
                    bfmirror(L) = btemp(LOUT)
                    chiMirror(L,1) = cVal(LOUT)
                    chiMirror(L,2) = cVal(1)
                 ELSE
                    !! Find upper mirror point
                    CALL DPCHFE(LOUT-nEquator+1,btemp(nEquator:LOUT),cVal(nEquator:LOUT), &
                                dCdB(nEquator:LOUT),INCFD,SKIP,1,bfMirror(L),chiMirror(L,1),IERR)
                    !! Find lower mirror point
                    CALL DPCHFE(nEquator,btemp(nEquator:1:-1),cVal(nEquator:1:-1), &
                                dCdB(nEquator:1:-1),INCFD,SKIP,1,bfMirror(L),chiMirror(L,2),IERR)
                 ENDIF
                 !! Fix any issues with the mirror points
                 if (chiMirror(L,2).lt.0.) chiMirror(L,2) = 0._dp ! Occasionally get a negative number when near 0
                 if (chiMirror(L,1).lt.chiMirror(L,2)) then ! Occasionally messes up for low pitch angles
                    chiMirror(L,1) = pi_d - chiMirror(L,2)
                    if (chiMirror(L,1).lt.0.) chiMirror(L,1) = chiMirror(L,2) + 0.001
                 endif

                 ! Evaluate I Integral, can use DQAG since I has no singularities
                 CALL DQAG(fIIntSub,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,3,valueIntegralI, &
                           ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

                 ! Evaluate H Integral, need to use DQAGS since it has end point singularities
                 CALL DQAGS(fHIntSub,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,valueIntegralH, &
                            ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

                 ! Evaluate HDens Integral
                 CALL DQAGS(fHDensIntSub,chiMirror(L,2),chiMirror(L,1),EPSABS,EPSREL,valueIntegralHDens, &
                            ABSERR,NEVAL,IERR,LIMIT,LENW,LAST,IWORK,WORK)

                 I_cart(i,j,L) = (length/(pi_d*r0(i,j))) * valueIntegralI/SQRT(Bfmirror(L))
                 H_cart(i,j,L) = (length/(pi_d*2*r0(i,j))) * valueIntegralH*SQRT(Bfmirror(L))
                 HDens_cart(i,j,L) = 1.E10_dp * valueIntegralHDens/valueIntegralH

                 if (scalingI(j,L).eq.0) scalingI(j,L) = I_cart(i-1,j,L)/I_cart(i,j,L)
                 if (scalingH(j,L).eq.0) scalingH(j,L) = H_cart(i-1,j,L)/H_cart(i,j,L)
                 if (scalingD(j,L).eq.0) scalingD(j,L) = HDens_cart(i-1,j,L)/HDens_cart(i,j,L)

                 I_cart(i,j,L) = I_cart(i,j,L)*scalingI(j,L)
                 H_cart(i,j,L) = H_cart(i,j,L)*scalingH(j,L)
                 HDens_cart(i,j,L) = HDens_cart(i,j,L)*scalingD(j,L)
              ENDDO
           endif
        enddo
     enddo
     I_cart(:,:,1) = 0._dp
     H_cart(:,:,2) = H_cart(:,:,3) ! Again, need to add actual 90 degree calculation -ME
     H_cart(:,:,1) = H_cart(:,:,2)
     HDens_cart(:,:,1) = HDens_cart(:,:,2)

       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       write(*,*) 'RAM h and I: ',stoptime-starttime

     ! Make sure all H and I integrals have been filled in correctly
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
           DO L = 2, NPA-1
              h_Cart_interp(i,j,L) = splint(PA(1:NPA),h_Cart(i,j,1:NPA),dyDummy1(1:NPA),PAbn(L))
              I_Cart_interp(i,j,L) = splint(PA(1:NPA),I_Cart(i,j,1:NPA),dyDummy2(1:NPA),PAbn(L))
           END DO
           ! Do not do the NPA inclusive in the not-a-knot interpolation above -> can lead to negative h,I(NPA)
           h_Cart_interp(i,j,NPA) = h_Cart_interp(i,j,NPA-1)
           I_Cart_interp(i,j,NPA) = I_Cart_interp(i,j,NPA-1)
           h_Cart_interp(i,j,1) = h_Cart_interp(i,j,2)
           I_Cart_interp(i,j,1) = I_Cart_interp(i,j,2)
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

15   FORMAT(A, F6.1)
21   FORMAT(F8.2, F10.1, 3X, I2, 5X, 8(3X, E12.4))

  END SELECT Operational_or_research

  RETURN

END SUBROUTINE computehI

!==============================================================================
FUNCTION fIInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  INTEGER  :: ierr
  REAL(DP) :: fIInt
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE(nthe,chiVal(:),bf(:,j,k),dBdTheta(:),INCFD,SKIP,1,chi_local,bf_local,IERR)
  !C if (chi_local > chiVal(nthe)) STOP 'Problem in fScalarInt.'

  IF (IERR < 0) THEN  ! IERR = 0 normal return, IERR > 0 non-fatal error
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  fIInt = SQRT(MAX((bfMirror(L) - bf_local), 0._dp)) ! For function I(mu0)

  RETURN

END FUNCTION fIInt

!==============================================================================
FUNCTION fIIntSub(chi_local)
  USE nrtype,          ONLY: DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  INTEGER  :: ierr
  REAL(DP) :: fIIntSub
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE(LOUT,cVal(1:LOUT),btemp(1:LOUT),dBdC(1:LOUT),INCFD,SKIP,1,chi_local,bf_local,IERR)
  !C if (chi_local > chiVal(nthe)) STOP 'Problem in fScalarInt.'

  IF (IERR < 0) THEN  ! IERR = 0 normal return, IERR > 0 non-fatal error
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  fIIntSub = SQRT(MAX((bfMirror(L) - bf_local), 0._dp)) ! For function I(mu0)

  RETURN

END FUNCTION fIIntSub

!==============================================================================
FUNCTION fHInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbVariables, ONLY: bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  INTEGER  :: ierr
  REAL(DP) :: fHInt
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE(nthe,chiVal(:),bf(:,j,k),dBdTheta(:),INCFD,SKIP,1,chi_local,bf_local,IERR)
  !C if (chi_local > chiVal(nthe)) STOP 'Problem in fScalarInt.'

  IF (IERR < 0) THEN  ! IERR = 0 normal return, IERR > 0 non-fatal error
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  fHInt = SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) ! For function h(mu0)

  RETURN
END FUNCTION fHInt

!==============================================================================
FUNCTION fHIntSub(chi_local)
  USE nrtype,          ONLY: DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  INTEGER  :: ierr
  REAL(DP) :: fHIntSub
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE(LOUT,cVal(1:LOUT),btemp(1:LOUT),dBdC(1:LOUT),INCFD,SKIP,1,chi_local,bf_local,IERR)
  !C if (chi_local > chiVal(nthe)) STOP 'Problem in fScalarInt.'

  IF (IERR < 0) THEN  ! IERR = 0 normal return, IERR > 0 non-fatal error
     PRINT*, 'IERR = ', IERR
     STOP
  END IF

  fHIntSub = SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) ! For function h(mu0)

  RETURN
END FUNCTION fHIntSub

!==============================================================================
FUNCTION fHDensInt(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbGrids,     ONLY: nthe
  use ModScbSpline,    ONLY: splint
  use ModScbVariables, ONLY: x, y, z, bf, chiVal, bnormal, nThetaEquator
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: radius
  INTEGER  :: ierr
  REAL(DP) :: fHDensInt
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE (nthe,chiVal(:),bf(:,j,k),dBdTheta(:),INCFD,SKIP,1,chi_local,bf_local,IERR)
  radius = splint(chiVal(:), distance(:), dyDummyDens(:), chi_local)

  fHDensInt = hdens_rairden(radius) * SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp)) 

  RETURN

END FUNCTION fHDensInt

!==============================================================================
FUNCTION fHDensIntSub(chi_local)
  USE nrtype,          ONLY: DP
  use ModScbSpline,    ONLY: splint
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: chi_local
  REAL(DP) :: bf_local
  REAL(DP) :: radius
  INTEGER  :: ierr
  REAL(DP) :: fHDensIntSub
  LOGICAL  :: SKIP = .FALSE.

  CALL DPCHFE(LOUT,cVal(1:LOUT),btemp(1:LOUT),dBdC(1:LOUT),INCFD,SKIP,1,chi_local,bf_local,IERR)
  radius = splint(cVal(1:LOUT), dTemp(1:LOUT), dyDummyDensTemp(1:LOUT), chi_local)

  fHDensIntSub = hdens_rairden(radius) * SQRT(MAX(1._dp/(bfMirror(L) - bf_local), 0._dp))

  RETURN

END FUNCTION fHDensIntSub

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
