!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModRamScb
! Contains subroutines and functions necessary to couple RAM and SCB

  use nrtype, ONLY: DP


  implicit none

  INTEGER, allocatable :: INDEXPA(:,:,:,:) ! Index that maps the pitch angle on each point of a field line to the equatorial (RAM) one
  REAL(DP), allocatable :: flux3DEQ(:,:,:,:,:)

contains

!==============================================================================
subroutine ramscb_allocate

  use ModRamGrids, ONLY: NS, NE, NPA
  use ModScbGrids, ONLY: nthe, npsi, nzeta, nXRaw


  implicit none

  ALLOCATE(INDEXPA(nthe,npsi,nzeta,NPA), Flux3DEq(NS,npsi,nzeta,NE,NPA))
  INDEXPA = 0.0; Flux3DEq = 0.0

  return

end subroutine ramscb_allocate

!==============================================================================
subroutine ramscb_deallocate


  implicit none

  DEALLOCATE(INDEXPA, Flux3DEq)

end subroutine ramscb_deallocate

!==============================================================================
subroutine Compute3DFlux

  use ModRamVariables, ONLY: FLUX, MU, FFACTOR, FNHS, F2
  use ModRamGrids,     ONLY: nR, nT, nPa, nE, radout

  USE ModScbParams,    ONLY: Symmetric
  use ModScbGrids,     ONLY: nthe, npsi, nzeta
  use ModScbVariables, ONLY: bf, radGrid, angleGrid, radRaw, azimRaw, &
                             nThetaEquator, x, y

  use ModRamGSL,       ONLY: GSL_Interpolation_2D
  use ModScbFunctions, ONLY: locate

  use nrtype, ONLY: DP, pi_d, twopi_d


  implicit none

  integer :: i, j, k, L, iS, GSLErr
  real(DP) :: MuEq, radius, angle

  DO iS = 1,4
     DO I = 2, NR
        DO K = 2, NE
           DO L = 2, NPA
              DO J = 1, NT-1
                 FLUX(iS,I,J,K,L) = F2(iS,I,J,K,L)/FFACTOR(iS,I,K,L)/FNHS(I,J,L)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO j = 0, nR
     radRaw(j) = 1.75_dp + (radOut-1.75_dp) * REAL(j,DP)/REAL(nR,DP)
  END DO
  DO k = 1, nT
     azimRaw(k) = 24._dp * REAL(k-1,DP)/REAL(nT-1,DP)
  END DO

  DO k = 2, nzeta
     DO j = 1, npsi
        radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
        angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
             angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
             angle = - ASIN(y(nThetaEquator,j,k) / radius)
        radGrid(j,k) = radius
        angleGrid(j,k) = angle
     END DO
  END DO

  ! Always fill this matrix; it's used by RAM outputs
  ! Interpolate RAM flux on 3DEQ grid, for mapping
  DO iS = 1, 4 ! ions & electrons
     DO L = 1, NPA
        DO K = 1, NE
           !Cubic GSL interpolation
           CALL GSL_Interpolation_2D(radRaw(1:nR), azimRaw*pi_d/12.0, &
                                     FLUX(iS,1:nR,1:nT,K,L), &
                                     radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), &
                                     flux3DEQ(iS,1:npsi,2:nzeta,K,L),GSLerr)
        END DO
     END DO
  END DO
  ! Periodicity
  flux3DEQ(:,:,1,:,:) = flux3DEQ(:,:,nzeta,:,:)

  DO k = 2, nzeta
     DO j = 1, npsi
        ! If bf(theta) not strictly increasing; to weed out very small
        ! differences  
        ! Generally if needed it means have to increase number of theta
        ! points (or crowd them more near the equatorial plane)

        !! Do nThetaEquator -> nthe
        Monotonicity_up: DO
           i = nThetaEquator
           L = 0
           fixMonotonicity_up: DO WHILE (i < nthe)
              IF (bf(i+1,j,k) < bf(i,j,k)) THEN
                 bf(i+1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                 L = 1
              END IF
              i = i+1
           END DO fixMonotonicity_up
           IF (L == 0) EXIT Monotonicity_up
           CYCLE Monotonicity_up
        END DO Monotonicity_up

        !! Do nThetaEquator -> 1
        Monotonicity_down: DO
           i = nThetaEquator
           L = 0
           fixMonotonicity_down: DO WHILE (i > 1)
              IF (bf(i-1,j,k) < bf(i,j,k)) THEN
                 bf(i-1,j,k) = bf(i,j,k)*(1._dp+1.E-15_dp)
                 L = 1
              END IF
              i = i-1
           END DO fixMonotonicity_down
           IF (L == 0) EXIT Monotonicity_down
           CYCLE Monotonicity_down
        END DO Monotonicity_down

!!!!! Compute IndexPA to go along with Flux3D
        ! Define indexPA for 90 degree pitch angle (-1 for all distances
        ! along field line except equatorial plane)
        indexPA(:,:,:,1) = -1
        indexPA(nThetaEquator,:,:,1) = 1

        DO L = 2, NPA
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
        END DO
     END DO
  END DO
  
  ! Periodicity for indexPA
  indexPA(:,:,1,:) = indexPA(:,:,nzeta,:)

  return
end subroutine Compute3DFlux

!==============================================================================
SUBROUTINE computehI(iter)
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
                             RLZ, LZ, MU, DMU, WMU, MLT, PAbn, PA, DL1
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamStart, DT_hI, TimeRamElapsed, TOld
  use ModRamConst,     ONLY: RE, HMIN, ME, b0dip
  use ModRamParams,    ONLY: electric, IsComponent, NameBoundMag
  use ModRamCouple,    ONLY: SwmfPot_II
  use ModRamGrids,     ONLY: nR, nT, nPa, nE, RadiusMax, radout

  USE ModScbMain,      ONLY: rHour, HourDecimal, Day, HourDecShort, iElectric, &
                             prefixOut
  USE ModScbParams,    ONLY: Symmetric
  use ModScbGrids,     ONLY: nthe, npsi, nzeta, nXRaw, nAzimRAM
  use ModScbVariables, ONLY: bf, chiVal, nZetaMidnight, x, y, z, bnormal, PhiIono, &
                             radGrid, angleGrid, h_Cart, I_Cart, h_Cart_interp, &
                             i_Cart_interp, bZEq_Cart, flux_vol_cart, bZ, chi, &
                             hdens_Cart, radRaw, azimRaw, nThetaEquator, fluxVolume, &
                             thetaVal

  use ModRamGSL,       ONLY: GSL_Interpolation_2D, GSL_Interpolation_1D, &
                             GSL_Integration_hI
  use ModRamFunctions, ONLY: RamFileName, COSD, FUNT, FUNI
  use ModScbFunctions, ONLY: extap, locate
  use ModScbIO,        ONLY: trace, PARMOD, IOPT

  use nrtype,    ONLY: DP, pi_d, pio2_d
  USE ModIOUnit, ONLY: UNITTMP_


  implicit none

  INTEGER, INTENT(IN) :: iter
  INTEGER :: ii, j, iS, GSLerr
  
  ! Variables for timing
  integer :: time1, clock_rate, clock_max
  real(dp) :: starttime,stoptime

  ! Variables for SCB
  integer  :: statusIntegralI, statusIntegralH, statusIntegralHDens
  REAL(DP) :: DthI
  REAL(DP), ALLOCATABLE :: length(:,:), r0(:,:), BeqDip(:,:)
  REAL(DP), ALLOCATABLE :: distance(:,:,:), H_value(:,:,:), I_value(:,:,:), &
                           Hdens_value(:,:,:), yI(:,:,:), yH(:,:,:), yD(:,:,:), &
                           bfmirror(:,:,:)

  ! Variables for RAM
  logical  :: Extrap
  integer  :: nEquator
  REAL(DP) :: MUeq, wn, xo, xn, xp, yo, yn, yp, t0, t1, rt, tt, zt
  REAL(DP) :: r0_RAM, length_RAM
  REAL(DP), ALLOCATABLE :: BeqDip_Cart(:), bfMirror_RAM(:), yI_RAM(:), yH_RAM(:), &
                           yD_RAM(:), bbx(:), bby(:), bbz(:), bb(:), dd(:), &
                           cVal(:), xx(:), yy(:), zz(:), ScaleAt(:)
  REAL(DP), ALLOCATABLE :: BzeqDiff_Cart(:,:), BNESPrev(:,:)
  REAL(DP), ALLOCATABLE :: FNISPrev(:,:,:), BOUNISPrev(:,:,:)
  REAL(DP) :: scalingI, scalingH, scalingD, I_Temp, H_Temp, D_Temp

  ! Variables for Tracing
  INTEGER :: LMAX, LOUT
  INTEGER :: ID
  REAL(DP) :: x0, y0, z0, xe, ye, ze, xf, yf, zf, RIN
  REAL(DP) :: ER, DSMAX, RLIM, DIR, PS
  REAL(DP) :: PDyn, BzIMF, DIST, XMGNP, YMGNP, ZMGNP
  REAL(DP), DIMENSION(200) :: xtemp, ytemp, ztemp, bxtemp, bytemp, bztemp, dtemp

  integer, save :: i, k, L
  !$OMP THREADPRIVATE(i, k, L)

  clock_rate = 1000
  clock_max = 100000
  LMAX = 200

  ALLOCATE(length(npsi,nzeta+1), r0(npsi,nzeta+1), BeqDip(npsi,nzeta+1))
  length = 0.0; r0 = 0.0; BEqDip = 0.0
  ALLOCATE(distance(nthe,npsi,nzeta))
  distance = 0.0
  ALLOCATe(H_value(npsi,nzeta,NPA), I_value(npsi,nzeta,NPA), HDens_value(npsi,nzeta,NPA), &
           yI(npsi,nzeta,NPA), yH(npsi,nzeta,NPA), yD(npsi,nzeta,NPA), bfmirror(npsi,nzeta,NPA))
  H_value = 0.0; I_value = 0.0; HDens_Value = 0.0; yI = 0.0; yH = 0.0; yD = 0.0
  bfMirror = 0.0
  ALLOCATE(BeqDip_Cart(nR))
  BEqDip_Cart = 0.0
  ALLOCATE(ScaleAt(nT))
  ScaleAt = 0
  ALLOCATE(bfMirror_RAM(nPa), yI_RAM(nPA), yH_RAM(nPA), yD_RAM(nPA))
  bfMirror_RAM = 0.0; yI_RAM = 0.0; yH_RAM = 0.0; yD_RAM = 0.0
  ALLOCATE(bbx(nthe), bby(nthe), bbz(nthe), bb(nthe), dd(nthe), cVal(nthe), &
           xx(nthe), yy(nthe), zz(nthe))
  bbx = 0.0; bby = 0.0; bbz = 0.0; bb = 0.0; dd = 0.0; cVal = 0.0
  xx = 0.0; yy = 0.0; zz = 0.0
  ALLOCATE(BzeqDiff_Cart(nR,nT))
  BzEqDiff_Cart = 0.0
  ALLOCATE(BNESPrev(nR+1,nT))
  BNESPrev = 0.0
  ALLOCATE(FNISPrev(nR+1,nT,nPa),BOUNISPrev(nR+1,nT,nPa))
  FNISPrev = 0.0; BOUNISPrev = 0.0

  h_Cart = 0._dp
  I_Cart = 0._dp
  bZEq_Cart = 0._dp
  flux_vol_Cart = 0._dp
  bZEqDiff_Cart = 0._dp

  h_value = 0._dp
  hdens_value = 0._dp
  I_value = 0._dp

  DO j = 0, nR
     radRaw(j) = 1.75_dp + (radOut-1.75_dp) * REAL(j,DP)/REAL(nR,DP)
  END DO
  DO k = 1, nT
     azimRaw(k) = 24._dp * REAL(k-1,DP)/REAL(nT-1,DP)
  END DO

  Operational_or_research: SELECT CASE (NameBoundMag)
  CASE('DIPL') ! Dipole without SCB calculation
     ! If this is the first iteration we need to populate the h and I integrals
     if (iter.ne.0) return

     ! Start RAM Timing
     call system_clock(time1,clock_rate,clock_max)
     starttime=time1/real(clock_rate,dp)
     ! Loop over RAM grid
     do i = 1,nR
        do j = 1,1 ! Dipole is azimuthially symmetric, only need to calculate one MLT
           ! Calculate h and I at RAM grid point
           zt = MLT(j)*2._dp*pi_d/24._dp - pi_d
           t1 = dasin(dsqrt(1.0/LZ(i+1)))
           t0 = pi_d-t1
           do k=1,nthe
              tt = t0 + REAL(k-1,DP)/REAL(nthe-1,DP)*(t1-t0) 
              rt = LZ(i+1)*dsin(tt)**2
              xx(k) = (rt)*dcos(zt)*dsin(tt)
              yy(k) = (rt)*dsin(zt)*dsin(tt)
              zz(k) = (rt)*dcos(tt)
              dd(k) = rt ! Distance from center of earth
              bb(k) = (b0dip/bnormal)*dsqrt(1+3*dcos(tt)**2)/rt**3
           enddo

           cVal(1) = 0._dp
           length_RAM = 0._dp
           nEquator = (nthe+1)/2
           do k = 2,nthe
              cVal(k) = cVal(k-1) + SQRT((xx(k)-xx(k-1))**2 &
                                       + (yy(k)-yy(k-1))**2 &
                                       + (zz(k)-zz(k-1))**2)
              length_RAM = length_RAM + SQRT((xx(k)-xx(k-1))**2 &
                                           + (yy(k)-yy(k-1))**2 &
                                           + (zz(k)-zz(k-1))**2)
            enddo
           bZEq_Cart(i,j) = bb(nEquator)*bnormal !Since in a dipole Beq = Bzeq
           cVal(1:nthe) = pi_d*cVal(1:nthe)/cVal(nthe)
           r0_RAM = SQRT(xx(nEquator)**2+yy(nEquator)**2)

           bfmirror_RAM(:) = bb(nEquator)/(1._dp - mu**2)
           bfmirror_RAM(NPA) = bb(nthe)

           CALL GSL_Integration_hI(bfMirror_RAM(:),cVal(1:nthe),bb(1:nthe), &
                                   dd(1:nthe), yI_RAM(:), yH_RAM(:), yD_RAM(:))

           I_cart(i,j,:) = (length_RAM/(pi_d*r0_RAM)) * yI_RAM(:)/SQRT(Bfmirror_RAM(:))
           H_cart(i,j,:) = (length_RAM/(pi_d*2*r0_RAM)) * yH_RAM(:)*SQRT(Bfmirror_RAM(:))
           HDens_cart(i,j,:) = 1.E5_dp * yD_RAM(:)/yH_RAM(:) ! Re-normalize

           CALL GSL_Interpolation_1D('Cubic',PA(NPA:1:-1),h_Cart(i,j,NPA:1:-1),&
                                     PAbn(NPA-1:2:-1),h_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
           CALL GSL_Interpolation_1D('Cubic',PA(NPA:1:-1),I_Cart(i,j,NPA:1:-1),&
                                     PAbn(NPA-1:2:-1),I_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
        enddo
        do j = 2,nT ! Now apply the one MLT to the rest of the MLTs
          bZEq_Cart(i,j) = bZEq_Cart(i,1)
          I_cart(i,j,:) = I_cart(i,1,:)
          H_cart(i,j,:) = H_cart(i,1,:)
          HDens_cart(i,j,:) = HDens_cart(i,1,:)
          h_cart_interp(i,j,:) = h_cart_interp(i,1,:)
          I_cart_interp(i,j,:) = I_cart_interp(i,1,:)
        enddo
     enddo

     h_Cart_interp(:,:,NPA) = h_Cart_interp(:,:,NPA-1)
     I_Cart_interp(:,:,NPA) = I_Cart_interp(:,:,NPA-1)
     h_Cart_interp(:,:,1) = h_Cart_interp(:,:,2)
     I_Cart_interp(:,:,1) = I_Cart_interp(:,:,2)

     ! End RAM Timing
     call system_clock(time1,clock_rate,clock_max)
     stoptime=time1/real(clock_rate,dp)
     write(*,*) 'RAM h and I: ',stoptime-starttime

     ! Make sure all H and I integrals have been filled in correctly
     IF (MINVAL(h_Cart) < 0._dp .OR. MINVAL(I_Cart)<0._dp) THEN
        PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart)
        PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart)
        !   STOP
     END IF

     DO I=2,NR+1
        DO J=1,NT
           DO L=1,NPA
              ! SCB Variables -> RAM Variables
              FNHS(I,J,L)   = h_Cart(I-1,J,L)
              FNIS(I,J,L)   = i_Cart(I-1,J,L)
              BNES(I,J)     = bZEq_Cart(I-1,J)
              HDNS(I,J,L)   = hdens_Cart(I-1,J,L)
              BOUNHS(I,J,L) = h_Cart_interp(I-1,J,L)
              BOUNIS(I,J,L) = i_Cart_interp(I-1,J,L)
              !
              dIdt(I,J,L)   = 0._dp
              dIbndt(I,J,L) = 0._dp
           ENDDO
           IF (J.LT.8.or.J.GT.18) THEN
              DO L=15,2,-1
                 if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
                    FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
                 endif
              ENDDO
           ENDIF
           BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
           dBdt(I,J) = 0._dp
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

     DEALLOCATE(length,r0,BeqDip,distance,H_value,I_value,HDens_value,yI,yH,yD, &
                bfmirror,BeqDip_Cart,bfMirror_RAM,yI_RAM,yH_RAM,yD_RAM,bbx,bby, &
                bbz,bb,dd,cVal,xx,yy,zz,BzeqDiff_Cart,BNESPrev,FNISPrev,ScaleAt,&
                BOUNISPrev)

     RETURN 

  CASE default  ! Regular calculation of h, I, bounce-averaged charge xchange etc.
!!!!! BEGIN SCB INTEGRAl CALCULATION
     write(*,*) 'Calculating h and I integrals'

     ! Start timing
     call system_clock(time1,clock_rate,clock_max)
     starttime=time1/real(clock_rate,dp)

     ! Start parallel loop
!$OMP PARALLEL DO
     RhoLoop: DO j = 2,npsi
        PhiLoop: DO k = 2,nzeta
           r0(j, k) = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
           !if (r0(j,k).gt.7.25) then ! Don't do the exterior field lines outside RAM domain
           !   I_value(j,k,:) = 0._dp
           !   H_value(j,k,:) = 0._dp
           !   HDens_value(j,k,:) = 0._dp
           !   cycle PhiLoop
           !endif
           DO i = 1, nthe
              distance(i,j,k) = SQRT(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2) ! Distance from center of Earth
           END DO
           ! Compute length of field line (j,k)
           length(j,k) = 0._dp
           DO i = 1, nthe-1
              length(j,k) = length(j,k) + SQRT((x(i+1,j,k)-x(i,j,k))**2 + (y(i+1,j,k)-y(i,j,k))**2 + (z(i+1,j,k)-z(i,j,k))**2)
           END DO

           bfmirror(j,k,:) = bf(nThetaEquator,j,k)/(1._dp - mu**2)
           bfmirror(j,k,NPA) = bf(nthe,j,k)

           CALL GSL_Integration_hI(bfMirror(j,k,:),chiVal(:),bf(:,j,k),distance(:,j,k), &
                                   yI(j,k,:), yH(j,k,:), yD(j,k,:))
           do L = 2,NPA
              if ((yI(j,k,L).ne.yI(j,k,L)).or.(yI(j,k,L)+1.eq.yI(j,k,L))) yI(j,k,L) = yI(j,k-1,L)
              if ((yH(j,k,L).ne.yH(j,k,L)).or.(yH(j,k,L)+1.eq.yH(j,k,L))) yH(j,k,L) = yH(j,k-1,L)
              if ((yD(j,k,L).ne.yD(j,k,L)).or.(yD(j,k,L)+1.eq.yD(j,k,L))) yD(j,k,L) = yD(j,k-1,L)
           enddo
           I_value(j,k,:) = (length(j,k)/(pi_d*r0(j,k))) * yI(j,k,:)/SQRT(Bfmirror(j,k,:))
           H_value(j,k,:) = (length(j,k)/(pi_d*2*r0(j,k))) * yH(j,k,:)*SQRT(Bfmirror(j,k,:))
           HDens_value(j,k,:) = 1.E5_dp * yD(j,k,:)/yH(j,k,:) ! Re-normalize

        END DO PhiLoop
     END DO RhoLoop
!$OMP END PARALLEL DO

     ! Handle periodicity
     I_value(:,1,:) = I_value(:,nzeta,:)
     H_value(:,1,:) = H_value(:,nzeta,:)

     ! End Timing
     call system_clock(time1,clock_rate,clock_max)
     stoptime=time1/real(clock_rate,dp)
     write(*,*) 'SCB h and I: ',stoptime-starttime
!!!!! END SCB INTEGRAl CALCULATION

!!!!! BEGIN SCB -> RAM CONVERSION
     ! Check if RAM point is inside or outside SCB domain
     ! If inside then interpolate to the the point
     ! If outside then use SCB magnetic field model to calculate h and I
     DIR = -1.0
     DSMAX = 0.1
     ER = 0.001
     RLIM = 20.0
     ScaleAt = 0._dp
     scalingI = 0._dp
     scalingH = 0._dp
     scalingD = 0._dp

     ! Start RAM Timing
     call system_clock(time1,clock_rate,clock_max)
     starttime=time1/real(clock_rate,dp)

     ! Interpolation will be for B - Bdip
     beqdip(1:npsi,2:nzeta)=b0dip/(x(nThetaEquator,1:npsi,2:nzeta)**2+y(nThetaEquator,1:npsi,2:nzeta)**2)**1.5
     beqdip_Cart(1:nR) = b0dip/radRaw(1:nR)**3

     ! Loop over RAM grid
     do j = 1,nT
        do i = 1,nR
           xo = radRaw(i) * COS(azimRaw(j)*2._dp*pi_d/24._dp - pi_d)
           yo = radRaw(i) * SIN(azimRaw(j)*2._dp*pi_d/24._dp - pi_d)

           ! Convex hull calculation
           wn = 0 ! wn =/= 0 means the point is inside the SCB domain
           DO k = 1,nzeta
              yn = y(nThetaEquator,npsi,k)
              yp = y(nThetaEquator,npsi,k+1)
              xn = x(nThetaEquator,npsi,k)
              xp = x(nThetaEquator,npsi,k+1)
              if (yn.le.yo) then
                 if (yp.gt.yo) then
                    if (((xp-xn)*(yo-yn)-(yp-yn)*(xo-xn)).gt.0) wn = wn + 1
                 endif
              else
                 if (yp.le.yo) then
                    if (((xp-xn)*(yo-yn)-(yp-yn)*(xo-xn)).lt.0) wn = wn - 1
                 endif
              endif
           ENDDO
           if ((wn.ne.0).or.(NameBoundMag.eq.'DIPS')) then
              ! Interpolate from SCB -> RAM grid
              DO L = 2,NPA
                 CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                           h_Value(:,2:nzeta,L), xo, yo, h_Cart(i,j,L),GSLerr)
                 CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                           I_Value(:,2:nzeta,L), xo, yo, I_Cart(i,j,L),GSLerr)
                 CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                           hDens_Value(:,2:nzeta,L), xo, yo, hDens_Cart(i,j,L),GSLerr)
              ENDDO
              CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                        bZ(nThetaEquator,:,2:nzeta)*bnormal-beqdip(:,2:nzeta), &
                                        xo, yo, bZEqDiff_Cart(i,j),GSLerr)
              CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                        fluxVolume(:,2:nzeta)/bnormal, xo, yo, flux_vol_Cart(i,j),GSLerr)
              bZEq_Cart(i,j) = bZEqDiff_Cart(i,j) + beqdip_Cart(i)
           else
              Pdyn = PARMOD(1)
              BzIMF = PARMOD(4)
              call SHUETAL_MGNP_08(PDyn,-1.0_dp,BzIMF,xo,yo,0.0_dp,XMGNP,YMGNP,ZMGNP,DIST,ID)
              if ((ID.lt.0).or.(DIST.lt.DL1)) then
                 bZEq_Cart(i,j) = -1.0
                 I_cart(i,j,:) = 1.0
                 H_cart(i,j,:) = 1.0
                 HDens_cart(i,j,:) = 1.0
              else
                 ! Calculate h and I at RAM grid point
                 !! Get equatorial point
                 x0 = LZ(i+1) * COS(MLT(j)*2._dp*pi_d/24._dp - pi_d)
                 y0 = LZ(i+1) * SIN(MLT(j)*2._dp*pi_d/24._dp - pi_d)
                 z0 = 0._dp

                 !! Trace from equatorial point to pole then pole to other pole
                 CALL trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD,  &
                            xe,ye,ze,xtemp(:),ytemp(:),ztemp(:),LOUT,LMAX, &
                            bxtemp(:),bytemp(:),bztemp(:))
                 CALL trace(xe,ye,ze,-DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD, &
                            xf,yf,zf,xtemp(:),ytemp(:),ztemp(:),LOUT,LMAX, &
                            bxtemp(:),bytemp(:),bztemp(:))
                 ! We have a problem now, LOUT is often only ~30 which is to few
                 ! points to get a good integral solution. For now I will just fit
                 ! the new points to the same grid as SCB. Can look into changing
                 ! this later. -ME
                 dtemp(1) = 0._dp
                 do k = 2,LOUT
                    dtemp(k) = dtemp(k-1) + SQRT((xtemp(k)-xtemp(k-1))**2 &
                                                +(ytemp(k)-ytemp(k-1))**2 &
                                                +(ztemp(k)-ztemp(k-1))**2)
                 enddo
                 cVal(:) = thetaVal(:)*dtemp(LOUT)/pi_d
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),xtemp(1:LOUT),cVal(:),xx(:),GSLerr)
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),ytemp(1:LOUT),cVal(:),yy(:),GSLerr)
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),ztemp(1:LOUT),cVal(:),zz(:),GSLerr)
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),bxtemp(1:LOUT),cVal(:),bbx(:),GSLerr)
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),bytemp(1:LOUT),cVal(:),bby(:),GSLerr)
                 CALL GSL_Interpolation_1D('Cubic',dtemp(1:LOUT),bztemp(1:LOUT),cVal(:),bbz(:),GSLerr)

                 !! Get necessary information at each point along the trace
                 dd(:) = SQRT(xx(:)**2+yy(:)**2+zz(:)**2) ! Distance from center of earth
                 bb(:) = SQRT(bbx(:)**2+bby(:)**2+bbz(:)**2)/bnormal

                 length_RAM = 0._dp
                 do k = 2,nthe
                    length_RAM = length_RAM + SQRT((xx(k)-xx(k-1))**2 &
                                                 + (yy(k)-yy(k-1))**2 &
                                                 + (zz(k)-zz(k-1))**2)
                 enddo
                 bZEq_Cart(i,j) = bbz(nThetaEquator)
                 cVal(:) = thetaVal(:)
                 r0_RAM = SQRT(xx(nThetaEquator)**2+yy(nThetaEquator)**2)

                 bfmirror_RAM(:) = bb(nThetaEquator)/(1._dp - mu**2)
                 bfmirror_RAM(NPA) = bb(nthe)

                 CALL GSL_Integration_hI(bfMirror_RAM(:),cVal(:),bb(:),dd(:), &
                                         yI_RAM(:),yH_RAM(:),yD_RAM(:))

                 I_cart(i,j,:) = (length_RAM/(pi_d*r0_RAM)) * yI_RAM(:)/SQRT(Bfmirror_RAM(:))
                 H_cart(i,j,:) = (length_RAM/(pi_d*2*r0_RAM)) * yH_RAM(:)*SQRT(Bfmirror_RAM(:))
                 HDens_cart(i,j,:) = 1.E5_dp * yD_RAM(:)/yH_RAM(:) ! Re-normalize

                 if (ScaleAt(j).eq.0) ScaleAt(j) = i
              endif
           endif
        enddo
     enddo

     !! Scale based on outer SCB boundary
     !Moved to it's own loop in case we want to add different scaling options later -ME
     do j = 1,nT
        if (ScaleAt(j).ne.0) then
           ii = ScaleAt(j)
           do L = 2, NPA
              I_Temp = I_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1)) &
                                         /(Lz(ii-2)-Lz(ii-1)) &
                                         *(I_Cart(ii-2,j,L)-I_Cart(ii-1,j,L))
              if (I_Temp.le.0) then
                 scalingI = I_Cart(ii-1,j,L)/I_Cart(ii,j,L)
              else
                 scalingI = I_Temp/I_cart(ii,j,L)
              endif
              !if (scalingI.le.0) write(*,*) j,L,scalingI,I_cart(ii,j,L),I_cart(ii-1,j,L),I_cart(ii-2,j,L)

              H_Temp = H_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1)) &
                                         /(Lz(ii-2)-Lz(ii-1)) &
                                         *(H_Cart(ii-2,j,L)-H_Cart(ii-1,j,L))
              if (H_Temp.le.0) then
                 scalingH = H_Cart(ii-1,j,L)/H_Cart(ii,j,L)
              else
                 scalingH = H_Temp/H_cart(ii,j,L)
              endif
              !if (scalingH.le.0) write(*,*) j,L,scalingH,H_cart(ii,j,L),H_cart(ii-1,j,L),H_cart(ii-2,j,L)

              D_Temp = HDens_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1)) &
                                             /(Lz(ii-2)-Lz(ii-1)) &
                                             *(HDens_Cart(ii-2,j,L)-HDens_Cart(ii-1,j,L))
              if (D_Temp.le.0) then
                 scalingD = HDens_Cart(ii-1,j,L)/HDens_Cart(ii,j,L)
              else
                 scalingD = D_Temp/HDens_cart(ii,j,L)
              endif
              !if (scalingD.le.0) write(*,*) j,L,scalingD,HDens_cart(ii,j,L),HDens_cart(ii-1,j,L),HDens_cart(ii-2,j,L)

              do i = ii, nR
                 I_cart(i,j,L) = I_cart(i,j,L)*scalingI
                 H_cart(i,j,L) = H_cart(i,j,L)*scalingH
                 HDens_cart(i,j,L) = HDens_cart(i,j,L)*scalingD
              enddo
           enddo
        endif
     enddo

     I_cart(:,:,3) = 0.99*I_cart(:,:,4)
     I_cart(:,:,2) = 0.99*I_cart(:,:,3)
     I_cart(:,:,1) = 0._dp
     H_cart(:,:,3) = 0.99*H_cart(:,:,4)
     H_cart(:,:,2) = 0.99*H_cart(:,:,3)
     H_cart(:,:,1) = 0.99*H_cart(:,:,2)
     HDens_cart(:,:,3) = 0.999*HDens_cart(:,:,4)
     HDens_cart(:,:,2) = 0.999*HDens_cart(:,:,3)
     HDens_cart(:,:,1) = 0.999*HDens_cart(:,:,2)

     ! End RAM Timing
     call system_clock(time1,clock_rate,clock_max)
     stoptime=time1/real(clock_rate,dp)
     write(*,*) 'RAM h and I: ',stoptime-starttime
!!!!! END SCB -> RAM CONVERSION

     ! Make sure all H and I integrals have been filled in correctly
     !! First Check for Negatives
     IF (MINVAL(h_Cart) < 0._dp .OR. MINVAL(I_Cart)<0._dp .OR. MINVAL(HDens_Cart)<0._dp) THEN
        PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart), minloc(H_Cart)
        PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart), minloc(I_Cart)
        print*, 'computehI: minval(D) = ', MINVAL(HDens_Cart), minloc(HDens_Cart)
        do j = 1, nT
           do i = 2, nR
              do L = 1, nPa
                 if (h_Cart(i,j,L).lt.0) h_Cart(i,j,L) = h_Cart(i-1,j,L)
                 if (I_Cart(i,j,L).lt.0) I_Cart(i,j,L) = I_Cart(i-1,j,L)
                 if (HDens_Cart(i,j,L).lt.0) HDens_Cart(i,j,L) = HDens_Cart(i-1,j,L)
              enddo
           enddo
        enddo
        PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart), minloc(H_Cart)
        PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart), minloc(I_Cart)
        print*, 'computehI: minval(D) = ', MINVAL(HDens_Cart), minloc(HDens_Cart)
     END IF
     !! Now check for bad integrals
     IF (MAXVAL(h_Cart) > 3._dp) THEN
        PRINT*, 'computehI: maxval(h) = ', MAXVAL(h_Cart), maxloc(H_Cart)
        PRINT*, 'computehI: maxval(I) = ', MAXVAL(I_Cart), maxloc(I_Cart)
        do j = 1, nT
           do i = 1, nR
              do L = 2, nPa
                 if (h_Cart(i,j,L).gt.3) then 
                    h_Cart(i,j,L) = h_Cart(i,j,L-1)
                    hDens_Cart(i,j,L) = hDens_Cart(i,j,L-1)
                 endif
              enddo
           enddo
        enddo
     ENDIF

     ! Cubic GSL interpolation with natural boundaries to get h and I at muboun 
     DO j = 1, nT
        DO i = 1, nR
           CALL GSL_Interpolation_1D('Cubic',PA(NPA:1:-1),h_Cart(i,j,NPA:1:-1),&
                                     PAbn(NPA-1:2:-1),h_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
           if (GSLerr.ne.0) then
              write(*,*) "  ModRamScb: Issue calculating BOUNHS; i,j = ", i, j
           endif
           CALL GSL_Interpolation_1D('Cubic',PA(NPA:1:-1),I_Cart(i,j,NPA:1:-1),&
                                     PAbn(NPA-1:2:-1),I_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
           if (GSLerr.ne.0) then
              write(*,*) "  ModRamScb: Issue calculating BOUNIS; i,j = ", i, j
           endif
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
           if (bZEq_Cart(I-1,J).lt.0) then
              FNIS(I-1,J,:) = -1.0
              FNHS(I-1,J,:) = -1.0
              HDNS(I-1,J,:) = -1.0
              BOUNHS(I-1,J,:) = -1.0
              BOUNIS(I-1,J,:) = -1.0
              dIdt(I,J,:) = -1.0
              dIbndt(I,J,:) = -1.0
              BNES(I,J) = -1.0
              dBdt(I,J) = -1.0
           else
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
                 if ((DthI.eq.0).or.(FNISPrev(I,J,L).lt.0)) then
                    dIdt(I,J,L)   = 0._dp
                    dIbndt(I,J,L) = 0._dp
                 else
                    dIdt(I,J,L)   = (FNIS(I,J,L)-FNISPrev(I,J,L))/DThI
                    dIbndt(I,J,L) = (BOUNIS(I,J,L)-BOUNISPrev(I,J,L))/DThI
                 endif
              ENDDO
              IF (J.LT.8.or.J.GT.18) THEN
                 DO L=15,2,-1
                    if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
                       FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
                    endif
                 ENDDO
              ENDIF
              BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
              if ((DthI.eq.0).or.(BNESPrev(I,J).lt.0)) then
                 dBdt(I,J) = 0._dp
              else
                 dBdt(I,J) = (BNES(I,J)-BNESPrev(I,J))/DThI
              endif
           endif
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

     DEALLOCATE(length,r0,BeqDip,distance,H_value,I_value,HDens_value,yI,yH,yD, &
                bfmirror,BeqDip_Cart,bfMirror_RAM,yI_RAM,yH_RAM,yD_RAM,bbx,bby, &
                bbz,bb,dd,cVal,xx,yy,zz,BzeqDiff_Cart,BNESPrev,FNISPrev,ScaleAt,&
                BOUNISPrev)

  END SELECT Operational_or_research

  RETURN

END SUBROUTINE computehI

END MODULE ModRamScb
