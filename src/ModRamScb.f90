!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModRamScb
! Contains subroutines and functions necessary to couple RAM and SCB

  use nrtype, ONLY: DP

  implicit none

! Index that maps the pitch angle on each point of a field line to the equatorial (RAM) one
  INTEGER, allocatable :: INDEXPA(:,:,:,:)
  REAL(DP), allocatable :: flux3DEQ(:,:,:,:,:)

  contains

!==================================================================================================
  subroutine ramscb_allocate
    ! Allocates arrays needed for storing the 3D flux mapping  
    use ModRamGrids, ONLY: NS, NE, NPA
    use ModScbGrids, ONLY: nthe, npsi, nzeta
  
    implicit none
  
    ALLOCATE(INDEXPA(nthe,npsi,nzeta,NPA), Flux3DEq(NS,npsi,nzeta,NE,NPA))
    INDEXPA = 0.0; Flux3DEq = 0.0
  
    return
  
  end subroutine ramscb_allocate

!==================================================================================================
  subroutine ramscb_deallocate
    ! Deallocates allocated arrays
    implicit none
  
    DEALLOCATE(INDEXPA, Flux3DEq)
  
  end subroutine ramscb_deallocate

!==================================================================================================
  subroutine Compute3DFlux
    ! Convert the distribution function to a 3D flux
    use ModRamVariables, ONLY: FLUX, MU, FFACTOR, FNHS, F2
    use ModRamGrids,     ONLY: nR, nT, nPa, nE, radout, RadiusMin, nS
  
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables, ONLY: bf, radGrid, angleGrid, radRaw, azimRaw, &
                               nThetaEquator, x, y
  
    use ModRamGSL,       ONLY: GSL_Interpolation_2D
    use ModScbFunctions, ONLY: locate
  
    use nrtype, ONLY: DP, pi_d, twopi_d
  
    implicit none
  
    integer :: i, j, k, L, iS, GSLErr
    real(DP) :: MuEq, radius, angle
 
    ! Convert RAM normalized distribution function to 2D Flux
    DO iS = 1,nS
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

    ! RAM Grid 
    DO j = 0,nR
       radRaw(j) = RadiusMin + (radOut-RadiusMin) * REAL(j,DP)/REAL(nR,DP)
    END DO
    DO k = 1,nT
       azimRaw(k) = 24.0 * REAL(k-1,DP)/REAL(nT-1,DP)
    END DO
 
    ! SCB Grid
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
    DO iS = 1, nS ! ions & electrons
       DO L = 1, NPA
          DO K = 1, NE
             !Cubic GSL interpolation
             CALL GSL_Interpolation_2D(radRaw(1:nR), azimRaw*pi_d/12.0, FLUX(iS,1:nR,1:nT,K,L), &
                                       radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), &
                                       flux3DEQ(iS,1:npsi,2:nzeta,K,L),GSLerr)
          END DO
       END DO
    END DO
    ! Periodicity
    flux3DEQ(:,:,1,:,:) = flux3DEQ(:,:,nzeta,:,:)
  
    DO k = 2, nzeta
       DO j = 1, npsi
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

!==================================================================================================
  SUBROUTINE computehI(iter)
    ! Calculates the h and I integral and bounce averaged densities by
    ! converting the SCB coordinates to RAM coordinates and then computing the
    ! interals 
    use ModRamVariables, ONLY: FNHS, FNIS, BNES, HDNS, dBdt, dHdt, &
                               dIdt, dIbndt, BOUNHS, BOUNIS, EIR, EIP, flux_volume, &
                               LZ, MU, MLT, PAbn, PA, DL1, outsideMGNP, xRAM, yRAM, zRAM
    use ModRamTiming,    ONLY: TimeRamElapsed, TOld, TimeRamNow
    use ModRamConst,     ONLY: b0dip
    use ModRamParams,    ONLY: NameBoundMag, verbose, checkMGNP, integral_smooth, densityMode
    use ModRamGrids,     ONLY: nR, nT, nPa, radiusMax, radiusMin
  
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbParams,    ONLY: method, constTheta
    use ModScbVariables, ONLY: bf, chiVal, x, y, z, bnormal, h_Cart, I_Cart, h_Cart_interp, &
                               i_Cart_interp, bZEq_Cart, flux_vol_cart, bZ, &
                               hdens_Cart, radRaw, azimRaw, nThetaEquator, fluxVolume, &
                               thetaVal, psi, alfa
 
    use ModRamCouple,    ONLY: BLines_DIII, nPoints, IsClosed_II

    use ModRamGSL,       ONLY: GSL_Interpolation_2D, GSL_Interpolation_1D, &
                               GSL_Integration_hI, GSL_Smooth_1D, GSL_BounceAverage
    use ModRamFunctions, ONLY: RamFileName, FUNT, FUNI
    use ModScbFunctions, ONLY: extap, locate, get_dipole_lines
    use ModScbIO,        ONLY: trace, PARMOD, IOPT
    use gaussian_filter, only: gaussian_kernel, convolve
  
    use nrtype,    ONLY: DP, pi_d, twopi_d
  
    use ModTimeConvert, ONLY: n_day_of_year

    implicit none
  
    INTEGER, INTENT(IN) :: iter
    INTEGER :: i, L, ii, GSLerr
    
    ! Variables for timing
    integer :: time1, clock_rate, clock_max
    real(dp) :: starttime,stoptime
  
    ! Variables for SCB
    REAL(DP) :: DthI
    REAL(DP), ALLOCATABLE :: length(:,:), r0(:,:), BeqDip(:,:)
    REAL(DP), ALLOCATABLE :: distance(:,:,:), H_value(:,:,:), I_value(:,:,:), &
                             Hdens_value(:,:,:), yI(:,:,:), yH(:,:,:), yD(:,:,:), &
                             bfmirror(:,:,:), density(:,:,:)
  
    ! Variables for RAM
    integer  :: nEquator
    integer, ALLOCATABLE :: ScaleAt(:), outsideSCB(:,:)
    REAL(DP) :: t0, t1, rt, tt, zt
    REAL(DP), ALLOCATABLE :: bbx(:), bby(:), bbz(:), xx(:), yy(:), zz(:), cval(:)
    REAL(DP), ALLOCATABLE :: BNESPrev(:,:), FNISPrev(:,:,:), BOUNISPrev(:,:,:), &
                             FNHSPrev(:,:,:), BOUNHSPrev(:,:,:), dHbndt(:,:,:)
    REAL(DP), ALLOCATABLE :: bRAM(:,:,:)
    REAL(DP) :: scalingI, scalingH, scalingD, I_Temp, H_Temp, D_Temp
 
    REAL(DP), ALLOCATABLE :: kernel(:,:), output(:,:)
 
    ! Variables for Tracing
    INTEGER :: LMAX, LOUT, nSWMF
    INTEGER :: ID
    REAL(DP) :: x0, y0, z0, xe, ye, ze, xf, yf, zf
    REAL(DP) :: ER, DSMAX, RLIM, DIR
    REAL(DP) :: PDyn, BzIMF, DIST, XMGNP, YMGNP, ZMGNP
    REAL(DP), DIMENSION(1000) :: xtemp, ytemp, ztemp, bxtemp, bytemp, bztemp, dtemp
  
    integer, save :: j, k, wn
    REAL(DP), save :: xo, xn, xp, yo, yn, yp, psiRAM, alphaRAM
    !$OMP THREADPRIVATE(k, j, wn, xo, xn, xp, yo, yn, yp, psiRAM, alphaRAM)

    ! Don't need to run if using Dipole magnetic field boundary unless it is run from the initialization step
    if ((NameBoundMag.eq.'DIPL').and.(iter.ne.0)) return

    ! If we need to do any Geopack tracing, we need to make sure RECALC has been called
    call RECALC_08(TimeRamNow%iYear,n_day_of_year(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay), &
                   TimeRamNow%iHour,TimeRamNow%iMinute,TimeRamNow%iSecond,-400._dp,0._dp,0._dp)
    clock_rate = 1000
    clock_max = 100000
    LMAX = 1000
  
    !!! Initialize Allocatable Arrays
    ALLOCATE(ScaleAt(nT), outsideSCB(nR,nT), bbx(nthe), bby(nthe), bbz(nthe), xx(nthe), &
             yy(nthe), zz(nthe), cval(nthe), distance(nthe,nR,nT), length(nR,nT), &
             r0(nR,nT), bfMirror(nR,nT,nPa), yI(nR,nT,NPA), yH(nR,nT,NPA), yD(nR,nT,NPA), &
             BNESPrev(nR+1,nT), FNHSPrev(nR+1,nT,nPa),FNISPrev(nR+1,nT,nPa), &
             BOUNISPrev(nR+1,nT,nPa), BOUNHSPrev(nR+1,nT,nPa), dHbndt(nR+1,nT,nPa), &
             bRAM(nthe,nR,nT), density(nthe,nR,nT))
    !!!

    h_Cart = 0._dp; I_Cart = 0._dp
    bZEq_Cart = 0._dp; flux_vol_Cart = 0._dp
    outsideMGNP = 0; outsideSCB  = 0; ScaleAt = 0
    xRAM = 0._dp; yRAM = 0._dp; zRAM = 0._dp; bRAM = 0._dp

    ! Start timing
    call system_clock(time1,clock_rate,clock_max)
    starttime=time1/real(clock_rate,dp)

    SELECT CASE (NameBoundMag)
    CASE('DIPL') ! Dipole without SCB calculation (only run on initialization step)
       call get_dipole_lines(radiusMin,radiusMax,constTheta,nthe,nR,nT,xRAM,yRAM,zRAM,bRAM,.true.)
       bRAM = bRAM*(b0dip/bnormal)

    CASE default  ! Convert SCB field lines to RAM field lines
       ! Start timing
       if (verbose) write(*,'(1x,a)',ADVANCE='NO') 'Converting SCB field lines to RAM field lines'
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

  !$OMP PARALLEL DO
       do i = 1,nR
          do j = 1,nT
             ! Convex hull calculation
             xo = LZ(i+1) * COS(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             yo = LZ(i+1) * SIN(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             wn = 0 ! wn =/= 0 means the point is inside the SCB domain
             DO k = 1,nzeta
                yn = y(nThetaEquator,npsi-1,k)
                yp = y(nThetaEquator,npsi-1,k+1)
                xn = x(nThetaEquator,npsi-1,k)
                xp = x(nThetaEquator,npsi-1,k+1)
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
             if (abs(wn) > 0) then ! Boundary overlap
                ! Step 1: Get psi and alpha value for the equatorial point
                !! alpha should just be the MLT stuff rotated
                alphaRAM = MLT(j)*pi_d/12._dp + pi_d
                if (alphaRAM > twopi_d) alphaRAM = alphaRAM - twopi_d
                !! psi can be gotten from interpolation
                CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                          psi(nThetaEquator,:,2:nzeta), xo, yo, psiRAM, GSLerr)
                ! Step 2: Interpolate from old (psi, alfa) grid onto psiRAM, alphaRAM point
                do k = 1,nthe
                   CALL GSL_Interpolation_2D(psi(k,:,2:nzeta), alfa(k,:,2:nzeta), x(k,:,2:nzeta), &
                                             psiRAM, alphaRAM, xRAM(k,i,j), GSLerr)
                   CALL GSL_Interpolation_2D(psi(k,:,2:nzeta), alfa(k,:,2:nzeta), y(k,:,2:nzeta), &
                                             psiRAM, alphaRAM, yRAM(k,i,j), GSLerr)
                   CALL GSL_Interpolation_2D(psi(k,:,2:nzeta), alfa(k,:,2:nzeta), z(k,:,2:nzeta), &
                                             psiRAM, alphaRAM, zRAM(k,i,j), GSLerr)
                   CALL GSL_Interpolation_2D(psi(k,:,2:nzeta), alfa(k,:,2:nzeta), bf(k,:,2:nzeta), &
                                             psiRAM, alphaRAM, bRAM(k,i,j), GSLerr)

                enddo
             else
                outsideSCB(i,j) = 1
             endif
          enddo
       enddo
  !$OMP END PARALLEL DO

       ! End Timing
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       if (verbose) write(*,'(a,1x,F6.2,1x,a)') ': Completed in', stoptime-starttime, 'seconds'

       ! If a point on the RAM grid lies outside the SCB grid then we need to
       ! calculate the scaling parameters (if still inside the magnetopause) or
       ! set the outsideMGNP flag to 1 so we can track it in RAM
       do j = 1,nT
          do i = 1,nR
             xo = LZ(i+1) * COS(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             yo = LZ(i+1) * SIN(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             if (outsideSCB(i,j) == 1) then
                if (ScaleAt(j) == 0) ScaleAt(j) = i
                if (NameBoundMag.eq.'SWMF') then
                   outsideMGNP(i,j) = 1
                   cycle
                endif
                Pdyn = PARMOD(1)
                BzIMF = PARMOD(4)
                call SHUETAL_MGNP_08(PDyn,-1.0_dp,BzIMF,xo,yo,0.0_dp,XMGNP,YMGNP,ZMGNP,DIST,ID)
                if ((CheckMGNP).and.((ID.lt.0).or.(DIST.lt.DL1))) then
                   if (verbose) write(*,'(2x,a,3F6.2)') 'Point outside magnetopause', xo, yo, DIST
                   outsideMGNP(i,j) = 1
                   cycle
                else
                   !! Trace from equatorial point to pole then pole to other pole
                   CALL trace(xo,yo,0._dp,1.0_dp,xe,ye,ze,xtemp(:),ytemp(:),ztemp(:),LOUT,LMAX, &
                              bxtemp(:),bytemp(:),bztemp(:))
                   CALL trace(xe,ye,ze,-1.0_dp,xf,yf,zf,xtemp(:),ytemp(:),ztemp(:),LOUT,-LMAX, &
                              bxtemp(:),bytemp(:),bztemp(:))
                   ! If the tracing fails for some reason assume we are on open field lines
                   if (LOUT.ge.LMAX) then
                      outsideMGNP(i,j) = 1
                      cycle
                   endif
                endif
                dtemp(1) = 0._dp
                do k = 2,LOUT
                   dtemp(k) = dtemp(k-1) + SQRT((xtemp(k)-xtemp(k-1))**2 &
                                               +(ytemp(k)-ytemp(k-1))**2 &
                                               +(ztemp(k)-ztemp(k-1))**2)
                enddo

                cVal(:) = chiVal(:)*dtemp(LOUT)/pi_d
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),xtemp(1:LOUT),cVal(:),xRAM(:,i,j),GSLerr)
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),ytemp(1:LOUT),cVal(:),yRAM(:,i,j),GSLerr)
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),ztemp(1:LOUT),cVal(:),zRAM(:,i,j),GSLerr)
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),bxtemp(1:LOUT),cVal(:),bbx(:),GSLerr)
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),bytemp(1:LOUT),cVal(:),bby(:),GSLerr)
                CALL GSL_Interpolation_1D(dtemp(1:LOUT),bztemp(1:LOUT),cVal(:),bbz(:),GSLerr)
                bRAM(:,i,j) = SQRT(bbx(:)**2+bby(:)**2+bbz(:)**2)/bnormal
             endif
          enddo
       enddo
    END SELECT

    ! Calculations
    distance(:,:,:) = SQRT(xRAM(:,:,:)**2+yRAM(:,:,:)**2+zRAM(:,:,:)**2) ! Distance from center of earth

    ! Density to use for bounce averaging
    if (densityMode == "RAIRDEN") then
        density = 10**(13.326 - 3.6908*distance     &
                              + 1.1362*distance**2  &
                              - 0.16984*distance**3 &
                              + 0.009553*distance**4)
    endif

    !!!!! BEGIN INTEGRAl CALCULATION
    if (verbose) write(*,'(1x,a)',ADVANCE='NO') 'Calculating h and I integrals'
    ! Start timing
    call system_clock(time1,clock_rate,clock_max)
    starttime=time1/real(clock_rate,dp)

  !$OMP PARALLEL DO
    do i = 1, nR
       do j = 1, nT
          if (outsideMGNP(i,j) == 0) then
             length(i,j) = 0._dp
             do k = 2,nthe
                length(i,j) = length(i,j) + SQRT((xRAM(k,i,j)-xRAM(k-1,i,j))**2 &
                                               + (yRAM(k,i,j)-yRAM(k-1,i,j))**2 &
                                               + (zRAM(k,i,j)-zRAM(k-1,i,j))**2)
             enddo
             r0(i,j) = SQRT(xRAM(nThetaEquator,i,j)**2+yRAM(nThetaEquator,i,j)**2)
             if (abs(bRAM(nThetaEquator,i,j) - minval(bRAM(:,i,j))) > 1e-9) then
              if (2._dp*minval(bRAM(:,i,j))-bRAM(nThetaEquator,i,j) > 0._dp) then
                 bRAM(nThetaEquator,i,j) = 2._dp*minval(bRAM(:,i,j))-bRAM(nThetaEquator,i,j)
              else
                 bRAM(nThetaEquator,i,j) = minval(bRAM(:,i,j))-0.01
              endif
             endif

             bfmirror(i,j,1:NPA-1) = bRAM(nThetaEquator,i,j)/(1._dp -mu(1:NPA-1)**2)
             bfmirror(i,j,NPA)     = bRAM(nthe,i,j)

             CALL GSL_Integration_hI(bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yI(i,j,:), yH(i,j,:))
             CALL GSL_BounceAverage(density(:,i,j), bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yD(i,j,:))

             I_cart(i,j,:) = (length(i,j)/(pi_d*r0(i,j))) *yI(i,j,:)/SQRT(Bfmirror(i,j,:))
             H_cart(i,j,:) = (length(i,j)/(pi_d*2*r0(i,j))) *yH(i,j,:)*SQRT(Bfmirror(i,j,:))
             HDens_cart(i,j,:) = yD(i,j,:)/yH(i,j,:)
             bZEq_Cart(i,j) = bRAM(nThetaEquator,i,j)*bnormal

          end if
       enddo
    enddo
  !$OMP END PARALLEL DO

    ! End Timing
    call system_clock(time1,clock_rate,clock_max)
    stoptime=time1/real(clock_rate,dp)
    if (verbose) write(*,'(a,1x,F6.2,1x,a)') ': Completed in', stoptime-starttime, 'seconds'
    !!!!! END SCB INTEGRAl CALCULATION

    ! Scale based on outer SCB boundary. This is needed if the RAM grid has
    ! points that lay outside the SCB grid.
    scalingI = 0._dp
    scalingH = 0._dp
    scalingD = 0._dp
    do j = 2,nT
       if (ScaleAt(j).ne.0) then
          ii = ScaleAt(j)
          do L = 2, NPA
             I_Temp = I_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1))/(Lz(ii-2)-Lz(ii-1)) &
                                        *(I_Cart(ii-2,j,L)-I_Cart(ii-1,j,L))
             if (I_Temp.le.0) then
                scalingI = I_Cart(ii-1,j,L)/I_Cart(ii,j,L)
             else
                scalingI = I_Temp/I_cart(ii,j,L)
             endif
  
             H_Temp = H_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1))/(Lz(ii-2)-Lz(ii-1)) &
                                        *(H_Cart(ii-2,j,L)-H_Cart(ii-1,j,L))
             if (H_Temp.le.0) then
                scalingH = H_Cart(ii-1,j,L)/H_Cart(ii,j,L)
             else
                scalingH = H_Temp/H_cart(ii,j,L)
             endif
  
             D_Temp = HDens_cart(ii-1,j,L) + (Lz(ii)-Lz(ii-1))/(Lz(ii-2)-Lz(ii-1)) &
                                            *(HDens_Cart(ii-2,j,L)-HDens_Cart(ii-1,j,L))
             if (D_Temp.le.0) then
                scalingD = HDens_Cart(ii-1,j,L)/HDens_Cart(ii,j,L)
             else
                scalingD = D_Temp/HDens_cart(ii,j,L)
             endif

             do i = ii, nR
                if (outsideMGNP(i,j) == 0) then
                   ! If inside Magnetopause but outside SCB domain then scale h
                   ! and I integrals based on the last SCB radial point
                   I_cart(i,j,L) = I_cart(i,j,L)*scalingI
                   H_cart(i,j,L) = H_cart(i,j,L)*scalingH
                   HDens_cart(i,j,L) = HDens_cart(i,j,L)*scalingD
                   bZEq_Cart(i,j) = bZEq_Cart(i-1,j)
                else 
                   ! If outside Magnetopause then set h and I integrals to 
                   ! the previous radial point
                   I_cart(i,j,L) = I_cart(i-1,j,L)
                   H_cart(i,j,L) = H_cart(i-1,j,L)
                   HDens_cart(i,j,L) = HDens_cart(i-1,j,L)
                   bZEq_Cart(i,j) = bZEq_Cart(i-1,j)
                endif
             enddo
          enddo
       endif
    enddo
  
    ! Continuity across MLT of 0
    I_Cart(:,1,:) = I_Cart(:,nT,:)
    H_Cart(:,1,:) = H_Cart(:,nT,:)
    HDens_Cart(:,1,:) = HDens_Cart(:,nT,:)
    bZEq_Cart(:,1) = bZEq_Cart(:,nT)

    ! Near 90 degree pitch angle corrections
    I_cart(:,:,3) = 0.50*I_cart(:,:,4)
    I_cart(:,:,2) = 0.20*I_cart(:,:,3)
    I_cart(:,:,1) = 0._dp
    H_cart(:,:,3) = 0.99*H_cart(:,:,4)
    H_cart(:,:,2) = 0.99*H_cart(:,:,3)
    H_cart(:,:,1) = 0.99*H_cart(:,:,2)
    HDens_cart(:,:,3) = 0.999*HDens_cart(:,:,4)
    HDens_cart(:,:,2) = 0.999*HDens_cart(:,:,3)
    HDens_cart(:,:,1) = 0.999*HDens_cart(:,:,2)
  
    ! Make sure all H and I integrals have been filled in correctly
    !! First Check for Negatives
    IF (MINVAL(h_Cart) < 0._dp .OR. MINVAL(I_Cart)<0._dp .OR. MINVAL(HDens_Cart)<0._dp) THEN
       if (verbose) PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart), minloc(H_Cart)
       if (verbose) PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart), minloc(I_Cart)
       if (verbose) print*, 'computehI: minval(D) = ', MINVAL(HDens_Cart), minloc(HDens_Cart)
       do j = 1, nT
          do i = 2, nR
             do L = 1, nPa
                if (h_Cart(i,j,L).lt.0) h_Cart(i,j,L) = h_Cart(i-1,j,L)
                if (I_Cart(i,j,L).lt.0) I_Cart(i,j,L) = I_Cart(i-1,j,L)
                if (HDens_Cart(i,j,L).lt.0) HDens_Cart(i,j,L) = HDens_Cart(i-1,j,L)
             enddo
          enddo
       enddo
       if (verbose) PRINT*, 'computehI: minval(h) = ', MINVAL(h_Cart), minloc(H_Cart)
       if (verbose) PRINT*, 'computehI: minval(I) = ', MINVAL(I_Cart), minloc(I_Cart)
       if (verbose) print*, 'computehI: minval(D) = ', MINVAL(HDens_Cart), minloc(HDens_Cart)
    END IF
    !! Now check for bad integrals (integrals that are too large)
    do i = 1, nR
       do j = 1, nT
          do L = nPa-1,1,-1
             if (i_Cart(i,j,L).gt.i_Cart(i,j,L+1)) i_Cart(i,j,L) = 0.99*i_Cart(i,j,L+1)
             if (H_Cart(i,j,L).gt.H_Cart(i,j,L+1)) H_Cart(i,j,L) = 0.99*H_Cart(i,j,L+1)
             if (HDens_Cart(i,j,L).gt.HDens_Cart(i,j,L+1)) HDens_Cart(i,j,L) = 0.999*HDens_Cart(i,j,L+1)
          enddo
       enddo
    enddo
  
    ! Cubic GSL interpolation with natural boundaries to get h and I at muboun 
    DO j = 1, nT
       DO i = 1, nR
          CALL GSL_Interpolation_1D(PA(NPA:1:-1),h_Cart(i,j,NPA:1:-1),&
                                    PAbn(NPA-1:2:-1),h_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
          CALL GSL_Interpolation_1D(PA(NPA:1:-1),I_Cart(i,j,NPA:1:-1),&
                                    PAbn(NPA-1:2:-1),I_Cart_interp(i,j,NPA-1:2:-1),GSLerr)
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

    ! Smooth the h and I integrals and bounce averaged values using a gaussian filter
    if (integral_smooth) then
       allocate(output(nR,nT))
       do L = 2,NPA
          call gaussian_kernel(1.0, kernel)
          call convolve(h_Cart(:,:,L), kernel, output)
          h_Cart(:,:,L) = output

          call gaussian_kernel(1.0, kernel)
          call convolve(i_Cart(:,:,L), kernel, output)
          i_Cart(:,:,L) = output

          call gaussian_kernel(1.0, kernel)
          call convolve(h_Cart_interp(:,:,L), kernel, output)
          h_Cart_interp(:,:,L) = output

          call gaussian_kernel(1.0, kernel)
          call convolve(i_Cart_interp(:,:,L), kernel, output)
          i_Cart_interp(:,:,L) = output

          call gaussian_kernel(1.0, kernel)
          call convolve(hdens_Cart(:,:,L), kernel, output)
          hdens_Cart(:,:,L) = output
       ENDDO
       deallocate(output)
    endif

    ! Fill in the RAM variables and calculate all dt values
    DO I=2,NR+1
       DO J=1,NT
          BNESPrev(I,J) = BNES(I,J)
          DO L=1,NPA
             FNISPrev(I,J,L)   = FNIS(I,J,L)
             FNHSPrev(I,J,L)   = FNHS(I,J,L)
             BOUNISPrev(I,J,L) = BOUNIS(I,J,L)
             BOUNHSPrev(I,J,L) = BOUNHS(I,J,L)
             ! SCB Variables -> RAM Variables
             FNHS(I,J,L)   = h_Cart(I-1,J,L)
             FNIS(I,J,L)   = i_Cart(I-1,J,L)
             BNES(I,J)     = bZEq_Cart(I-1,J)
             HDNS(I,J,L)   = hdens_Cart(I-1,J,L)
             BOUNHS(I,J,L) = h_Cart_interp(I-1,J,L)
             BOUNIS(I,J,L) = i_Cart_interp(I-1,J,L)
             !
             if (abs(DthI).le.1e-9) then
                dIdt(I,J,L)   = 0._dp
                dHdt(I,J,L)   = 0._dp
                dIbndt(I,J,L) = 0._dp
                dHbndt(I,J,L) = 0._dp
             else
                dIdt(I,J,L)   = (FNIS(I,J,L) - FNISPrev(I,J,L))/DThI
                dHdt(I,J,L)   = (FNHS(I,J,L) - FNHSPrev(I,J,L))/DThI
                dIbndt(I,J,L) = (BOUNIS(I,J,L) - BOUNISPrev(I,J,L))/DThI
                dHbndt(I,J,L) = (BOUNHS(I,J,L) - BOUNHSPrev(I,J,L))/DthI
             endif
             !
          ENDDO
          BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
          if (abs(DthI).le.1e-9) then
             dBdt(I,J) = 0._dp
          else
             dBdt(I,J) = (BNES(I,J) - BNESPrev(I,J))/DThI
          endif
       ENDDO
    ENDDO
    DO J=1,NT ! use dipole B at I=1
       BNES(1,J) = 0.32/LZ(1)**3/1.e4
       dBdt(1,J) = 0._dp
       EIR(1,J)  = 0._dp
       EIP(1,J)  = 0._dp
       DO L=1,NPA
          FNHS(1,J,L)   = FNHS(2,J,L)
          FNIS(1,J,L)   = FNIS(2,J,L)
          BOUNHS(1,J,L) = BOUNHS(2,J,L)
          BOUNIS(1,J,L) = BOUNIS(2,J,L)
          HDNS(1,J,L)   = HDNS(2,J,L)
          dIdt(1,J,L)   = 0._dp
          dHdt(1,J,L)   = 0._dp
          dIbndt(1,J,L) = 0._dp
       ENDDO
    ENDDO
    TOld = TimeRamElapsed

    ! Check for NaN's
    do i = 2, nR+1
       do j = 1, nT
          do L = 1, nPa
             if (isnan(FNIS(I,J,L))) FNIS(I,J,L) = FNIS(I-1,J,L)
             if (isnan(FNHS(I,J,L))) FNHS(I,J,L) = FNHS(I-1,J,L)
             if (isnan(BOUNIS(I,J,L))) BOUNIS(I,J,L) = BOUNIS(I-1,J,L)
             if (isnan(BOUNHS(I,J,L))) BOUNHS(I,J,L) = BOUNHS(I-1,J,L)
             if (isnan(HDNS(I,J,L))) HDNS(I,J,L) = HDNS(I-1,J,L)
             if (isnan(dIdt(I,J,L))) dIdt(I,J,L) = 0._dp
             if (isnan(dIbndt(I,J,L))) dIbndt(I,J,L) = 0._dp
          enddo
       enddo
    enddo

    DEALLOCATE(length,r0,distance,yI,yH,yD,bfmirror,bbx,bby,bbz,cVal,xx,yy,zz, &
               BNESPrev,FNISPrev,ScaleAt,BOUNISPrev,FNHSPrev,BOUNHSPrev,dHbndt)
    DEALLOCATE(bRAM,outsideSCB,density)
 
    RETURN
  
  END SUBROUTINE computehI

!==================================================================================================
  SUBROUTINE computehI_test(iter)
    ! Calculates the h and I integral and bounce averaged densities by
    ! computing them on the SCB grid and then interpolating them onto the RAM
    ! grid
    use ModRamVariables, ONLY: FNHS, FNIS, BNES, HDNS, ODNS, NDNS, dBdt, dHdt, &
                               dIdt, dIbndt, BOUNHS, BOUNIS, EIR, EIP, flux_volume, Phi, &
                               LZ, MU, MLT, PAbn, PA, DL1, outsideMGNP, xRAM, yRAM, zRAM
    use ModRamTiming,    ONLY: TimeRamElapsed, TOld, TimeRamNow
    use ModRamConst,     ONLY: b0dip, RE
    use ModRamParams,    ONLY: NameBoundMag, verbose, checkMGNP, integral_smooth, densityMode
    use ModRamGrids,     ONLY: nR, nT, nPa, radiusMax, radiusMin
  
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbParams,    ONLY: method, constTheta
    use ModScbVariables, ONLY: bf, chiVal, x, y, z, bnormal, bZ, &
                               radRaw, azimRaw, nThetaEquator, fluxVolume, &
                               thetaVal, psi, alfa
 
    use ModRamCouple,    ONLY: BLines_DIII, nPoints, IsClosed_II

    use ModRamGSL,       ONLY: GSL_Interpolation_2D, GSL_Interpolation_1D, &
                               GSL_Integration_hI, GSL_Smooth_1D, GSL_BounceAverage
    use ModRamFunctions, ONLY: RamFileName, FUNT, FUNI
    use ModScbFunctions, ONLY: extap, locate, get_dipole_lines
    use ModScbIO,        ONLY: trace, PARMOD, IOPT
    use gaussian_filter, only: gaussian_kernel, convolve
  
    use nrtype,    ONLY: DP, pi_d, twopi_d, pio2_d
  
    use ModTimeConvert, ONLY: n_day_of_year

    implicit none
  
    INTEGER, INTENT(IN) :: iter
    INTEGER :: j, L, ii, GSLerr
    
    ! Variables for timing
    integer :: time1, clock_rate, clock_max
    real(dp) :: starttime,stoptime
  
    ! Variables for SCB
    REAL(DP) :: DthI
    REAL(DP), ALLOCATABLE :: BeqDip(:,:)
    REAL(DP), ALLOCATABLE :: distance(:,:,:)
  
    ! Variables for RAM
    integer  :: nEquator
    integer, ALLOCATABLE :: ScaleAt(:), outsideSCB(:,:)
    REAL(DP) :: t0, t1, rt, tt, zt, radius, alt, lat, lon, D(9), T(2)
    REAL(DP), ALLOCATABLE :: bbx(:), bby(:), bbz(:), xx(:), yy(:), zz(:), cval(:)
    REAL(DP), ALLOCATABLE :: BNESPrev(:,:), FNISPrev(:,:,:), BOUNISPrev(:,:,:), &
                             FNHSPrev(:,:,:), BOUNHSPrev(:,:,:), dHbndt(:,:,:)
    REAL(DP), ALLOCATABLE :: bRAM(:,:,:), HDen(:,:,:), ODen(:,:,:), NDen(:,:,:), &
                             I_cart(:,:,:), h_cart(:,:,:), HDens_cart(:,:,:), &
                             ODens_cart(:,:,:), NDens_cart(:,:,:), bZEq_cart(:,:)
    REAL(DP) :: scalingI, scalingH, scalingD, I_Temp, H_Temp, D_Temp
 
    REAL(DP), ALLOCATABLE :: kernel(:,:), output(:,:)
 
    ! Variables for Tracing
    INTEGER :: LMAX, LOUT, nSWMF
    INTEGER :: ID
    REAL(DP) :: x0, y0, z0, xe, ye, ze, xf, yf, zf
    REAL(DP) :: ER, DSMAX, RLIM, DIR
    REAL(DP) :: PDyn, BzIMF, DIST, XMGNP, YMGNP, ZMGNP
    REAL(DP), DIMENSION(1000) :: xtemp, ytemp, ztemp, bxtemp, bytemp, bztemp, dtemp
  
    integer, save :: i, k, wn
    REAL(DP), save :: xo, xn, xp, yo, yn, yp, psiRAM, alphaRAM, length, r0
    REAL(DP), allocatable :: bfMirror(:,:,:), yI(:,:,:), yH(:,:,:), yD(:,:,:)
    !$OMP THREADPRIVATE(k, i, wn, xo, xn, xp, yo, yn, yp, psiRAM, alphaRAM)
    !$OMP THREADPRIVATE(length, r0)

    ! Don't need to run if using Dipole magnetic field boundary unless it is run from the initialization step
    if ((NameBoundMag.eq.'DIPL').and.(iter.ne.0)) return
    call RECALC_08(TimeRamNow%iYear,n_day_of_year(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay), &
                   TimeRamNow%iHour,TimeRamNow%iMinute,TimeRamNow%iSecond,-400._dp,0._dp,0._dp)
    clock_rate = 1000
    clock_max = 100000
    LMAX = 1000
  
    !!! Initialize Allocatable Arrays
    ALLOCATE(I_cart(npsi,nzeta,nPa), h_cart(npsi,nzeta,nPa), bZEq_Cart(npsi,nzeta), &
             HDens_Cart(npsi,nzeta,nPa), ODens_Cart(npsi,nzeta,nPa), NDens_cart(npsi,nzeta,nPa))
    ALLOCATE(outsideSCB(nR,nT), BNESPrev(nR+1,nT), FNHSPrev(nR+1,nT,nPa),FNISPrev(nR+1,nT,nPa), &
             BOUNISPrev(nR+1,nT,nPa), BOUNHSPrev(nR+1,nT,nPa))
    !!!

    BNESPrev   = BNES
    FNISPrev   = FNIS
    FNHSPrev   = FNHS
    BOUNISPrev = BOUNIS
    BOUNHSPrev = BOUNHS

    ! Start timing
    call system_clock(time1,clock_rate,clock_max)
    starttime=time1/real(clock_rate,dp)

    SELECT CASE (NameBoundMag)
    CASE('DIPL') ! Dipole without SCB calculation (only run on initialization step)i
       ALLOCATE(HDen(nthe,nR,nT), ODen(nthe,nR,nT), NDen(nthe,nR,nT), distance(nthe,nR,nT))
       ALLOCATE(bfMirror(nR,nT,nPa), yI(nR,nT,nPa), yH(nR,nT,nPa), yD(nR,nT,nPa))
       ALLOCATE(bRAM(nthe,nR,nT))
       call get_dipole_lines(radiusMin,radiusMax,constTheta,nthe,nR,nT,xRAM,yRAM,zRAM,bRAM,.true.)
       bRAM = bRAM*(b0dip/bnormal)

       ! Density to use for bounce averaging
       distance = SQRT(xRAM(:,:,:)**2+yRAM(:,:,:)**2+zRAM(:,:,:)**2)
       select case(trim(densityMode))
       case ("RAIRDEN")
           HDen = 10**(13.326 - 3.6908*distance     &
                              + 1.1362*distance**2  &
                              - 0.16984*distance**3 &
                              + 0.009553*distance**4)
        case("MSIS")
           call METERC(.true.)
           do i = 1, nR
              do j = 2, nT
                 do k = 1, nthe
                    radius = distance(k,i,j)
                    alt = (radius*RE - RE)/1000._dp
                    if (alt < 80._dp) alt = 80._dp
                    lat = (ACOS(zRAM(k,i,j)/radius)-pio2_d)*180._dp/pi_d
                    lon = (ATAN2(yRAM(k,i,j),xRAM(k,i,j)))*180._dp/pi_d
                    CALL GTD7(62190,12*3600,alt,lat,lon,23.,150.,150.,4.,48,D,T)
                    ODen(k,i,j) = D(2)
                    HDen(k,i,j) = D(7)
                    NDen(k,i,j) = D(8)
                 enddo
              enddo
           enddo
           HDen(:,:,1) = HDen(:,:,nT)
           ODen(:,:,1) = ODen(:,:,nT)
           NDen(:,:,1) = NDen(:,:,nT)
        case default
           call CON_STOP("Unrecognized density mode")
        end select

  !$OMP PARALLEL DO
        do j = 2, nT
           do i = 1, nR
              length = 0._dp
              do k = 2,nthe
                 length = length + SQRT((xRAM(k,i,j)-xRAM(k-1,i,j))**2 &
                                      + (yRAM(k,i,j)-yRAM(k-1,i,j))**2 &
                                      + (zRAM(k,i,j)-zRAM(k-1,i,j))**2)
              enddo
              r0 = SQRT(xRAM(nThetaEquator,i,j)**2+yRAM(nThetaEquator,i,j)**2)
    
              bfmirror(i,j,1:NPA-1) = bRAM(nThetaEquator,i,j)/(1._dp-mu(1:NPA-1)**2)
              bfmirror(i,j,NPA)     = bRAM(nthe,i,j)
    
              CALL GSL_Integration_hI(bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yI(i,j,:), yH(i,j,:))
              CALL GSL_BounceAverage(HDen(:,i,j), bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yD(i,j,:))
    
              FNIS(i+1,j,:) = (length/(pi_d*r0))*yI(i,j,:)/SQRT(Bfmirror(i,j,:))
              FNHS(i+1,j,:) = (length/(pi_d*2*r0))*yH(i,j,:)*SQRT(Bfmirror(i,j,:))
              HDNS(i+1,j,:) = yD(i,j,:)/yH(i,j,:)
              BNES(i+1,j) = bRAM(nThetaEquator,i,j)*bnormal
    
              if (trim(densityMode) == "MSIS") then
                 CALL GSL_BounceAverage(ODen(:,i,j), bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yD(i,j,:))
                 ODNS(i+1,j,:) = yD(i,j,:)/yH(i,j,:)
                 CALL GSL_BounceAverage(NDen(:,i,j), bfMirror(i,j,:), chiVal(:), bRAM(:,i,j), yD(i,j,:))
                 NDNS(i+1,j,:) = yD(i,j,:)/yH(i,j,:)
              endif
           enddo
        enddo
  !$OMP END PARALLEL DO

    CASE default  ! Convert SCB field lines to RAM field lines
       ALLOCATE(HDen(nthe,npsi,nzeta), ODen(nthe,npsi,nzeta), NDen(nthe,npsi,nzeta), &
                distance(nthe,npsi,nzeta))
       ALLOCATE(bfMirror(npsi,nzeta,nPa), yI(npsi,nzeta,nPa), yH(npsi,nzeta,nPa), yD(npsi,nzeta,nPa))

       ! Start timing
       if (verbose) write(*,'(1x,a)',ADVANCE='NO') 'Calculating h and I integrals'
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

       ! Density to use for bounce averaging
       distance = SQRT(x(:,:,:)**2+y(:,:,:)**2+z(:,:,:)**2)
       select case(trim(densityMode))
       case ("RAIRDEN")
           HDen = 10**(13.326 - 3.6908*distance     &
                              + 1.1362*distance**2  &
                              - 0.16984*distance**3 &
                              + 0.009553*distance**4)
        case("MSIS")
           call METERC(.true.)
           do i = 1, npsi
              do j = 2, nzeta
                 do k = 1, nthe
                    radius = distance(k,i,j)
                    alt = (radius*RE - RE)/1000._dp
                    if (alt < 80._dp) alt = 80._dp
                    lat = (ACOS(z(k,i,j)/radius)-pio2_d)*180._dp/pi_d
                    lon = (ATAN2(y(k,i,j),x(k,i,j)))*180._dp/pi_d
                    CALL GTD7(62190,12*3600,alt,lat,lon,23.,150.,150.,4.,48,D,T)
                    ODen(k,i,j) = D(2)
                    HDen(k,i,j) = D(7)
                    NDen(k,i,j) = D(8)
                 enddo
              enddo
           enddo
           HDen(:,:,1) = HDen(:,:,nzeta)
           ODen(:,:,1) = ODen(:,:,nzeta)
           NDen(:,:,1) = NDen(:,:,nzeta)
        case default
           call CON_STOP("Unrecognized density mode")
        end select

  !$OMP PARALLEL DO
       do j = 2,nzeta
          do i = 1,npsi
             bfmirror(i,j,1:NPA-1) = minval(bf(:,i,j))/(1._dp -mu(1:NPA-1)**2)
             bfmirror(i,j,NPA)     = bf(nthe,i,j)

             CALL GSL_Integration_hI(bfMirror(i,j,:), chiVal(:), bf(:,i,j), yI(i,j,:), yH(i,j,:))
             CALL GSL_BounceAverage(HDen(:,i,j), bfMirror(i,j,:), chiVal(:), bf(:,i,j), yD(i,j,:))

             length = 0._dp
             do k = 2,nthe
                length = length + SQRT((x(k,i,j)-x(k-1,i,j))**2 &
                                     + (y(k,i,j)-y(k-1,i,j))**2 &
                                     + (z(k,i,j)-z(k-1,i,j))**2)
             enddo
             r0 = SQRT(x(nThetaEquator,i,j)**2+y(nThetaEquator,i,j)**2)

             I_cart(i,j,:) = (length/(pi_d*r0))*yI(i,j,:)/SQRT(Bfmirror(i,j,:))
             H_cart(i,j,:) = (length/(pi_d*2*r0))*yH(i,j,:)*SQRT(Bfmirror(i,j,:))
             HDens_cart(i,j,:) = yD(i,j,:)/yH(i,j,:)
             bZEq_Cart(i,j) = bf(nThetaEquator,i,j)*bnormal

             if (trim(densityMode) == "MSIS") then
                CALL GSL_BounceAverage(ODen(:,i,j), bfMirror(i,j,:), chiVal(:), bf(:,i,j), yD(i,j,:))
                ODens_cart(i,j,:) = yD(i,j,:)/yH(i,j,:)
                CALL GSL_BounceAverage(NDen(:,i,j), bfMirror(i,j,:), chiVal(:), bf(:,i,j), yD(i,j,:))
                NDens_cart(i,j,:) = yD(i,j,:)/yH(i,j,:)
             endif
          enddo
       enddo
  !$OMP END PARALLEL DO
       ! Continuity across MLT of 0
       I_Cart(:,1,:) = I_Cart(:,nzeta,:)
       H_Cart(:,1,:) = H_Cart(:,nzeta,:)
       HDens_Cart(:,1,:) = HDens_Cart(:,nzeta,:)
       ODens_Cart(:,1,:) = ODens_Cart(:,nzeta,:)
       NDens_Cart(:,1,:) = NDens_Cart(:,nzeta,:)
       bZEq_Cart(:,1) = bZEq_Cart(:,nzeta)

       ! Near 90 degree pitch angle corrections
       I_cart(:,:,3) = 0.50*I_cart(:,:,4)
       I_cart(:,:,2) = 0.20*I_cart(:,:,3)
       I_cart(:,:,1) = 0._dp
       H_cart(:,:,3) = 0.99*H_cart(:,:,4)
       H_cart(:,:,2) = 0.99*H_cart(:,:,3)
       H_cart(:,:,1) = 0.99*H_cart(:,:,2)
       HDens_cart(:,:,3) = 0.999*HDens_cart(:,:,4)
       HDens_cart(:,:,2) = 0.999*HDens_cart(:,:,3)
       HDens_cart(:,:,1) = 0.999*HDens_cart(:,:,2)
       ODens_cart(:,:,3) = 0.999*ODens_cart(:,:,4)
       ODens_cart(:,:,2) = 0.999*ODens_cart(:,:,3)
       ODens_cart(:,:,1) = 0.999*ODens_cart(:,:,2)
       NDens_cart(:,:,3) = 0.999*NDens_cart(:,:,4)
       NDens_cart(:,:,2) = 0.999*NDens_cart(:,:,3)
       NDens_cart(:,:,1) = 0.999*NDens_cart(:,:,2)

       ! End Timing
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       if (verbose) write(*,'(a,1x,F6.2,1x,a)') ': Completed in', stoptime-starttime, 'seconds'

  !$OMP PARALLEL DO
       do j = 2, nT
          do i = 1, nR
             xo = LZ(i+1) * COS(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             yo = LZ(i+1) * SIN(MLT(j)*2._dp*pi_d/24._dp - pi_d)
             wn = 0 ! wn =/= 0 means the point is inside the SCB domain
             DO k = 1,nzeta
                yn = y(nThetaEquator,npsi-1,k)
                yp = y(nThetaEquator,npsi-1,k+1)
                xn = x(nThetaEquator,npsi-1,k)
                xp = x(nThetaEquator,npsi-1,k+1)
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
             if (abs(wn) > 0) then ! Boundary Overlap
                do k = 1, nPa
                   CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                             I_Cart(:,2:nzeta,k), xo, yo, FNIS(i+1,j,k), GSLerr)
                   CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                             h_Cart(:,2:nzeta,k), xo, yo, FNHS(i+1,j,k), GSLerr)
                   CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                             HDens_Cart(:,2:nzeta,k), xo, yo, HDNS(i+1,j,k), GSLerr)
                   if (trim(densityMode) == "MSIS") then
                      CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                                ODens_Cart(:,2:nzeta,k), xo, yo, ODNS(i+1,j,k), GSLerr)
                      CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                                NDens_Cart(:,2:nzeta,k), xo, yo, NDNS(i+1,j,k), GSLerr)
                   endif
                enddo
                CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta), y(nThetaEquator,:,2:nzeta), &
                                          bZeq_Cart(:,2:nzeta), xo, yo, BNES(i+1,j), GSLerr)
             else
                outsideSCB(i,j) = 1
             endif
          enddo
       enddo
  !$OMP END PARALLEL DO
    END SELECT
    FNIS(:,1,:) = FNIS(:,nT,:)
    FNHS(:,1,:) = FNHS(:,nT,:)
    HDNS(:,1,:) = HDNS(:,nT,:)
    ODNS(:,1,:) = ODNS(:,nT,:)
    NDNS(:,1,:) = NDNS(:,nT,:)
    BNES(:,1) = BNES(:,nT)


    ! Make sure all H and I integrals have been filled in correctly
    !! First Check for Negatives
    IF (MINVAL(FNHS) < 0._dp .OR. MINVAL(FNIS)<0._dp .OR. MINVAL(HDNS)<0._dp) THEN
       if (verbose) PRINT*, 'computehI: minval(h) = ', MINVAL(FNHS), minloc(FNHS)
       if (verbose) PRINT*, 'computehI: minval(I) = ', MINVAL(FNIS), minloc(FNIS)
       if (verbose) print*, 'computehI: minval(D) = ', MINVAL(HDNS), minloc(HDNS)
       do j = 1, nT
          do i = 2, nR+1
             do L = 1, nPa
                if (FNHS(i,j,L).lt.0) FNHS(i,j,L) = FNHS(i-1,j,L)
                if (FNIS(i,j,L).lt.0) FNIS(i,j,L) = FNIS(i-1,j,L)
                if (HDNS(i,j,L).lt.0) HDNS(i,j,L) = HDNS(i-1,j,L)
             enddo
          enddo
       enddo
    END IF
  
    ! Cubic GSL interpolation with natural boundaries to get h and I at muboun 
    DO j = 2, nT
       DO i = 2, nR+1
          CALL GSL_Interpolation_1D(PA(NPA:1:-1), FNHS(i,j,NPA:1:-1),&
                                    PAbn(NPA-1:2:-1), BOUNHS(i,j,NPA-1:2:-1), GSLerr)
          CALL GSL_Interpolation_1D(PA(NPA:1:-1), FNIS(i,j,NPA:1:-1),&
                                    PAbn(NPA-1:2:-1), BOUNIS(i,j,NPA-1:2:-1), GSLerr)
          ! Do not do the NPA inclusive in the not-a-knot interpolation above -> can lead to negative h,I(NPA)
          BOUNHS(i,j,NPA) = BOUNHS(i,j,NPA-1)
          BOUNIS(i,j,NPA) = BOUNIS(i,j,NPA-1)
          BOUNHS(i,j,1)   = BOUNHS(i,j,2)
          BOUNIS(i,j,1)   = BOUNIS(i,j,2)
       END DO
    END DO
    BOUNHS(:,1,:) = BOUNHS(:,nT,:)
    BOUNIS(:,1,:) = BOUNIS(:,nT,:)
 
    ! Update h and I values for RAM (note that the output of the hI files
    ! now uses RAM variables so the numbers will be different
    DthI = TimeRamElapsed-TOld

    if (integral_smooth) then
       allocate(output(nR,nT))
       do L = 2,NPA
          call gaussian_kernel(1.0, kernel)
          call convolve(FNHS(2:nR+1,:,L), kernel, output)
          FNHS(2:nR+1,:,L) = output
          FNHS(2:nR+1,1,L) = FNHS(2:nR+1,nT,L)

          call gaussian_kernel(1.0, kernel)
          call convolve(FNIS(2:nR+1,:,L), kernel, output)
          FNIS(2:nR+1,:,L) = output
          FNIS(2:nR+1,1,L) = FNIS(2:nR+1,nT,L)

          call gaussian_kernel(1.0, kernel)
          call convolve(BOUNHS(2:nR+1,:,L), kernel, output)
          BOUNHS(2:nR+1,:,L) = output
          BOUNHS(2:nR+1,1,L) = BOUNHS(2:nR+1,nT,L)

          call gaussian_kernel(1.0, kernel)
          call convolve(BOUNIS(2:nR+1,:,L), kernel, output)
          BOUNIS(2:nR+1,:,L) = output
          BOUNIS(2:nR+1,1,L) = BOUNIS(2:nR+1,nT,L)
       ENDDO
       deallocate(output)
    endif

    DO I=2,NR+1
       DO J=1,NT
          DO L=1,NPA
             if (abs(DthI).le.1e-9) then
                dIdt(I,J,L)   = 0._dp
                dIbndt(I,J,L) = 0._dp
             else
                dIdt(I,J,L)   = (FNIS(I,J,L) - FNISPrev(I,J,L))/DThI
                dIbndt(I,J,L) = (BOUNIS(I,J,L) - BOUNISPrev(I,J,L))/DThI
             endif
          ENDDO
          BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
          if (abs(DthI).le.1e-9) then
             dBdt(I,J) = 0._dp
          else
             dBdt(I,J) = (BNES(I,J) - BNESPrev(I,J))/DThI
          endif
       ENDDO
    ENDDO
    DO J=1,NT ! use dipole B at I=1
       BNES(1,J) = 0.32/LZ(1)**3/1.e4
       dBdt(1,J) = 0._dp
       EIR(1,J)  = 0._dp
       EIP(1,J)  = 0._dp
       DO L=1,NPA
          FNHS(1,J,L)   = FNHS(2,J,L)
          FNIS(1,J,L)   = FNIS(2,J,L)
          BOUNHS(1,J,L) = BOUNHS(2,J,L)
          BOUNIS(1,J,L) = BOUNIS(2,J,L)
          HDNS(1,J,L)   = HDNS(2,J,L)
          ODNS(1,J,L)   = ODNS(2,J,L)
          NDNS(1,J,L)   = NDNS(2,J,L)
          dIdt(1,J,L)   = 0._dp
          dIbndt(1,J,L) = 0._dp
       ENDDO
    ENDDO
    TOld = TimeRamElapsed

    ! Check for NaN's
    do i = 2, nR+1
       do j = 1, nT
          do L = 1, nPa
             if (isnan(FNIS(I,J,L))) FNIS(I,J,L) = FNIS(I-1,J,L)
             if (isnan(FNHS(I,J,L))) FNHS(I,J,L) = FNHS(I-1,J,L)
             if (isnan(BOUNIS(I,J,L))) BOUNIS(I,J,L) = BOUNIS(I-1,J,L)
             if (isnan(BOUNHS(I,J,L))) BOUNHS(I,J,L) = BOUNHS(I-1,J,L)
             if (isnan(HDNS(I,J,L))) HDNS(I,J,L) = HDNS(I-1,J,L)
             if (isnan(dIdt(I,J,L))) dIdt(I,J,L) = 0._dp
             if (isnan(dIbndt(I,J,L))) dIbndt(I,J,L) = 0._dp
          enddo
       enddo
    enddo

    DEALLOCATE(HDen, ODen, NDen, distance)
    DEALLOCATE(bfMirror, yI, yH, yD)
    DEALLOCATE(I_cart, h_cart, bZEq_Cart, HDens_Cart, ODens_Cart, NDens_cart)
    DEALLOCATE(outsideSCB, BNESPrev, FNHSPrev, FNISPrev, BOUNISPrev, BOUNHSPrev)

    RETURN
  
  END SUBROUTINE computehI_test

END MODULE ModRamScb
