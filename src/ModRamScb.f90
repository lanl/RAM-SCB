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
  
    use ModRamGrids, ONLY: NS, NE, NPA
    use ModScbGrids, ONLY: nthe, npsi, nzeta
  
    implicit none
  
    ALLOCATE(INDEXPA(nthe,npsi,nzeta,NPA), Flux3DEq(NS,npsi,nzeta,NE,NPA))
    INDEXPA = 0.0; Flux3DEq = 0.0
  
    return
  
  end subroutine ramscb_allocate

!==================================================================================================
  subroutine ramscb_deallocate
  
    implicit none
  
    DEALLOCATE(INDEXPA, Flux3DEq)
  
  end subroutine ramscb_deallocate

!==================================================================================================
  subroutine Compute3DFlux
  
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
  
    use ModRamVariables, ONLY: gdLon, gdLat, FNHS, FNIS, BNES, HDNS, dBdt, dHdt, &
                               dIdt, dIbndt, BOUNHS, BOUNIS, EIR, EIP, flux_volume, &
                               LZ, MU, MLT, PAbn, PA, DL1, outsideMGNP
    use ModRamTiming,    ONLY: TimeRamElapsed, TOld
    use ModRamConst,     ONLY: b0dip
    use ModRamParams,    ONLY: NameBoundMag, verbose, checkMGNP, integral_smooth
    use ModRamGrids,     ONLY: nR, nT, nPa, radiusMax, radiusMin
  
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbParams,    ONLY: method, constTheta
    use ModScbVariables, ONLY: bf, chiVal, x, y, z, bnormal, h_Cart, I_Cart, h_Cart_interp, &
                               i_Cart_interp, bZEq_Cart, flux_vol_cart, bZ, &
                               hdens_Cart, radRaw, azimRaw, nThetaEquator, fluxVolume, &
                               thetaVal, psi, alfa
 
    use ModRamCouple,    ONLY: BLines_DIII, nPoints, IsClosed_II

    use ModRamGSL,       ONLY: GSL_Interpolation_2D, GSL_Interpolation_1D, &
                               GSL_Integration_hI, GSL_Smooth_1D
    use ModRamFunctions, ONLY: RamFileName, FUNT, FUNI
    use ModScbFunctions, ONLY: extap, locate, get_dipole_lines
    use ModScbIO,        ONLY: trace, PARMOD, IOPT
    use gaussian_filter, only: gaussian_kernel, convolve
  
    use nrtype,    ONLY: DP, pi_d, twopi_d
  
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
                             bfmirror(:,:,:)
  
    ! Variables for RAM
    integer  :: nEquator
    integer, ALLOCATABLE :: ScaleAt(:), outsideSCB(:,:)
    REAL(DP) :: t0, t1, rt, tt, zt
    REAL(DP), ALLOCATABLE :: bbx(:), bby(:), bbz(:), xx(:), yy(:), zz(:), cval(:)
    REAL(DP), ALLOCATABLE :: BNESPrev(:,:), FNISPrev(:,:,:), BOUNISPrev(:,:,:), &
                             FNHSPrev(:,:,:), BOUNHSPrev(:,:,:), dHbndt(:,:,:)
    REAL(DP), ALLOCATABLE :: xRAM(:,:,:), yRAM(:,:,:), zRAM(:,:,:), bRAM(:,:,:)
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

    clock_rate = 1000
    clock_max = 100000
    LMAX = 1000
  
    !!! Initialize Allocatable Arrays
    ALLOCATE(ScaleAt(nT), outsideSCB(nR,nT), bbx(nthe), bby(nthe), bbz(nthe), xx(nthe), &
             yy(nthe), zz(nthe), cval(nthe), distance(nthe,nR,nT), length(nR,nT), &
             r0(nR,nT), bfMirror(nR,nT,nPa), yI(nR,nT,NPA), yH(nR,nT,NPA), yD(nR,nT,NPA), &
             BNESPrev(nR+1,nT), FNHSPrev(nR+1,nT,nPa),FNISPrev(nR+1,nT,nPa), &
             BOUNISPrev(nR+1,nT,nPa), BOUNHSPrev(nR+1,nT,nPa), dHbndt(nR+1,nT,nPa), &
             xRAM(nthe,nR,nT),yRAM(nthe,nR,nT),zRAM(nthe,nR,nT),bRAM(nthe,nR,nT))
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
       write(*,'(1x,a)',ADVANCE='NO') 'Converting SCB field lines to RAM field lines'
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)

       ! If using SWMF without SCB calculation just take directly from grids
       if ((NameBoundMag.eq.'SWMF').and.(method == 3)) then
          xRAM(:,:,:) = x(:,1:nR,:)
          yRAM(:,:,:) = y(:,1:nR,:)
          zRAM(:,:,:) = z(:,1:nR,:)
          bRAM(:,:,:) = bf(:,1:nR,:)
       endif
   
  !$OMP PARALLEL DO
       do i = 1,nR
          do j = 2,nT
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
       write(*,'(a,1x,F6.2,1x,a)') ': Completed in', stoptime-starttime, 'seconds'

       ! If a point on the RAM grid lies outside the SCB grid then we need to
       ! calculate the scaling parameters (if still inside the magnetopause) or
       ! set the outsideMGNP flag to 1 so we can track it in RAM
       do j = 2,nT
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

    !!!!! BEGIN INTEGRAl CALCULATION
    write(*,'(1x,a)',ADVANCE='NO') 'Calculating h and I integrals'
    ! Start timing
    call system_clock(time1,clock_rate,clock_max)
    starttime=time1/real(clock_rate,dp)

  !$OMP PARALLEL DO
    do i = 1, nR
       do j = 2, nT
          if (outsideMGNP(i,j) == 0) then
             distance(:,i,j) = SQRT(xRAM(:,i,j)**2+yRAM(:,i,j)**2+zRAM(:,i,j)**2) ! Distance from center of earth
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

             CALL GSL_Integration_hI(bfMirror(i,j,:),chiVal(:),bRAM(:,i,j),distance(:,i,j),&
                                     yI(i,j,:),yH(i,j,:),yD(i,j,:))

             I_cart(i,j,:) = (length(i,j)/(pi_d*r0(i,j))) *yI(i,j,:)/SQRT(Bfmirror(i,j,:))
             H_cart(i,j,:) = (length(i,j)/(pi_d*2*r0(i,j))) *yH(i,j,:)*SQRT(Bfmirror(i,j,:))
             HDens_cart(i,j,:) = 1.E5_dp * yD(i,j,:)/yH(i,j,:) !Re-normalize
             bZEq_Cart(i,j) = bRAM(nThetaEquator,i,j)*bnormal

          end if
       enddo
    enddo
  !$OMP END PARALLEL DO

    ! End Timing
    call system_clock(time1,clock_rate,clock_max)
    stoptime=time1/real(clock_rate,dp)
    write(*,'(a,1x,F6.2,1x,a)') ': Completed in', stoptime-starttime, 'seconds'
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
    DEALLOCATE(xRAM,yRAM,zRAM,bRAM,outsideSCB)
 
    RETURN
  
  END SUBROUTINE computehI

END MODULE ModRamScb
