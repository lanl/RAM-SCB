MODULE ModScbInterp
  ! Contains subroutines used for interpolation in the SCB code
  
  implicit none
  
  contains

!============================================================================== 
  SUBROUTINE Interpolation_natgrid_2D_EField(r_local, azim_local, Epot_l, EpotCart_l)

    ! gives values on polar grid, using interpolation
    ! .. Use Statements ..
    USE nrtype,          ONLY: DP, pi_d
    USE ModScbVariables, ONLY: x, y, z, nThetaEquator
    USE ModScbMain,      ONLY: prefixOut

    use ModRamGrids, ONLY: NPA

    IMPLICIT NONE

    integer, parameter :: NN = 9

    REAL(DP), INTENT(IN) :: r_local(:), azim_local(:)
    REAL(DP), INTENT(IN) :: Epot_l(:,:)
    REAL(DP), INTENT(OUT) :: EpotCart_l(:,:)
    INTEGER :: nxl, nyl, npsil, nzetal, npal

    INTEGER :: iTriangle(SIZE(r_local), SIZE(azim_local))
    INTEGER :: nTotal

    REAL(DP) :: xMax, yMax, zMax, xMin, yMin, zMin, dx, dy, dz, bet, del, end_time
    INTEGER :: ix, ier, i, j, k, iTotal, iTotal2, ierr, idealerr
    INTEGER :: nq, nr, nw

    REAL(DP) :: eps, eq, eqx, eqy, eqz, p2y, p2z, q1, rmax, xtemp, pMin
    REAL(DP) :: triangle(2,3)
    INTEGER :: l

    ! .. Intrinsic Functions ..
    INTRINSIC KIND, REAL, SUM
    ! .. Parameters ..
    INTEGER, PARAMETER :: wp = KIND(1.0D0)
    ! .. Local Scalars ..
    INTEGER :: m, numberTriangles
    REAL (DP) :: q, qx, qy, qz, u, v, w, rad, radSqLoc
    REAL(DP) :: center(2)

    REAL(DP), ALLOCATABLE :: xScatter(:), yScatter(:), radSqScatter(:), EpotScatter(:)
    REAL(DP) :: coordPoint(2)

    INTEGER :: jy, it, isInTriangle,iDim, itMin, inBigSphere

!
    REAL(DP), ALLOCATABLE :: distance(:)
    REAL(DP), DIMENSION(NN) :: xNear, yNear, bZNear, fluxVolNear, EPotNear
    REAL(DP), DIMENSION(NN,NPA) :: hNear, INear, hDensNear
    INTEGER, DIMENSION(1) :: iTemp
!
    npsil  = SIZE(Epot_l,1)
    nzetal = SIZE(Epot_l,2) + 1 ! Since only 2:nzeta is being passed
    nxl    = SIZE(EpotCart_l,1)
    nyl    = SIZE(EpotCart_l,2)

    nTotal = npsil*(nzetal-1)

    ALLOCATE(distance(nTotal), stat=ierr)
    ALLOCATE(xScatter(nTotal), stat = ierr)
    ALLOCATE(yScatter(nTotal), stat = ierr)
    ALLOCATE(EpotScatter(nTotal), stat=ierr)

    iTotal = 0
    j = 1
    j_loop: DO WHILE (j <= npsil)
       k = 2
       k_Loop:  DO WHILE(k <= nzetal)
          iTotal = iTotal+1
          xScatter(iTotal) = x(nThetaEquator,j,k)
          yScatter(iTotal) = y(nThetaEquator,j,k)
          EPotScatter(iTotal) = Epot_l(j,k-1)
          k = k+1
       END DO k_Loop
       j = j+1
    END DO j_loop

    DO ix = 1, nxl
       DO jy = 1, nyl
          coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
          coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
          distance = (xScatter - coordPoint(1))**2 + (yScatter - coordPoint(2))**2
          do i = 1,NN
             iTemp = minloc(distance)
             xNear(i)       = xScatter(iTemp(1))
             yNear(i)       = yScatter(iTemp(1))
             EPotNear(i) = EPotScatter(iTemp(1))
             distance(iTemp(1)) = 999999.9
          end do
          CALL DSPNT2D(NN,xNear,yNear,EPotNear,1,coordPoint(1),coordPoint(2),EPotCart_l(ix,jy),ierr)
       END DO
    END DO

    DEALLOCATE(EPotScatter, xScatter, yScatter, distance)
    return

  END SUBROUTINE Interpolation_natgrid_2D_EField

!==============================================================================
  SUBROUTINE Interpolation_natgrid_2D(r_local, azim_local, hFlux_l, IFlux_l, bZEqFlux_l, &
       flux_vol_l, hdensFlux_l, hCart_l, ICart_l, bZEqCart_l, flux_vol_Cart_l, hdensCart_l)
  
    ! gives values on polar grid, using interpolation
  
    ! .. Use Statements ..
    USE nrtype,          ONLY: DP, pi_d
    USE ModScbVariables, ONLY: x, y, z, nThetaEquator
    USE ModScbMain,      ONLY: prefixOut

    use ModRamGrids, ONLY: NPA
 
    IMPLICIT NONE
  
    integer, parameter :: NN = 9

    REAL(DP), INTENT(IN) :: r_local(:), azim_local(:)
    REAL(DP), INTENT(IN) :: hFlux_l(:,:,:), IFlux_l(:,:,:), bZEqFlux_l(:,:), &
         flux_vol_l(:,:), hdensFlux_l(:,:,:)
    REAL(DP), INTENT(OUT) :: hCart_l(:,:,:), ICart_l(:,:,:), bZEqCart_l(:,:), &
         hdensCart_l(:,:,:), flux_vol_Cart_l(:,:)
    INTEGER :: nxl, nyl, npsil, nzetal, npal
  
    INTEGER :: iTriangle(SIZE(r_local), SIZE(azim_local))
    INTEGER :: nTotal 
  
    REAL(DP) :: xMax, yMax, zMax, xMin, yMin, zMin, dx, dy, dz, bet, del, end_time
    INTEGER :: ix, ier, i, j, k, iTotal, iTotal2, ierr, idealerr
    INTEGER :: nq, nr, nw
  
    REAL(DP) :: eps, eq, eqx, eqy, eqz, p2y, p2z, q1, rmax, xtemp, pMin
    REAL(DP) :: triangle(2,3)
    INTEGER :: l
  
    ! .. Intrinsic Functions ..
    INTRINSIC KIND, REAL, SUM
    ! .. Parameters ..
    INTEGER, PARAMETER :: wp = KIND(1.0D0)
    ! .. Local Scalars ..
    INTEGER :: m, numberTriangles
    REAL (DP) :: q, qx, qy, qz, u, v, w, rad, radSqLoc
    REAL(DP) :: center(2)
  
    REAL(DP), ALLOCATABLE :: xScatter(:), yScatter(:), radSqScatter(:), bZScatter(:)
    REAL(DP), ALLOCATABLE :: hScatter(:,:), IScatter(:,:), hdensScatter(:,:), fluxVolScatter(:)
    REAL(DP) :: coordPoint(2)
  
    INTEGER :: jy, it, isInTriangle,iDim, itMin, inBigSphere
 
!
    REAL(DP), ALLOCATABLE :: distance(:)
    REAL(DP), DIMENSION(NN) :: xNear, yNear, bZNear, fluxVolNear, EPotNear
    REAL(DP), DIMENSION(NN,NPA) :: hNear, INear, hDensNear
    INTEGER, DIMENSION(1) :: iTemp
!
    npsil = SIZE(hFlux_l,1)
    nzetal = SIZE(hFlux_l,2) + 1 ! Since only 2:nzeta is being passed
    npal = SIZE(hFlux_l,3)
    nxl = SIZE(hCart_l,1)
    nyl = SIZE(hCart_l,2)
  
    nTotal = npsil*(nzetal-1)

    if (.not.ALLOCATED(distance)) ALLOCATE(distance(nTotal), stat=ierr)

    IF(.NOT.ALLOCATED(xScatter))   ALLOCATE(xScatter(nTotal), stat = ierr)
    IF(.NOT.ALLOCATED(yScatter))   ALLOCATE(yScatter(nTotal), stat = ierr)

    IF(.NOT.ALLOCATED(radSqScatter))   ALLOCATE(radSqScatter(nTotal), stat = ierr)  
    IF(.NOT.ALLOCATED(hScatter))   ALLOCATE(hScatter(nTotal, NPAL), stat = ierr)
    IF(.NOT.ALLOCATED(hDensScatter))   ALLOCATE(hDensScatter(nTotal, NPAL), stat = ierr)
    IF(.NOT.ALLOCATED(IScatter))   ALLOCATE(IScatter(nTotal, NPAL), stat = ierr)
    IF(.NOT.ALLOCATED(bZScatter))   ALLOCATE(bZScatter(nTotal), stat = ierr)  
    IF (.NOT. ALLOCATED(fluxVolScatter)) ALLOCATE(fluxVolScatter(nTotal), stat = ierr)
  
    iTotal = 0
    j = 1
    j_loop: DO WHILE (j <= npsil)     
       k = 2
       k_Loop:  DO WHILE(k <= nzetal) 
          iTotal = iTotal+1
          xScatter(iTotal) = x(nThetaEquator,j,k)
          yScatter(iTotal) = y(nThetaEquator,j,k)
          radSqScatter(iTotal) = xScatter(iTotal)**2 + yScatter(iTotal)**2
          bZScatter(iTotal) = bZEqFlux_l(j,k-1)
          fluxVolScatter(iTotal) = flux_vol_l(j,k-1)
          Pitch_angle_loop: DO L = 1, NPAL
             hScatter(iTotal, L) = hFlux_l(j,k-1, L)
             IScatter(iTotal, L) = IFlux_l(j,k-1, L)
             hdensScatter(iTotal, L) = hdensFlux_l(j,k-1, L)
          END DO Pitch_angle_loop
          k = k+1
       END DO k_Loop
       j = j+1
    END DO j_loop

    DO ix = 1, nxl
       DO jy = 1, nyl
          coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
          coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
!!!!!!!!!!!!!!
! ME My basic attempt at a nearest neighbor search
          distance = (xScatter - coordPoint(1))**2 &
                    +(yScatter - coordPoint(2))**2
          do i = 1,NN
             iTemp = minloc(distance)
             xNear(i)       = xScatter(iTemp(1))
             yNear(i)       = yScatter(iTemp(1))
             bZNear(i)      = bZScatter(iTemp(1))
             fluxVolNear(i) = fluxVolScatter(iTemp(1))
             hNear(i,:)     = hScatter(iTemp(1),:)
             INear(i,:)     = IScatter(iTemp(1),:)
             hDensNear(i,:) = hDensScatter(iTemp(1),:)
             distance(iTemp(1)) = 999999.9
          end do
          CALL DSPNT2D(NN,xNear,yNear,bZNear,1,coordPoint(1),coordPoint(2),bZEqCart_l(ix,jy),ierr)
          CALL DSPNT2D(NN,xNear,yNear,fluxVolNear,1,coordPoint(1),coordPoint(2),flux_vol_Cart_l(ix,jy),ierr)
          do L = 1,NPAL
             CALL DSPNT2D(NN,xNear,yNear,hNear(:,L),1,coordPoint(1),coordPoint(2),hCart_l(ix,jy,L),ierr)
             CALL DSPNT2D(NN,xNear,yNear,INear(:,L),1,coordPoint(1),coordPoint(2),ICart_l(ix,jy,L),ierr)
             CALL DSPNT2D(NN,xNear,yNear,hDensNear(:,L),1,coordPoint(1),coordPoint(2),hDensCart_l(ix,jy,L),ierr)
          end do
!!!!!!!!!!!!!!
       END DO
    END DO
  
    if (ALLOCATED(distance)) DEALLOCATE(distance, stat=ierr)
  
    IF(ALLOCATED(xScatter))  DEALLOCATE(xScatter, stat = ierr)
    IF (ALLOCATED(yScatter))  DEALLOCATE(yScatter, stat = ierr)
    IF(ALLOCATED(radSqScatter))  DEALLOCATE(radSqScatter, stat = ierr)
  
    IF(ALLOCATED(hScatter))  DEALLOCATE(hScatter, stat = ierr)
    IF(ALLOCATED(IScatter))  DEALLOCATE(IScatter, stat = ierr)
    IF(ALLOCATED(bZScatter))  DEALLOCATE(bZScatter, stat = ierr)
    IF(ALLOCATED(hdensScatter))  DEALLOCATE(hdensScatter, stat = ierr)
  
    IF (ALLOCATED(fluxVolScatter)) DEALLOCATE(fluxVolScatter, stat = ierr)
  
    RETURN
  
  END SUBROUTINE Interpolation_natgrid_2D
  
!==============================================================================
  SUBROUTINE Interpolation_natgrid_Euler_2D(r_local, azim_local, alpha_l, beta_l, alphaCyl_l, betaCyl_l)
  
    ! gives values on polar grid, using interpolation
  
    ! .. Use Statements ..
    USE nrtype, ONLY : DP, pi_d, twopi_d
    USE ModScbVariables, ONLY : x, y, z, nThetaEquator, nZetaMidnight, alphaVal
  
    IMPLICIT NONE
  
    REAL(DP), INTENT(IN) :: r_local(:), azim_local(:)
    REAL(DP), INTENT(IN) :: alpha_l(:), beta_l(:)
    REAL(DP), INTENT(OUT) :: alphaCyl_l(:,:), betaCyl_l(:,:)
    INTEGER :: nxl, nyl, npsil, nzetal, npal
  
    INTEGER, ALLOCATABLE :: iTriangle(:, :)
    INTEGER :: nTotal 
  
    REAL(DP) :: xMax, yMax, zMax, xMin, yMin, zMin, dx, dy, dz, bet, del, end_time
    INTEGER :: ix, ier, i, j, k, iTotal, iTotal2, ierr, idealerr
    INTEGER :: nq, nr, nw
  
    REAL(DP) :: eps, eq, eqx, eqy, eqz, p2y, p2z, q1, rmax, xtemp, pMin, radius, angle
    REAL(DP) :: triangle(2,3)
    INTEGER :: l
  
    REAL, ALLOCATABLE :: rvec(:)
    INTEGER, ALLOCATABLE :: indexOrder(:)
    ! .. Intrinsic Functions ..
    INTRINSIC KIND, REAL, SUM
    ! .. Parameters ..
    INTEGER, PARAMETER :: wp = KIND(1.0D0)
    ! .. Local Scalars ..
    INTEGER :: m, numberTriangles
    REAL (DP) :: q, qx, qy, qz, u, v, w, rad, radSqLoc
    REAL(DP) :: center(2)
  
    INTEGER, ALLOCATABLE :: index1(:), index2(:), index3(:)
    REAL(DP), ALLOCATABLE :: xScatter(:), yScatter(:), radSqScatter(:), alphaScatter(:), betaScatter(:)
  
    INTEGER :: jy, it, isInTriangle,iDim, itMin, inBigSphere
    LOGICAL :: degen
  
    npsil = SIZE(alpha_l)
    nzetal = SIZE(beta_l) 
  
    nxl = SIZE(alphaCyl_l,1)
    nyl = SIZE(alphaCyl_l,2)
  
  !SZ  PRINT*, 'INE: npsil, nzetal, nxl, nyl = ', npsil, nzetal, nxl, nyl
  
    nTotal = 2*npsil*(nzetal+1) + 1
  
    IF (.NOT. ALLOCATED(iTriangle)) ALLOCATE(iTriangle(SIZE(r_local), SIZE(azim_local)), STAT = ierr)
  
    IF (.NOT. ALLOCATED(xScatter)) ALLOCATE(xScatter(nTotal), stat = ierr)
    IF (.NOT. ALLOCATED(yScatter)) ALLOCATE(yScatter(nTotal), stat = ierr)
    IF (.NOT. ALLOCATED(radSqScatter)) ALLOCATE(radSqScatter(nTotal), stat = ierr)
  
    IF (.NOT. ALLOCATED(alphaScatter)) ALLOCATE(alphaScatter(nTotal), stat = ierr)
    IF (.NOT. ALLOCATED(betaScatter)) ALLOCATE(betaScatter(nTotal), stat = ierr)  
  
    iTotal = 0
  
    j_loop: DO j = 1, npsil
       k_Loop:  DO k = 2, nzetal
          iTotal = iTotal+1
          radius = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
          angle = ASIN(y(nThetaEquator,j,k) / radius) 
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) > 0) .AND. (y(nThetaEquator,j,k) <= 0)) &
               angle = twopi_d + ASIN(y(nThetaEquator,j,k) / radius)
  
          !C        if (k == 2 .or. k == 3) print*, 'k, j, angle: ', k, j, real(angle)
          IF (k == 2 .AND. angle > pi_d) then
          !   print*, 'j, k, angle = ', j, k, angle
             angle = angle - twopi_d  ! Enforce small angle at k=2
          !   print*, 'Interp_natgrid_Euler: int enforced 1.'
          end IF
          !C IF (k == 2 .AND. angle < 0) angle = 0.
  
          IF (k == nzetal+1 .AND. angle < pi_d) then
             angle = angle + twopi_d  ! Enforce large angle at k=nzetp
          !   print*, 'Inter_natgrid_Euler: int enforced 2.'
          end IF
  
          yScatter(iTotal) = angle !C x(nThetaEquator,j,k)
  
          xScatter(iTotal) = SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) !C y(nThetaEquator,j,k)
          radSqScatter(iTotal) = xScatter(iTotal)**2 + yScatter(iTotal)**2
  
          alphaScatter(iTotal) = alpha_l(j)        
          betaScatter(iTotal) = beta_l(k-1)
       END DO k_Loop
  
       ! Add duplicate points at angle > 2*pi
       do k = 2, nzetal/4
          iTotal = iTotal+1
          radius = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
          angle = ASIN(y(nThetaEquator,j,k) / radius)  
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) > 0) .AND. (y(nThetaEquator,j,k) <= 0)) &
               angle = twopi_d + ASIN(y(nThetaEquator,j,k) / radius)
  
          if (angle > pi_d) angle = angle - twopi_d
  
          angle = angle + twopi_d
          yScatter(iTotal) = angle !C x(nThetaEquator,j,k)
          xScatter(iTotal) = SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) !C y(nThetaEquator,j,k)
          alphaScatter(iTotal) = alpha_l(j)
          betaScatter(iTotal) = beta_l(k-1) + twopi_d
       end do
  
        ! Add duplicate points at angle < 0
       do k = nzetal+1-nzetal/4, nzetal
          iTotal = iTotal+1
          radius = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
          angle = ASIN(y(nThetaEquator,j,k) / radius)  
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
               angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) > 0) .AND. (y(nThetaEquator,j,k) <= 0)) &
               angle = twopi_d + ASIN(y(nThetaEquator,j,k) / radius)
          angle = angle - twopi_d
          yScatter(iTotal) = angle 
          xScatter(iTotal) = SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) !C y(nThetaEquator,j,k)
          alphaScatter(iTotal) = alpha_l(j)
          betaScatter(iTotal) = beta_l(k-1) - twopi_d
       end do
    END DO j_loop
  
    CALL NNSETI('IGR', 1)  ! Non-linear interpolation
    !C  call NNSETI('EXT', 0) ! No extrapolation outside the convex hull
    !C  call NNSETR('NUL', -10.) ! Value when ext not allowed
  
    CALL NATGRIDD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), alphaScatter(1:iTotal), &
         nxl, nyl, r_local(1:nxl), azim_local(1:nyl), alphaCyl_l(1:nxl,1:nyl), IER)
  
    IF (IER .NE. 0) THEN
       WRITE (6,510) IER
  510  FORMAT('Error return from NATGRIDD = ',I3)
    END IF
  
    !C  CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), alphaScatter(1:iTotal))
    !C DO ix = 1, nxl
    !C   DO jy = 1, nyl
    !C      coordPoint(1) = r_local(ix) * COS(azim_local(jy))
    !C      coordPoint(2) = r_local(ix) * SIN(azim_local(jy))
    !C      !C print*, 'coord: ', coordPoint(1), coordPoint(2)
    !C      CALL NNPNTD(coordPoint(1),coordPoint(2),alphaCyl_l(ix,jy))
    !C   END DO
    !C END DO
    !C CALL NNPNTENDD()
  
    !SZ CALL NNSETI('ADF', 1) ! Outputs triangulation file (default name 'nnalg.dat'), can be seen with nnalg
    !SZ CALL NNSETC('ALG', 'triangulation_'//ST3//'.dat')
  
    CALL NATGRIDD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), betaScatter(1:iTotal), &
         nxl, nyl, r_local(1:nxl), azim_local(1:nyl), betaCyl_l(1:nxl,1:nyl), IER)
    IF (IER .NE. 0) THEN
       WRITE (6,510) IER
    END IF
  
    ! Revert to default values for speed in hRAM
    CALL NNSETI('IGR', 0)
    CALL NNSETI('ADF', 0)
  
  
    !C  CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), betaScatter(1:iTotal))
    !C  DO ix = 1, nxl
    !C     DO jy = 1, nyl
    !C        coordPoint(1) = r_local(ix) * COS(azim_local(jy))
    !C        coordPoint(2) = r_local(ix) * SIN(azim_local(jy))
    !C        CALL NNPNTD(coordPoint(1),coordPoint(2),betaCyl_l(ix,jy))
    !C     END DO
    !C  END DO
    !C  CALL NNPNTENDD() 
  
    IF(ALLOCATED(xScatter))  DEALLOCATE(xScatter, stat = ierr)
    IF (ALLOCATED(yScatter))  DEALLOCATE(yScatter, stat = ierr)
    IF(ALLOCATED(radSqScatter))  DEALLOCATE(radSqScatter, stat = ierr)
  
    IF(ALLOCATED(alphaScatter))  DEALLOCATE(alphaScatter, stat = ierr)
    IF(ALLOCATED(betaScatter))  DEALLOCATE(betaScatter, stat = ierr)
  
  
    RETURN
  
  END SUBROUTINE Interpolation_natgrid_Euler_2D
  
!==============================================================================
  SUBROUTINE Interpolation_natgrid_Shellbuild(quantitySWMF, lat_local, lon_local, iPoint, nlat, nlon, quantityInterp)
  
    ! gives values on polar grid, using interpolation
  
    ! .. Use Statements ..
    USE nrtype, ONLY : SP, DP, pi_d
    USE ModRamCouple, ONLY : nRadSWMF, nLonSWMF, xSWMF, ySWMF, &
                             zSWMF, LatSWMF, LonSWMF
  
    IMPLICIT NONE
  
    REAL(DP), INTENT(IN) :: lat_local(:,:), lon_local(:)
    REAL(DP), INTENT(IN) :: quantitySWMF(:,:)
    REAL(DP), INTENT(OUT) :: quantityInterp(nlat, nlon)
    INTEGER, INTENT(IN) :: iPoint, nlat, nlon
  
    INTEGER :: ix, j, jy, k, iTotal, ierr, nTotal
  
    real(dp) :: LatInner(nLonSWMF)
    REAL(DP), ALLOCATABLE :: latScatter(:), lonScatter(:), quantityScatter(:)
    !-----------------------------------------------------------------------
  
    ! Set total number of points, including extra to provide for continuity 
    ! across 360 to 0 degree boundary and inner boundary ghost cells.
    nTotal = (nRadSWMF+1)*(nLonSWMF+20)
  
    ! Create inner latitude values that map inside of inner boundary:
    LatInner = 2.*LatSWMF(1,:) - LatSWMF(2,:)
  
    ! Allocate variables to hold unordered points:
    IF (.NOT. ALLOCATED(latScatter)) ALLOCATE(latScatter(nTotal), stat = ierr)
    IF (.NOT. ALLOCATED(lonScatter)) ALLOCATE(lonScatter(nTotal), stat = ierr)
    IF (.NOT. ALLOCATED(quantityScatter)) ALLOCATE(quantityScatter(nTotal), &
         stat = ierr)  
  
    ! Begin tracking the current point:
    iTotal = 0
  
    ! Unpack grid points into 1D vector:
    k_Loop:  DO k = 1, nLonSWMF-1 
       ! Add "safety layer" at inner radius:
       iTotal = iTotal + 1
       latScatter(iTotal) = LatInner(k)
       lonScatter(iTotal) = LonSWMF(1,k)
       quantityScatter(iTotal) = quantitySWMF(1,k)
       ! Add real points outside of inner ring:
       j_loop: DO j = 1, nRadSWMF
          iTotal = iTotal+1
          latScatter(iTotal) = LatSWMF(j,k)
          lonScatter(iTotal) = LonSWMF(j,k)
          quantityScatter(iTotal) = quantitySWMF(j,k)
       END DO j_loop
    END DO k_Loop
  
    ! Add another set of points for 1 and nLonSWMF.  Small negative longitude 
    ! will ensure that neighbors are found for 0 degrees longitude if field 
    ! lines are curved, and no extrapolation is needed)
    DO k = nLonSWMF-10, nLonSWMF-1 ! 10 longitudes
        ! Add "safety layer" at inner radius:
       iTotal = iTotal + 1
       latScatter(iTotal) = LatInner(k)
       lonScatter(iTotal) = LonSWMF(1,k)
       quantityScatter(iTotal) = quantitySWMF(1,k)
       DO j = 1, nRadSWMF     
          iTotal = iTotal+1
          latScatter(iTotal) = LatSWMF(j,k)
          lonScatter(iTotal) = LonSWMF(j,k) - 360.
          quantityScatter(iTotal) = quantitySWMF(j,k)
       END DO
    END DO
    DO k = 1, 10 ! 10 longitudes
       ! Add "safety layer" at inner radius:
       iTotal = iTotal + 1
       latScatter(iTotal) = LatInner(k)
       lonScatter(iTotal) = LonSWMF(1,k)
       quantityScatter(iTotal) = quantitySWMF(1,k)
       DO j = 1, nRadSWMF     
          iTotal = iTotal+1
          latScatter(iTotal) = LatSWMF(j,k)
          lonScatter(iTotal) = LonSWMF(j,k) + 360.
          quantityScatter(iTotal) = quantitySWMF(j,k)
       END DO
    END DO
  
    ! Initialize delaunay triangulation:
    CALL NNPNTINITD(iTotal, latScatter, lonScatter, quantityScatter)
  
    ! Interpolate:
    DO ix = 1, nlat
       DO jy = 1, nlon
          CALL NNPNTD(lat_local(ix,jy),lon_local(jy),quantityInterp(ix,jy))
       END DO
    END DO
  
    ! Clean up interpolator:
    CALL NNPNTENDD()
  
    ! Deallocate arrays.
    IF (ALLOCATED(latScatter)) DEALLOCATE(latScatter, stat = ierr)
    IF (ALLOCATED(lonScatter)) DEALLOCATE(lonScatter, stat = ierr)
    IF (ALLOCATED(quantityScatter)) DEALLOCATE(quantityScatter, stat = ierr)
  
    RETURN
  
  
  END SUBROUTINE Interpolation_natgrid_Shellbuild
  
  END MODULE ModScbInterp
