SUBROUTINE Interpolation_natgrid_Euler_2D(r_local, azim_local, alpha_l, beta_l, alphaCyl_l, betaCyl_l)

  ! gives values on polar grid, using interpolation

  ! .. Use Statements ..
  USE nrtype, ONLY : DP, pi_d, twopi_d
  USE Module1, ONLY : x, y, z, nThetaEquator, nZetaMidnight, rank, numProc, alphaVal, iST3
  USE mpi

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
  CHARACTER*4 :: ST3

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

  write(ST3, '(I4.4)') iST3

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


