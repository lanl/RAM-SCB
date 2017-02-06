!******************************************************************************
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

SUBROUTINE Interpolation_Delaunay_Euler_2D(r_local, azim_local, alpha_l, beta_l, alphaCyl_l, betaCyl_l)

  ! gives values on polar grid, using interpolation

  ! .. Use Statements ..
  USE nrtype, ONLY : DP, pi_d, twopi_d
  USE Module1, ONLY : x, y, z, nThetaEquator, nZetaMidnight, rank, numProc, alphaVal
  USE mpi
  USE ModIoUnit, ONLY: UNITTMP_
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
  REAL(DP), ALLOCATABLE :: xTriangle1(:), xTriangle2(:), xTriangle3(:), &
       yTriangle1(:), yTriangle2(:), yTriangle3(:)

  REAL(DP), ALLOCATABLE :: centerVector(:,:), radSqVector(:)
  REAL(DP), ALLOCATABLE :: dSq(:)

  REAL(DP) :: A(2), B(2), C(2), D(2), coordPoint(2)
  REAL(DP) :: baryCoord(3)

  REAL(DP), PARAMETER :: tiny = 1.E-2_dp

  INTEGER :: jy, it, isInTriangle,iDim, itMin, inBigSphere

  LOGICAL :: degen

  npsil = SIZE(alpha_l)
  nzetal = SIZE(beta_l)

  nxl = SIZE(alphaCyl_l,1)
  nyl = SIZE(alphaCyl_l,2)

  !C PRINT*, 'nxl, nyl = ', nxl, nyl

  nTotal = npsil*(nzetal) + 1

  IF (.NOT. ALLOCATED(iTriangle)) ALLOCATE(iTriangle(SIZE(r_local), SIZE(azim_local)), STAT = ierr)

  IF (.NOT. ALLOCATED(xScatter)) ALLOCATE(xScatter(nTotal), stat = ierr)
  IF (.NOT. ALLOCATED(yScatter)) ALLOCATE(yScatter(nTotal), stat = ierr)
  IF (.NOT. ALLOCATED(radSqScatter)) ALLOCATE(radSqScatter(nTotal), stat = ierr)

  IF (.NOT. ALLOCATED(alphaScatter)) ALLOCATE(alphaScatter(nTotal), stat = ierr)
  IF (.NOT. ALLOCATED(betaScatter)) ALLOCATE(betaScatter(nTotal), stat = ierr)  

  iTotal = 0

  j = 1
  j_loop: DO WHILE (j <= npsil)     
     k = 2
     k_Loop:  DO WHILE(k < nzetal+2) 
        iTotal = iTotal+1
        radius = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
        angle = ASIN(y(nThetaEquator,j,k) / radius) !C + pi_d
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
             angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
             angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
        IF ((x(nThetaEquator,j,k) > 0) .AND. (y(nThetaEquator,j,k) <= 0)) &
             angle = twopi_d + ASIN(y(nThetaEquator,j,k) / radius)
        xScatter(iTotal) = angle !C alphaVal(k) !C x(nThetaEquator,j,k)
        yScatter(iTotal) = SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) !C y(nThetaEquator,j,k)
        radSqScatter(iTotal) = xScatter(iTotal)**2 + yScatter(iTotal)**2

        alphaScatter(iTotal) = alpha_l(j)        
        betaScatter(iTotal) = beta_l(k-1)
        k = k+1
     END DO k_Loop
     j = j+1
  END DO j_loop

  ! Add a point in the center for proper triangulation
  ! iTotal = iTotal + 1
  ! xScatter(iTotal) = 0._dp
  ! yScatter(iTotal) = 0._dp
  ! radSqScatter(iTotal) = 0._dp
  ! alphaScatter(iTotal) = alpha_l(1,nZetaMidnight) 
  ! betaScatter(iTotal) = beta_l(1,nZetaMidnight) 

  ! Output points to file for Delaunay triangulation
  OPEN(UNITTMP_, file='/tmp/Delaunay_input', action='write', status='replace')
  WRITE(UNITTMP_, '(A, I6)') '2  Produced by 3D code.    '
  WRITE(UNITTMP_,'(I6)') iTotal
  WRITE(UNITTMP_, '(F7.3, 1X, F7.3)') (xScatter(i), yScatter(i), i = 1, iTotal)
  CLOSE(UNITTMP_)

  ! Delaunay triangulation uses qdelaunay external program 

  OPEN(UNITTMP_, file='/tmp/Delaunay_output', action='read', status='old')
  READ(UNITTMP_, '(I6)') numberTriangles
  ALLOCATE(index1(numberTriangles), index2(numberTriangles), index3(numberTriangles))
  READ(UNITTMP_, *) (index1(i), index2(i), index3(i), i=1,numberTriangles)
  CLOSE(UNITTMP_)

  ! The changes below because the output from qdelaunay is C-style, i.e. vectors from index 0 to n-1 !!!
  index1 = index1 + 1
  index2 = index2 + 1
  index3 = index3 + 1

  ALLOCATE(xTriangle1(numberTriangles),xTriangle2(numberTriangles),xTriangle3(numberTriangles),&
       yTriangle1(numberTriangles),yTriangle2(numberTriangles),yTriangle3(numberTriangles))

  ALLOCATE(dSq(numberTriangles), STAT = ierr)

  DO i = 1, numberTriangles
     ! Coordinates of triangle i
     xTriangle1(i) = xScatter(index1(i))
     xTriangle2(i) = xScatter(index2(i))
     xTriangle3(i) = xScatter(index3(i))

     yTriangle1(i) = yScatter(index1(i))
     yTriangle2(i) = yScatter(index2(i))
     yTriangle3(i) = yScatter(index3(i))
  END DO

  ! Search procedure to find out in which triangle each point of the Cartesian grid is
  X_loop:  DO ix = 1, nxl
     Y_loop:     DO jy = 1, nyl
        !C print*, 'ix, jy = ', ix, jy
        iTriangle(ix,jy) = 0     
        Triangles_loop: DO it = 1, numberTriangles
           A(1) = xTriangle1(it)
           A(2) = yTriangle1(it)
           B(1) = xTriangle2(it)
           B(2) = yTriangle2(it)
           C(1) = xTriangle3(it)
           C(2) = yTriangle3(it)
           triangle(1,1) = A(1)
           triangle(2,1) = A(2)
           triangle(2,2) = B(1)
           triangle(2,2) = B(2)
           triangle(1,3) = C(1)
           triangle(2,3) = C(2)
           coordPoint(1) = azim_local(jy) !C r_local(ix) * COS(azim_local(jy))
           coordPoint(2) = r_local(ix) !C r_local(ix) * SIN(azim_local(jy))

           ! Find the barycentric coordinates of point wrt triangle it
           CALL barycoord_mine_2D(a, b, c, coordPoint, baryCoord, degen)

           IF(baryCoord(1) > -tiny .AND. baryCoord(2) > -tiny .AND. baryCoord(3) > -tiny) THEN
              ! Point inside triangle
              iTriangle(ix,jy) = it  
              ! linear interpolation to find value of function at the respective point
              alphaCyl_l(ix,jy) = baryCoord(1)*alphaScatter(index1(it)) + baryCoord(2)*alphaScatter(index2(it)) + &
                   baryCoord(3)*alphaScatter(index3(it)) 
              betaCyl_l(ix,jy) = baryCoord(1)*betaScatter(index1(it)) + baryCoord(2)*betaScatter(index2(it)) + &
                   baryCoord(3)*betaScatter(index3(it)) 
              EXIT Triangles_loop
           END IF
        END DO Triangles_loop

        IF (iTriangle(ix, jy) == 0) THEN 
           ! point is not inside any triangle

           dSq = 0.0_dp
           coordPoint(1) = azim_local(jy) !C r_local(ix) * COS(azim_local(jy))
           coordPoint(2) = r_local(ix) !C r_local(ix) * SIN(azim_local(jy))
           radSqLoc = coordPoint(1)**2 + coordPoint(2)**2

           Triangles_loop_2: DO it = 1, numberTriangles
              A(1) = xTriangle1(it)
              A(2) = yTriangle1(it)
              B(1) = xTriangle2(it)
              B(2) = yTriangle2(it)
              C(1) = xTriangle3(it)
              C(2) = yTriangle3(it)
              ! Try nearest point instead of nearest triangle center
              dSq(it) = (coordPoint(1)-A(1))**2 + (coordPoint(2)-A(2))**2 
           END DO Triangles_loop_2

           ! Index of triangle closest to point
           itMin = SUM(MINLOC(dSq))

           ! Here extrapolate the ratio B**2/Bdip**2, and then find alpha
           alphaCyl_l(ix,jy) = (radSqScatter(index1(itMin))/radSqLoc) * alphaScatter(index1(itMin)) 
           ! Beta (equiv. to azimuthal angle) found by simple nearest neighbor 
           betaCyl_l(ix,jy) = betaScatter(index1(itMin)) 

2002       CONTINUE
        END IF
     END DO Y_loop
  END DO X_loop


  DEALLOCATE(iTriangle, STAT = ierr)
  DEALLOCATE(index1, index2, index3)
  DEALLOCATE(xTriangle1,xTriangle2,xTriangle3,yTriangle1,yTriangle2,yTriangle3)

  DEALLOCATE(dSq, STAT = ierr)

  DEALLOCATE(xScatter, stat = ierr)
  DEALLOCATE(yScatter, stat = ierr)
  DEALLOCATE(radSqScatter, stat = ierr)

  DEALLOCATE(alphaScatter, stat = ierr)
  DEALLOCATE(betaScatter, stat = ierr)


  RETURN

END SUBROUTINE Interpolation_Delaunay_Euler_2D


