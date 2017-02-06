SUBROUTINE mapthe1d(const_theta)
  USE nrtype, ONLY : SP, DP, PI_D
  USE Module_points
  use ModRamFunctions
  USE ModRamCouple, ONLY : nPoints, nRadSWMF, nLonSWMF, &
       xSWMF, ySWMF, zSWMF, LatSWMF, LonSWMF

  IMPLICIT NONE

  INTERFACE splint_interface
     FUNCTION splint(xa,ya,y2a,x)
       USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
       USE nr, ONLY: locate
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: splint
       INTEGER(I4B) :: khi,klo,n
       REAL(DP) :: a,b,h
     END FUNCTION splint
  END INTERFACE

  INTERFACE spline_interface
     SUBROUTINE spline(x,y,yp1,ypn,y2)
       USE nrtype; USE nrutil, ONLY : assert_eq
       USE nr, ONLY : tridag
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(DP), INTENT(IN) :: yp1,ypn
       REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
       INTEGER(I4B) :: n
       REAL(DP), DIMENSION(SIZE(x)) :: a,b,c,r
     END SUBROUTINE spline
  END INTERFACE

  REAL(DP), INTENT(IN)   :: const_theta
  REAL(DP), DIMENSION(:), ALLOCATABLE :: distance, xOld, yOld, zOld, distance2derivsX, &
       & distance2derivsY, distance2derivsZ
  REAL(DP), DIMENSION(:), ALLOCATABLE :: theval, thetaVal, chiVal
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lonSWMFEarth
  INTEGER :: ierr, i, ii, j, k
  LOGICAL, EXTERNAL :: isNaN

  !   move theta coordinates along one line

  ALLOCATE(distance(nthe))
  ALLOCATE(xOld(nthe))
  ALLOCATE(yOld(nthe))
  ALLOCATE(zOld(nthe))
  ALLOCATE(distance2derivsX(nthe))
  ALLOCATE(distance2derivsY(nthe))
  ALLOCATE(distance2derivsZ(nthe))
  ALLOCATE(theval(nthe))
  ALLOCATE(thetaVal(nthe))
  ALLOCATE(chiVal(nthe))

  ierr = 0

  ! print*, 'mapthe1d: nthe = ', nthe

  IF (.NOT. ALLOCATED(LatSWMF)) ALLOCATE(LatSWMF(nRadSWMF,nLonSWMF), stat = ierr)
  IF (.NOT. ALLOCATED(LonSWMF)) ALLOCATE(LonSWMF(nRadSWMF,nLonSWMF), stat = ierr)
  IF (.NOT. ALLOCATED(LonSWMFEarth)) ALLOCATE(LonSWMFEarth(nRadSWMF,nLonSWMF), stat = ierr)

  DO k = 1, nLonSWMF-1
     DO j = 1, nRadSWMF
        LatSWMF(j,k) = ASIND(zSWMF(nthe,j,k) / &
             SQRT(xSWMF(nthe,j,k)**2+ySWMF(nthe,j,k)**2+&
             zSWMF(nthe,j,k)**2))
        LonSWMF(j,k) = ATAN2D(ySWMF(nthe,j,k), xSWMF(nthe,j,k))
        LonSWMFEarth(j,k) = ATAN2D(ySWMF(nPoints,j,k), xSWMF(nPoints,j,k))
        IF (LonSWMF(j,k) < 0. .AND. k > nLonSWMF/4) LonSWMF(j,k) = LonSWMF(j,k) + 360.
        IF (LonSWMF(j,k) < 90. .AND. k > 3*nLonSWMF/4) LonSWMF(j,k) = LonSWMF(j,k) + 360.
        if (LonSWMFEarth(j,k) < 0.) LonSWMFEarth(j,k) = LonSWMFEarth(j,k) + 360.
        !if (iDebug) WRITE(*,'(A, 2I4, 3F12.5)') 'mapthe1d: k, j, LatI, LonI, LonE = ', &
        !     k, j, LatSWMF(j,k), LonSWMF(j,k), LonSWMFEarth(j,k)

        !       if (j == 14) print*, 'mapthe1d: j, k, x(1,j,k) = ', j, k, xSWMF(1,j,k)
        distance(1) = 0._dp
        !        PRINT*, 'mapthe1d: j, k, 2, x, y, z: ', j, k, REAL(xSWMF(2,j,k),sp), &
        !             REAL(ySWMF(2,j,k),sp), REAL(zSWMF(2,j,k),sp)
        xOld(:) = xSWMF(:,j,k)
        yOld(:) = ySWMF(:,j,k)
        zOld(:) = zSWMF(:,j,k)

        DO  i = 2, nthe
           distance(i) = distance(i-1) + SQRT((xSWMF(i,j,k)-xSWMF(i-1,j,k))**2 +  &
                (ySWMF(i,j,k)-ySWMF(i-1,j,k))**2 +(zSWMF(i,j,k)-zSWMF(i-1,j,k))**2)
           !C         if (j == 14 .and. k == 4) PRINT*, 'mapthe1d: j, k, i, x, y, z, distance, diff: ', j, k, i, REAL(xOld(i),sp), &
           !C              REAL(yOld(i),sp), REAL(zOld(i),sp), REAL(distance(i),sp), real(distance(i)-distance(i-1),sp)
        END DO

        DO  ii = 1, nthe
           theval(ii) = distance(nthe) * REAL(ii - 1,dp) / REAL(nthe - 1,dp)
        END DO

        thetaVal = theval * pi_d / distance(nthe)
        chiVal = (thetaVal + const_theta * SIN(2.*thetaVal)) * distance(nthe)/pi_d

        ! define chival here

        !C        PRINT*, 'mapthe1d: j, k, x(1), y(1), z(1) : ', j, k, xSWMF(1,j,k), ySWMF(1,j,k), zSWMF(1,j,k)

        ! "Natural" splines
        CALL spline(distance(1:nthe), xOld(1:nthe), 1.E31_DP, 1.E31_DP, distance2derivsX(1:nthe))
        CALL spline(distance(1:nthe), yOld(1:nthe), 1.E31_DP, 1.E31_DP, distance2derivsY(1:nthe))
        CALL spline(distance(1:nthe), zOld(1:nthe), 1.E31_DP, 1.E31_DP, distance2derivsZ(1:nthe))

        DO i = 2, nthe-1
           xSWMF(i,j,k) = splint(distance(1:nthe), xOld(1:nthe), distance2derivsX(1:nthe), chiVal(i))
           ySWMF(i,j,k) = splint(distance(1:nthe), yOld(1:nthe), distance2derivsY(1:nthe), chiVal(i))
           zSWMF(i,j,k) = splint(distance(1:nthe), zOld(1:nthe), distance2derivsZ(1:nthe), chiVal(i))
        END DO


     END DO
     !C if (k == 4) print*, 'mapthe1d: xEq(:,k) = ', xSWMF(nPoints,:,k)
  END DO

  DEALLOCATE(distance, xOld, yOld, zOld, distance2derivsX, distance2derivsY, distance2derivsZ, &
       theval, thetaVal, chiVal)

  RETURN

END SUBROUTINE mapthe1d



