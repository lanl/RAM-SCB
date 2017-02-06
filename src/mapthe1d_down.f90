SUBROUTINE mapthe1d_down(const_theta)
  USE nrtype, ONLY : SP, DP, PI_D
  USE Module_points
  USE ModRamCouple, ONLY : nPoints, nLonSWMF, &
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
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: distance, xOld, yOld, zOld, distance2derivsX, &
       & distance2derivsY, distance2derivsZ
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: theval, thetaVal, chiVal
  INTEGER :: ierr, i, ii, j, k
  LOGICAL, EXTERNAL :: isNaN

  ALLOCATE(distance(ntheOld))
  ALLOCATE(xOld(ntheOld))
  ALLOCATE(yOld(ntheOld))
  ALLOCATE(zOld(ntheOld))
  ALLOCATE(distance2derivsX(ntheOld))
  ALLOCATE(distance2derivsY(ntheOld))
  ALLOCATE(distance2derivsZ(ntheOld))

  ALLOCATE(theval(nthe))
  ALLOCATE(thetaVal(nthe))
  ALLOCATE(chiVal(nthe))

  !   move theta coordinates along one line

  ierr = 0

  DO k = 2, nzeta
     DO j = 1, npsi
  !      print*, 'mt1d_down: j, k = ', j, k
        distance(1) = 0._dp
        xOld(:) = xTsyg(:,j,k)
        yOld(:) = yTsyg(:,j,k)
        zOld(:) = zTsyg(:,j,k)
        DO  i = 2, ntheOld
           distance(i) = distance(i-1) + SQRT((xTsyg(i,j,k)-xTsyg(i-1,j,k))**2 +  &
                (yTsyg(i,j,k)-yTsyg(i-1,j,k))**2 +(zTsyg(i,j,k)-zTsyg(i-1,j,k))**2)
        END DO

        ! New coordinates
        DO  ii = 1, nthe
           theval(ii) = distance(ntheOld) * REAL(ii - 1,dp) / REAL(nthe - 1,dp)
        END DO
        thetaVal = theval * pi_d / distance(ntheOld)
        chiVal = (thetaVal + const_theta * SIN(2.*thetaVal)) * distance(ntheOld)/pi_d ! dimensions of distance

        ! print*, 'mt1d_down, here0: j, k, ntheOld = ', j, k, ntheOld, size(distance), size(xOld), size(distance2derivsX), &
        !     kind(distance), kind(xOld), kind(distance2derivsX)
!        print*, 'mt1d_down: distance = ', distance

        ! "Natural" splines
        CALL spline(distance(1:ntheOld), xOld(1:ntheOld), 1.E31_DP, 1.E31_DP, distance2derivsX(1:ntheOld))
        !        print*, '1'
        CALL spline(distance(1:ntheOld), yOld(1:ntheOld), 1.E31_DP, 1.E31_DP, distance2derivsY(1:ntheOld))
        !        print*, '2'
        CALL spline(distance(1:ntheOld), zOld(1:ntheOld), 1.E31_DP, 1.E31_DP, distance2derivsZ(1:ntheOld))
        !        print*, '3'

        ! print*, 'mt1d_down, here1: j, k, ntheOld = ', j, k

        DO i = 2, nthe-1
           xTsygNew(i,j,k) = splint(distance(1:ntheOld), xOld(1:ntheOld), distance2derivsX(1:ntheOld), chiVal(i))
           yTsygNew(i,j,k) = splint(distance(1:ntheOld), yOld(1:ntheOld), distance2derivsY(1:ntheOld), chiVal(i))
           zTsygNew(i,j,k) = splint(distance(1:ntheOld), zOld(1:ntheOld), distance2derivsZ(1:ntheOld), chiVal(i))
        END DO

        ! print*, 'mt1d_down, here: j, k = ', j, k

        xTsygNew(1,j,k) = xTsyg(1,j,k)
        xTsygNew(nthe,j,k) = xTsyg(ntheOld,j,k)
        yTsygNew(1,j,k) = yTsyg(1,j,k)
        yTsygNew(nthe,j,k) = yTsyg(ntheOld,j,k)
        zTsygNew(1,j,k) = zTsyg(1,j,k)
        zTsygNew(nthe,j,k) = zTsyg(ntheOld,j,k)

        !     PRINT*, 'j, k, R1, RN1 = ', j, k, &
        !          SQRT(xTsyg(1,j,k)**2 + yTsyg(1,j,k)**2 + zTsyg(1,j,k)**2), &
        !          SQRT(xTsygNew(1,j,k)**2 + yTsygNew(1,j,k)**2 + zTsygNew(1,j,k)**2)

     END DO
  END DO

  DEALLOCATE(distance)
  DEALLOCATE(xOld)
  DEALLOCATE(yOld)
  DEALLOCATE(zOld)
  DEALLOCATE(distance2derivsX)
  DEALLOCATE(distance2derivsY)
  DEALLOCATE(distance2derivsZ)
  DEALLOCATE(theval)
  DEALLOCATE(thetaVal)
  DEALLOCATE(chiVal)

  RETURN

END SUBROUTINE mapthe1d_down



