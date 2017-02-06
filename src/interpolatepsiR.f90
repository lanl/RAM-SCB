SUBROUTINE InterpolatePsiR   

  ! Interpolates the new values of psi at the new locations xnew on midnight 
  ! equator

  USE nrtype
  USE Module1

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

  REAL(DP), DIMENSION(npsi) :: radEqmid, psival1D, psi2Deriv
  INTEGER :: ialloc, ierr, j

  IF (iAzimOffset == 1) THEN
     radEqmid(1:npsi) = SQRT(x(nThetaEquator, 1:npsi, nZetaMidnight)**2 + &
          y(nThetaEquator, 1:npsi, nZetaMidnight)**2)
  ELSE IF (iAzimOffset == 2) THEN  
     radEqmid(1:npsi) = SQRT(x(nThetaEquator, 1:npsi, kMax)**2 + &
          y(nThetaEquator, 1:npsi, kMax)**2)
  END IF

  psival1D(1:npsi) = psival(1:npsi)

  ! print*, 'I am here in intpsi.'; call flush(6)

  ! do j = 1, npsi
  !   print*, 'j, x(j), radEqmid = ', j, x(nThetaEquator,j,nZetaMidnight), radEqmid(j)
  ! end do
  
  CALL spline(radEqmid, psival1D, 1.E31_dp, 1.E31_dp, psi2Deriv)

  DO j = 2, npsi-1
     psival(j) = splint(radEqmid, psival1D, psi2Deriv, radEqmidNew(j))
  END DO


  RETURN

END SUBROUTINE InterpolatePsiR


