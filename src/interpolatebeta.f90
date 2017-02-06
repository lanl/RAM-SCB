SUBROUTINE InterpolateBeta

  ! Interpolates the new values of psi at the new locations xnew on midnight 
  ! equator

  USE nrtype
  USE Module1

  IMPLICIT NONE

  interface Spline_1D
     SUBROUTINE Spline_1D_periodic(x_1, f, x, func, ierrDomain)
       USE nrtype, ONLY : DP
       USE EZspline_obj ! import the modules
       USE EZspline  
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP 
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1 ! independent variable
       REAL(r8), DIMENSION(:), INTENT(IN) :: f
       REAL(r8), DIMENSION(:), INTENT(OUT) :: func   ! interpolated values 
       INTEGER, INTENT(OUT) :: ierrDomain
       REAL(r8), DIMENSION(:), INTENT(IN)    :: x  ! grid of points for output
     end SUBROUTINE Spline_1D_periodic
  end interface Spline_1D

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
  END INTERFACE splint_interface

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
  END INTERFACE spline_interface

  REAL(DP), DIMENSION(nzetp) :: radEqmid, alphaval1D, alpha2Deriv
  INTEGER :: ialloc, ierr, k

  
  radEqmid(1) = 0._dp

  DO k = 2, nzeta
     radEqmid(k) = radEqmid(k-1) + SQRT((x(nThetaEquator, jMax, k)-x(nThetaEquator,jMax,k-1))**2 + &
          (y(nThetaEquator, jMax, k)-y(nThetaEquator,jMax,k-1))**2)
  END DO

  alphaval1D(1:nzeta) = alphaval(1:nzeta)

  ! do j = 1, npsi
  !   print*, 'j, x(j), radEqmid = ', j, x(nThetaEquator,j,nZetaMidnight), radEqmid(j)
  ! end do

  !C CALL spline(radEqmid(1:nzeta), alphaval1D(1:nzeta), 1.E31_dp, 1.E31_dp, alpha2Deriv(1:nzeta))
  !C DO k = 2, nzeta-1
  !C   alphaval(k) = splint(radEqmid(1:nzeta), alphaval1D(1:nzeta), alpha2Deriv(1:nzeta), radEqmidNew(k))
  !C END DO


  CALL Spline_1D(radEqmid(1:nzeta), alphaval1D(1:nzeta), radEqmidNew(1:nzeta), alphaval(1:nzeta), ierr)

  ! Periodicity
  alphaval(1) = alphaval(nzeta) - 2.*pi_d
  alphaval(nzetp) = alphaval(2) + 2.*pi_d


  RETURN

END SUBROUTINE InterpolateBeta


