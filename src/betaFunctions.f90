SUBROUTINE betaFunctions

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTERFACE Spline_1D
     SUBROUTINE Spline_derivs_1D_zeta(x1, f, func, derivsX1)
       USE nrtype, ONLY : DP
       USE EZspline_obj ! import the modules
       USE EZspline  

       IMPLICIT NONE

       INTEGER, PARAMETER :: r8 = DP 
       REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8

       REAL(r8), DIMENSION(:), INTENT(IN) :: x1 ! independent variable
       REAL(r8), DIMENSION(:), INTENT(IN) :: f
       REAL(r8), DIMENSION(:), INTENT(OUT) :: func, derivsX1 ! interpolated values
       INTEGER n1, ier, BCS1(2)
       TYPE(EZspline1_r8) :: spline_o ! 1-D EZspline object
       INTEGER k1
       REAL(r8), DIMENSION(:), ALLOCATABLE :: z1
     END SUBROUTINE Spline_derivs_1D_zeta
  END INTERFACE

  REAL(DP) :: betaValue(nzetp)
  INTEGER :: k

  CALL Spline_1D(zetaVal(1:nzeta), alphaval(1:nzeta), betaValue(1:nzeta), fzet(1:nzeta))
  CALL Spline_1D(zetaVal(1:nzeta), f(1:nzeta), betaValue(1:nzeta), fzetp(1:nzeta))
 
  fzet(1) = fzet(nzeta)
  fzetp(1) = fzetp(nzeta)
  fzet(nzetp) = fzet(2)
  fzetp(nzetp) = fzetp(2)

  RETURN

END SUBROUTINE betaFunctions



