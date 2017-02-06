SUBROUTINE psiFunctions

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTERFACE Spline_1D
     SUBROUTINE Spline_derivs_1D(x1, f, func, derivsX1)
       USE EZspline_obj ! import the modules
       USE EZspline  
       USE nrtype, ONLY : DP
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
     END SUBROUTINE Spline_derivs_1D
  END INTERFACE

  REAL(DP) :: psiValue(npsi)
  INTEGER :: j

  CALL Spline_1D(rhoVal(1:npsi), psival(1:npsi), psiValue, f(1:npsi))
  CALL Spline_1D(rhoVal(1:npsi), f(1:npsi), psiValue, fp(1:npsi))

 !C OPEN(10, file = 'psi_values', action = 'write', status = 'replace')
 !C DO j = 1,npsi
 !C    WRITE(10, '(I3, 1X, 3G14.7)') j, psival(j), f(j), fp(j)
 !C END DO
 !C CLOSE(10)
 

  RETURN
END SUBROUTINE psiFunctions



