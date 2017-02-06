FUNCTION qtrap(func,a,b)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : trapzd
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP) :: qtrap
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: JMAX=20
  REAL(DP), PARAMETER :: EPS=1.0e-6_DP, UNLIKELY=-1.0e30_DP
  REAL(DP) :: olds
  INTEGER(I4B) :: j
  olds=UNLIKELY
  DO j=1, JMAX
     CALL trapzd(func,a,b,qtrap,j)
     ! print*, 'qtrap: j, qtrap = ', j, qtrap
     IF (ABS(qtrap-olds) < EPS*ABS(olds)) RETURN
     IF (qtrap == 0.0 .AND. olds == 0.0 .AND. j > 6) RETURN
     olds=qtrap
  END DO
  CALL nrerror('qtrap: too many steps')
END FUNCTION qtrap
