FUNCTION locate(xx,x)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  INTEGER(I4B) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=SIZE(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  DO
     IF (ju-jl <= 1) EXIT
     jm=(ju+jl)/2
     IF (ascnd .EQV. (x >= xx(jm))) THEN
        jl=jm
     ELSE
        ju=jm
     END IF
  END DO
  IF (x == xx(1)) THEN
     locate=1
  ELSE IF (x == xx(n)) THEN
     locate=n-1
  ELSE
     locate=jl
  END IF
END FUNCTION locate
