SUBROUTINE trapzd(func,a,b,s,n)
  USE nrtype; USE nrutil, ONLY : arth
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP), INTENT(INOUT) :: s
  INTEGER(I4B), INTENT(IN) :: n
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  REAL(DP) :: del,fsum
  INTEGER(I4B) :: it
  IF (n == 1) THEN
     s = 0.5_dp*(b-a)*SUM(func( (/ a, b /) ))
  ELSE
     it = 2**(n-2)
     del=(b-a)/REAL(it, dp)
     fsum = SUM(func(arth(a+0.5_dp*del,del,it)))
     s = 0.5_dp*(s+del*fsum)
  END IF

  RETURN
END SUBROUTINE trapzd
