FUNCTION qromb(func,a,b)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : polint, trapzd
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a, b
  REAL(DP) :: qromb
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: JMAX=50, JMAXP=JMAX+1, K=2, KM=K-1 ! With K = 2 this is equivalent to Simpson rule
  REAL(DP), PARAMETER :: EPS = 1.0e-3_dp
  REAL(DP), DIMENSION(JMAXP) :: h, s
  REAL(DP) :: dqromb
  INTEGER(I4B) :: j
  h(1) = 1.0_dp
  DO j = 1, JMAX
     CALL trapzd(func,a,b,s(j),j)
     ! print*, 'qromb: j, s(j) = ', j, s(j)
     IF (j >= K) THEN
        CALL polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
        IF (ABS(dqromb) <= EPS*ABS(qromb)) RETURN
     END IF
     s(j+1)=s(j)
     h(j+1)=0.25_dp*h(j)
  END DO
  CALL nrerror('qromb: too many steps')
END FUNCTION qromb
