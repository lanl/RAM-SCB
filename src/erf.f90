FUNCTION erf_s(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: erf_s
  erf_s=gammp(0.5_DP,x**2)
  IF (x < 0.0) erf_s=-erf_s
END FUNCTION erf_s


FUNCTION erf_v(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(SIZE(x)) :: erf_v
  erf_v=gammp(SPREAD(0.5_DP,1,SIZE(x)),x**2)
  WHERE (x < 0.0) erf_v=-erf_v
END FUNCTION erf_v
