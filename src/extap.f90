SUBROUTINE extap(x1,x2,x3,x4)

  USE nrtype
  IMPLICIT NONE
  REAL(DP), INTENT(IN)     :: x1, x2, x3
  REAL(DP), INTENT(OUT)    :: x4
  REAL(DP)                 :: ddx1, ddx2, ddx, pm

  x4 = 3. * x3 - 3. * x2 + x1
  ddx1 = x3 - x2
  ddx2 = x2 - x1
  ddx = x4 - x3
  pm = ddx * ddx1
  IF(pm > 0.) RETURN
  IF(ddx2 == 0) x4 = 2. * x3 - x2
  IF(ddx2 == 0) RETURN
  ddx = (ddx1 * ddx1) / ddx2
  x4 = x3 + ddx
  RETURN
END SUBROUTINE extap
