FUNCTION pofpsiCS(psi0, psi1, psi2, p1, p2, delta, psi)

  USE nrtype, ONLY : dp, pi_d

  IMPLICIT NONE

  INTERFACE erf
     FUNCTION erf_s(x)
       USE nrtype
       USE nr, ONLY : gammp
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: erf_s
     END FUNCTION erf_s
     FUNCTION erf_v(x)
       USE nrtype
       USE nr, ONLY : gammp
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: erf_v
     END FUNCTION erf_v
  END INTERFACE

  REAL(DP), INTENT(IN)  :: psi0, psi1, psi2, p1, p2, delta
  REAL(DP), INTENT(IN)  :: psi(:)
  REAL(DP), DIMENSION(SIZE(psi))  :: pofpsiCS

  REAL(DP), DIMENSION(SIZE(psi)) :: x
  REAL(DP)              :: x1, x2, er1, er2, const, p0
  REAL(DP)              :: sqrtPi
  INTEGER               :: ierr, ierrdeal

  sqrtPi = SQRT(pi_d)
  x1 = (psi1-psi0)/delta
  x2 = (psi2-psi0)/delta
  x = (psi - psi0)/delta
  er1 = erf(x1)
  er2 = erf(x2)

  const = (p2-p1) / (delta * (0.5_dp*sqrtPi * er2 + x2 - 0.5_dp*sqrtPi * er1 - x1))
  p0 = ((p1*x2 - p2*x1) + 0.5_dp*sqrtPi * (p1*er2 - p2*er1)) / (0.5_dp*sqrtPi*(er2-er1) + x2-x1)
!  pofpsiCS = const * delta * (0.5_dp*sqrtPi*erf(x) + x) + p0

pofpsiCS = - delta * 1./9.255 * (0.5*sqrtPi * (erf(x) - er2) + 0.1*(x-x2)) + p2

END FUNCTION pofpsiCS
