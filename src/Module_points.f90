MODULE Module_points
  USE nrtype, ONLY : SP, DP
  use Module1, ONLY: npsi, nzeta
  IMPLICIT NONE

  INTEGER :: last
  INTEGER :: ilat, ilon
  REAL(DP) :: xx01(100000), yy01(100000), zz01(100000)
  REAL(DP), ALLOCATABLE :: xTsyg(:,:,:), yTsyg(:,:,:), zTsyg(:,:,:), &
	xTsygNew(:,:,:), yTsygNew(:,:,:), zTsygNew(:,:,:)

  INTEGER        :: nthe, ntheOld, iQuiet, iDay, iTime
  REAL(DP)       :: rlatval(500), rlonval(500), psival(500), rlon(500), &
       rlat(500), funcModBeta(500), funcModBetaLarge(5000)
  ! Modifies the alpha Euler potential
  REAL(DP)       :: radiusStart

  REAL(DP) :: xpsiin ,xpsiout

END MODULE Module_points
