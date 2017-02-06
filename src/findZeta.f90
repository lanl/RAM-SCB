SUBROUTINE findZeta

  ! Obtains equidistant azimuthal points on the equatorial plane, at a radial location where the plasma gradients are strong

  USE nrtype
  USE Module1, ONLY: iAzimOffset, jMax, x, y, &
       rank, nThetaEquator, npsi, nzeta, nzetp, nZetaMidnight, radEqMidNew, pper, bsq

  IMPLICIT NONE

  REAL(DP) :: deltaR, beta_plasma
  REAL(DP) :: radius(nzetp)
  INTEGER :: j, k

  beta_plasma = 0._dp
  DO j = 1, npsi
      IF (MAXVAL(pper(:,j,2:nzeta)/bsq(:,j,2:nzeta)) > beta_plasma) THEN
        beta_plasma = MAXVAL(pper(:,j,2:nzeta)/bsq(:,j,2:nzeta))
        jMax = j
     END IF
  END DO

  ! PRINT*, 'rank: findZeta: jMax = ', rank, jMax; CALL FLUSH(6)

  radius(1) = 0._dp
  DO k = 2, nzeta
     radius(k) = radius(k-1) + SQRT((x(nThetaEquator,jMax,k)-x(nThetaEquator,jMax,k-1))**2 + &
          (y(nThetaEquator,jMax,k)-y(nThetaEquator,jMax,k-1))**2)
  END DO

  deltaR = (radius(nzeta) - radius(1)) / REAL(nzeta-1,dp)

  DO k = 2, nzeta
     radEqMidNew(k) = radius(1) + REAL(k-1,dp)*deltaR
  END DO

  radEqMidNew(1) = radius(1)

  RETURN

END SUBROUTINE findZeta
