SUBROUTINE findR

  ! Obtains equidistant x-s on the equatorial plane, at a local time where the plasma gradients are strong (e.g. between midnight and dusk
  ! during a storm main phase)

  USE nrtype
  USE Module1, ONLY: iAzimOffset, kMax, x, y, &
	rank, nThetaEquator, npsi, nzeta, nZetaMidnight, radEqMidNew

  IMPLICIT NONE

  REAL(DP) :: deltaR, distConsecFluxSqOld, distConsecFluxSq
  REAL(DP) :: radius(npsi)
  INTEGER :: j, k

  distConsecFluxSq = 0._dp
  distConsecFluxSqOld = 0._dp

  IF (iAzimOffset == 2) THEN
     DO k = 2, nzeta
        do j = 2, npsi
           distConsecFluxSq = (x(nThetaEquator,j,k) - x(nThetaEquator,j-1,k))**2 + &
                (y(nThetaEquator,j,k) - y(nThetaEquator,j-1,k))**2 
           ! distConsecFluxSq = x(nThetaEquator,j,k)**2 - x(nThetaEquator,j-1,k)**2 + &
           !     y(nThetaEquator,j,k)**2 - y(nThetaEquator,j-1,k)**2 
           IF (distConsecFluxSq > distConsecFluxSqOld) THEN
              distConsecFluxSqOld = distConsecFluxSq
              kMax = (3*nzeta)/8 !C k
           END IF
        END DO
     end DO
     !C PRINT*, 'rank: findR: kMax = ', rank, kMax; CALL FLUSH(6)
     DO j = 1, npsi
        radius(j) = SQRT(x(nThetaEquator,j,kMax)**2 + &
             y(nThetaEquator,j,kMax)**2)
     END DO
  ELSE IF (iAzimOffset == 1) THEN
     DO j = 1, npsi
        radius(j) = SQRT(x(nThetaEquator,j,nZetaMidnight)**2 + &
             y(nThetaEquator,j,nZetaMidnight)**2)
     END DO
  END IF

  deltaR = (radius(npsi) - radius(1)) / REAL(npsi-1,dp)

  DO j = 2, npsi
     radEqMidNew(j) = radius(1) + REAL(j-1,dp)*deltaR
  END DO

  radEqMidNew(1) = radius(1)

  RETURN

END SUBROUTINE findR
