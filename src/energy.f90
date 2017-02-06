SUBROUTINE energy

  USE nrtype, ONLY : DP 
  USE Module1

  IMPLICIT NONE

  INTEGER               :: i, j, k, iplx
  REAL(DP)              ::  magneticEnergy, thermalEnergy, totalEnergy, volumeTotal 

  magneticEnergy = 0.0_dp
  thermalEnergy = 0.0_dp
  totalEnergy = 0.0_dp
  volumeTotal = 0.0_dp

  IF (isotropy == 1) THEN 
     DO i = 2, nthe-1
        DO j = 2, npsi-1
           DO k = 2, nzeta
              IF ((ABS(x(i,j,k)) > 4._dp)) THEN  
                 ! Only for |X| > 5 R_E, closer distances not considered 
                 magneticEnergy = magneticEnergy + jacobian(i,j,k) * dr * dpPrime * dt * &
                      & 0.5_dp * bsq(i,j,k)
                 thermalEnergy = thermalEnergy + jacobian(i,j,k) * dr * dpPrime * dt * &
                      & 1.5_dp * pressure3D(i,j,k)
                 volumeTotal = volumeTotal + jacobian(i,j,k) * dr * dpPrime * dt
              END IF
           END DO
        END DO
     END DO
  END IF

  totalEnergy = magneticEnergy + thermalEnergy

  WRITE(*, '(A, 1X, E12.3)') 'Magnetic energy = ', magneticEnergy
  WRITE(*, '(A, 1X, E12.3)') 'Thermal energy = ', thermalEnergy
  WRITE(*, '(A, 1X, E12.3)') 'Total energy = ', totalEnergy
  WRITE(*, '(A, 1X, E12.3)') 'Total volume = ', volumeTotal

  

  RETURN

END SUBROUTINE energy



