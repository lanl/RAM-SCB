SUBROUTINE dps_general
  ! Calculates the depression on the Earth's surface - generalized Dessler-Parker-Sckopke relationship (Siscoe, 1970)
  ! Uses anisotropic pressure

  USE nrtype, ONLY : DP, pi_d
  USE Module1

  IMPLICIT NONE

  INTEGER               :: i, j, k, iplx
  REAL(DP)              ::  magneticEnergy(nthe), magneticEnergyInsideGeo(nthe), &
       magneticEnergyDipole, thermalEnergy(nthe), thermalEnergyInsideGeo(nthe), &
       rsq, totalEnergy, volumeTotal 

  magneticEnergy = 0.0_dp
  magneticEnergyInsideGeo = 0.0_dp
  magneticEnergyDipole = 0.0_dp
  thermalEnergy = 0.0_dp
  thermalEnergyInsideGeo = 0.0_dp
  totalEnergy = 0.0_dp
  volumeTotal = 0.0_dp

  DO i = 2, nthe-1
     DO j = 2, npsi
        DO k = 2, nzeta
           rsq = x(i,j,k)**2 + y(i,j,k)**2
           magneticEnergy(i) = magneticEnergy(i) + jacobian(i,j,k) * dr * dpPrime * dt * &
                & (bsq(i,j,k)*bnormal**2) / (2._dp*mu0) * 1.E-18_dp*REarth**3 ! In Joules now
           IF (rsq < 6.6_dp**2) magneticEnergyInsideGeo(i) = magneticEnergyInsideGeo(i) + &
                jacobian(i,j,k) * dr * dpPrime * dt * &
                & (bsq(i,j,k)*bnormal**2) / (2._dp*mu0) * 1.E-18_dp*REarth**3  
           thermalEnergy(i) = thermalEnergy(i) + jacobian(i,j,k) * dr * dpPrime * dt * &
                & (pper(i,j,k) + 0.5_dp*ppar(i,j,k))*pnormal * 1.E-9_dp * REarth**3 ! In Joules now
           IF (rsq < 6.6_dp**2) thermalEnergyInsideGeo(i) = thermalEnergyInsideGeo(i) + &
                jacobian(i,j,k) * dr * dpPrime * dt * &
                & (pper(i,j,k) + 0.5_dp*ppar(i,j,k))*pnormal * 1.E-9_dp * REarth**3
           volumeTotal = volumeTotal + jacobian(i,j,k) * dr * dpPrime * dt 
        END DO
     END DO
     !C WRITE(*, '(A, 1X, I3, 1X, E12.3)') 'i, Magnetic energy = ', i, magneticEnergy(i)
     !C WRITE(*, '(A, 1X, I3, 1X, E12.3)') 'i, Thermal energy = ', i, thermalEnergy(i)
  END DO

  magneticEnergyDipole = 4._dp*pi_d/(3._dp*mu0) * BEarth**2 * REarth**3 ! In J

  totalEnergy = SUM(magneticEnergy) + SUM(thermalEnergy)

!C  PRINT*, ' '; CALL FLUSH(6)
  !  WRITE(*, '(A, 1X, E12.3)') 'Total magnetic energy = ', SUM(magneticEnergy); CALL FLUSH(6)
  !  WRITE(*, '(A, 1X, E12.3)') 'Magnetic dipole energy = ', magneticEnergyDipole; CALL FLUSH(6)
  !  WRITE(*, '(A, 1X, E12.3)') 'Total thermal energy = ', SUM(thermalEnergy); CALL FLUSH(6)
  !  WRITE(*, '(A, 1X, E12.3)') 'Total energy = ', totalEnergy; CALL FLUSH(6)
  !  WRITE(*, '(A, 1X, F10.3)') 'Total volume of the computational cavity = ', volumeTotal, ' RE^3'; CALL FLUSH(6)
  DstDPS = 1.3_dp * (-BEarth) * (2._dp*SUM(thermalEnergy))/(3._dp*magneticEnergyDipole) * 1.E9_dp
  DstDPSInsideGeo = 1.3_dp * (-BEarth) * (2._dp*SUM(thermalEnergyInsideGeo))/(3._dp*magneticEnergyDipole) * 1.E9_dp
  WRITE(*, '(A, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, A)') 'DstDPS, DstDPSGeo, DstBiot, DstBiotGeo = ', real(DstDPS), &
       real(DstDPSInsideGeo), real(DstBiot), real(DstBiotInsideGeo), ' nT' ; print*, ' '; CALL FLUSH(6)! 1.3 factor due to currents induced in the Earth 
  !C WRITE(iUnitDst, *) 'Hour, DstDPS, DstDPSGeo, DstBiot, DstBiotGeo: ', real(rHour), real(Dst), real(DstInsideGeo), &
  !C     real(DstBiot), real(DstBiotInsideGeo)  ! Written in run_ramscb
  CALL flush(559)
  

  RETURN

END SUBROUTINE dps_general



