REAL FUNCTION pressureRad(radius)

  USE nrtype
  USE Module1

  IMPLICIT NONE

  REAL(DP) :: radius, LargeA, LA, dPdRRight, dPdR2Right, Acap, Bcap, Ccap, Dcap, &
       Acap2, Bcap2, Ccap2
  REAL(DP) :: m, n, pUp, pDown, x1, x2, delta1, delta2, pressureSK, pUp2, pDown2

  ! This subroutine outputs pressure as a function of radial distance r in the
  ! equatorial plane, midnight 

  ! Spence - Kivelson empirical formula
  IF (iCorrectedPressure /= 1) THEN
     !     IF (ABS(radius) > 3.5)
     pressureRad = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))

     LargeA = pressurequot / pnormal * 89.*EXP(-0.59*3.5) + 8.9 * 3.5**(-1.53) ! value at x = 3
     !     IF (ABS(radius) <= 3.5) pressureRad = pressurequot/ pnormal*LargeA**(((ABS(radius) - 2.5)))
     ! Took this out for now

     RETURN 
  END IF

  ! GOTO 500
  ! "corrected" Spence-Kivelson pressure for test runs, with a smaller grad p at large distances

  IF (iCorrectedPressure == 1) THEN
     ! This is the profile used for obtaining thin CS (FLR calculations)
     pressureSK = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))
     m = 4.
     n = 2.
     ! pUp = m * pressureSK
     pUp = 10. / pnormal   
     pDown = n * pressureSK
     x1 = 7.5
     x2 = 9.25
     delta1 = 3.0_dp
     delta2 = 3.0_dp
     pUp2 = pUp / 2. * (1. + TANH((x1 - ABS(radius)) / delta1))
     pDown2 = pDown / 2. * (1. + TANH((ABS(radius) - x2)/delta2))
     pressureRad = pUp2 + pDown2
  END IF


500 CONTINUE

  go to 1000

  IF (iCorrectedPressure == 1) THEN
     ! This is the profile used for obtaining thin CS (Annales paper)
     pressureSK = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))
     m = 4.
     n = 2.
     !     n = 1.
     ! pUp = m * pressureSK
     pUp = 10. / pnormal   
     pDown = n * pressureSK
     x1 = 8.   
     x2 = 9.25
     delta1 = 2.5_dp
     delta2 = 2.5_dp

     pUp2 = pUp / 2. * (1. + TANH((x1 - ABS(radius)) / delta1))
     pDown2 = pDown / 2. * (1. + TANH((ABS(radius) - x2)/delta2))
     pressureRad = pUp2 + pDown2

     ! Try this, for P(\psi, \alpha) functional case
     ! pressureRad = 2._dp * pressureSK
  END IF


  IF (iCorrectedPressure == 1) THEN
     ! This is the profile used for obtaining thin CS - do not TOUCH !!!
     pressureSK = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))
     m = 3.
     n = 1.
     !     pUp = m * pressureSK
     pUp = 10. * (20. - ABS(radius))**3/16.**3 / pnormal
     pDown = n * pressureSK

     x1 = 10.   
     x2 = 12.
     !     delta1 = 1.5_dp
     delta1 = 2._dp
     delta2 = 1.5_dp

     pUp2 = pUp / 2. * (1. + TANH((x1 - ABS(radius)) / delta1)) 
     pDown2 = pDown / 2. * (1. + TANH((ABS(radius) - x2)/delta2))
     pressureRad = pUp2 + pDown2

     ! New pressure profile, with only one tanh
     !  delta = 0.5
     !  pressureRad = pDown + 0.5 * (pUp -pDown) * (1. + tanh((x1 - abs(radius)) / delta)); 

  END IF

1000 CONTINUE

  RETURN 

END FUNCTION pressureRad


