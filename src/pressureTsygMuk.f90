FUNCTION pressureTsygMuk(xEqGsm, yEqGsm)

  ! Tsyganenko-Mukai statistical pressure model
  ! from Tsyganenko, N. and T. Mukai, Tail plasma sheet models derived from Geotail particle data
  ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 108, NO. A3, 1136, doi:10.1029/2002JA009707, 2003

  ! Table 1 
  !Average hT i = 3.795 keV hNi = 0.625 cm 3 hPi = 0.229 nPa
  !Model rms 
  ! deviation
  ! 1.422 0.342 0.0427
  ! Correlation 0.708 0.567 0.955
  ! Equation number (4)        (5)         (6)

  ! PSW=1.94E-6*xn_pd*vel**2 

  USE nrtype
  USE Module1, ONLY : pdynGlobal, byimfGlobal, bzimfGlobal, pnormal

  IMPLICIT NONE

  REAL(DP) :: xEqGsm, yEqGsm, pressureTsygMuk
  REAL(DP), PARAMETER :: A1 = 0.057, A2 = 0.524, A3 = 0.0908, A4 = 0.527, &
       A5 = 0.078, A6 = -4.422, A7 = -1.533, A8 = -1.217, A9 = 2.54, &
       A10 = 0.32, A11 = 0.754, A12 = 1.048, A13 = -0.074, A14 = 1.015 
  REAL(DP) :: Bperp, rhoStar, rhoStar10, PSWStar, F, FStar, phi, theta, presAt10RE
  REAL(DP), EXTERNAL :: pressureRad

  rhoStar = 0.1_dp * SQRT(xEqGsm**2+yEqGsm**2) ! 1/10RE * rho

  phi = ATAN2(xEqGsm, yEqGsm) ! azimuthal angle, atan(y/x)

  PSWStar = pdynGlobal/3._dp

  theta = ATAN2(byimfGlobal,bzimfGlobal) ! IMF clock angle
  IF (theta < 0.0) theta = theta + twopi_d
  Bperp = SQRT(byimfGlobal**2 + bzimfGlobal**2) ! BSW perpendicular to Sun-Earth axis
  F = ABS(Bperp)*SQRT(SIN(0.5_dp*theta))
  FStar = 0.2_dp*F

!C  print*, 'pdyn, theta, F = ', pdynGlobal, theta, F

  pressureTsygMuk = A1*rhoStar**A6 + A2*PSWStar**A11*rhoStar**A7 + A3*FStar**A12*rhoStar**A8 + &
       (A4*PSWStar**A13*EXP(-A9*rhoStar) + A5*FStar**A14*EXP(-A10*rhoStar))*(SIN(phi))**2
  pressureTsygMuk = pressureTsygMuk/pnormal !In computational units

  ! Continuation to R < 10 RE (rhoStar < 1.)
  rhoStar10 = 1._dp
  presAt10RE = A1*rhoStar10**A6 + A2*PSWStar**A11*rhoStar10**A7 + A3*FStar**A12*rhoStar10**A8 + &
       (A4*PSWStar**A13*EXP(-A9*rhoStar10) + A5*FStar**A14*EXP(-A10*rhoStar10))*(SIN(phi))**2
  presAt10RE = presAt10RE/pnormal
  IF (rhoStar < 1._dp) pressureTsygMuk = presAt10RE/pressureRad(10._dp)  * &
       pressureRad(10._dp*rhoStar)*(1._dp-exp(1.*rhoStar-1.)) + pressureTsygMuk*exp(1.*rhoStar-1)



  RETURN 

END FUNCTION pressureTsygMuk


