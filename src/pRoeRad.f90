FUNCTION pRoeRad(radius)

  USE nrtype, ONLY : DP

  IMPLICIT NONE

  REAL(DP) :: radius, pRoeRad

  ! This subroutine outputs pressure as a function of radial distance r in the
  ! equatorial plane

  ! Fit of Roeder pressure
  !  pRoeRad = 1220.*exp(-3.8*radius) + 89.*exp(-0.6*radius)
  pRoeRad = 8.4027*exp(-1.7845*radius) + 90.3150*exp(-0.7659*radius)

  RETURN 

END FUNCTION pRoeRad


