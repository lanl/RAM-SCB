FUNCTION anisEq(radius)

  USE nrtype
  USE Module1

  IMPLICIT NONE

  REAL(DP), INTENT(IN)  :: radius
  real(DP) :: anisEq
  REAL(DP) :: aLui(6)

  ! This subroutine outputs pressure anisotropy on the equator anisEq = P_perp/P_par 
  ! as a function of radial distance r using the Lui (1994) empirical formula

  aLui(1) = -0.410554
  aLui(2) = -9.94369
  aLui(3) = -86.9877
  aLui(4) = -504.066
  aLui(5) = -1110.73
  aLui(6) = -847.912

  anisEq = aLui(1) + aLui(2)*1./radius + aLui(3)*1./((radius))**2 + &
       aLui(4)*1./((radius))**3 + aLui(5)*1./((radius))**4 + &
       aLui(6)*1./((radius))**5


  RETURN
END FUNCTION anisEq


