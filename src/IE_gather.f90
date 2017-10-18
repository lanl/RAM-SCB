
subroutine IE_gather

!  use ModMpi
  use ModIonosphere
  use ModKind
  implicit none

  real(real8_) :: localVar(2*IONO_nTheta-1, IONO_nPsi)
  integer :: iError, iSize
  !--------------------------------------------------------------------------
  IONO_phi = -1.0e32
  IONO_IonNumFlux = -1.0e32
  IONO_Joule = -1.0e32
  IONO_Jr = -1.0e32
  IONO_Ave_E = -1.0e32
  IONO_Eflux = -1.0e32

  IONO_phi(1:IONO_nTheta,:) = IONO_NORTH_Phi
  IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_NORTH_IonNumFlux
  IONO_Joule(1:IONO_nTheta,:) = IONO_NORTH_Joule
  IONO_Jr(1:IONO_nTheta,:) = IONO_NORTH_Jr
  IONO_Ave_E(1:IONO_nTheta,:) = IONO_NORTH_Ave_E
  IONO_Eflux(1:IONO_nTheta,:) = IONO_NORTH_EFlux
  
  IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Phi
  IONO_IonNumFlux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_IonNumFlux
  IONO_Joule(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Joule
  IONO_Jr(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Jr
  IONO_Ave_E(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Ave_E
  IONO_Eflux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_EFlux
  
end subroutine IE_gather
