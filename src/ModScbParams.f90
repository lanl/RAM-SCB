!==============================================================================
module ModScbParams

  use ModScbMain, ONLY: DP

  implicit none
  save

  integer :: iDumpRAMFlux = 0 ! Writes RAM flux mapped along 3D field lines in NetCDF format

  real(DP) :: decreaseConvAlphaMin = 5e-1,  &
              decreaseConvAlphaMax = 5e-1,  &
              decreaseConvPsiMin   = 5e-1,  &
              decreaseConvPsiMax   = 5e-1,  &
              blendAlphaInit       = 0.20,  &
              blendPsiInit         = 0.20,  &
              blendMin             = 0.01,  &
              blendMax             = 0.5

end Module ModScbParams
!==============================================================================
