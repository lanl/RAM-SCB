module ModScbParams

  use ModScbMain, ONLY: DP

  implicit none
  save

  integer :: psiChange = 1
  integer :: theChange = 1

  integer :: iDumpRAMFlux = 0 ! Writes RAM flux mapped along 3D field lines in NetCDF format
  integer :: MinSCBIterations = 1

  real(DP) :: decreaseConvAlphaMin = 5e-1,  &
              decreaseConvAlphaMax = 5e-1,  &
              decreaseConvPsiMin   = 5e-1,  &
              decreaseConvPsiMax   = 5e-1,  &
              blendAlphaInit       = 0.20,  &
              blendPsiInit         = 0.20,  &
              blendMin             = 0.01,  &
              blendMax             = 0.5,   &
              InConAlpha           = 1e-6,  &
              InConPsi             = 1e-6

 !integer  :: method = 1 ! Direct matrix inversion calculation of Euler Potentials
  integer  :: method = 2 ! Iterative SOR calculation of Euler Potentials
 !integer  :: method = 3 ! No SCB calculation

  integer  :: isotropy = 0 ! Anisotropic pressure case 
 !integer  :: isotropy = 1 ! Isotropic pressure case

  integer  :: iAMR = 0 ! Mesh refinement in magnetic flux, so that one has equidistant magnetic flux surfaces; improves convergence a lot

  integer  :: iReduceAnisotropy = 0 ! No change in anisotropy 
 !integer  :: iReduceAnisotropy = 1 ! Change anisotropy to marginally mirror-stable

  integer  :: iOuterMethod = 1 ! Picard iteration
 !integer  :: iOuterMethod = 2 ! Newton iteration  

  integer  :: iWantAlphaExtrapolation = 1 ! Extrapolate alpha (beta) on the first/last flux surface
 !integer  :: iWantAlphaExtrapolation = 0 ! Do not extrapolat alpha (beta) on  !the first/last flux surface

 !integer  :: iAzimOffset = 2 ! Equidistance sought for most problematic local time
  integer  :: iAzimOffset = 1 ! Equidistance maintained at most stretched MLT

  integer  :: isSORDetailNeeded = 0 ! No details for the inner SOR iterations 
 !integer  :: isSORDetailNeeded = 1 ! Details about inner SOR iterations

  integer  :: isEnergDetailNeeded = 1 ! Dst computation (DPS formula with thermal energy inside domain)
 !integer  :: isEnergDetailNeeded = 0

  integer  :: isFBDetailNeeded = 0 ! Does not compute global force imbalance
 !integer  :: isFBDetailNeeded = 1 ! Computes global force imbalance after SCB calculation

  integer  :: isPressureDetailNeeded = 0 ! Does not output SCB pressure map
 !integer  :: isPressureDetailNeeded = 1 ! Outputs SCB pressure map

  integer  :: iLossCone = 1 ! Filled loss cone
 !integer  :: iLossCone = 2 ! More realistic, empty loss cone for RAM  !computations (M. Liemohn's formalism, Liemohn, 2004)

end Module ModScbParams
