MODULE ModScbMain

  USE nrtype,     ONLY: pi_d, SP, DP ! Definition of single/double precision type
  USE ModRamMain, ONLY: PathRamOut, PathScbOut, PathScbIn

  IMPLICIT NONE
  SAVE

  REAL(DP), PARAMETER :: mu0 = 4._dp*pi_d*1.E-7_dp, &
                         BEarth = 0.31_dp*1.E-4_dp, &
                         REarth = 6.4_dp*1.E6_dp 

!!! Intilization
  integer  :: method               = 2,     &
              iCountPressureCall   = 0,     &
              iElectric            = 0
              
  REAL(SP) :: rHour


  CHARACTER(LEN=2)   :: Day, DayP1, HourInteger, HourIntegerP1, MinChar
  CHARACTER(Len=5)   :: HourDecimal, HourDecShort
  CHARACTER(LEN=4)   :: suffixFilePress
  CHARACTER(LEN=500) :: fileNamePressure
  CHARACTER(LEN=200) :: prefixPres
  CHARACTER(LEN=2)   :: Hour

  integer :: tsygcorrect = 0

  integer :: isotropy = 0 ! Anisotropic pressure case 

  integer :: iInducedE = 0 ! No calculation/writing of induced E-fields    
  integer :: iConvE = 1 ! No calculation/writing of convective E-fields (if zero)

  ! For runs with RAM pressures, only the following choices should be used

  integer :: iReduceAnisotropy = 0 ! No change in anisotropy 
  !iReduceAnisotropy = 1 ! Change anisotropy to marginally mirror-stable

  integer :: iOuterMethod = 1 ! DO NOT CHANGE FOR RAM-SCB ! Picard iteration  ! Newton has not been benchmarked for anisotropic pressure
  !iOuterMethod = 2 ! Newton iteration  

  integer :: iInterpMethod = 1 ! DO NOT CHANGE FOR RAM-SCB ! PSPLINE interpolation

  !method = 2 ! DO NOT CHANGE FOR RAM-SCB 
  ! SOR computation of PDE; faster than direct method for N larger than about 41
  ! ! method = 3 means no calculation
  ! method NOW SET IN INIT_RAMSCB.F90; DEFAULTS TO 2.

  ! iWantAlphaExtrapolation = 0 ! Do not extrapolate alpha (beta) on the
  ! first/last flux surface
  integer :: iWantAlphaExtrapolation = 1 ! Extrapolate alpha (beta) on the first/last flux surface

  integer :: numit = 200 ! DO NOT CHANGE FOR RAM-SCB Max. # of iterations; the actual number is much lower (code exits on convergence) and depends on convergence speed 

  integer :: iConvergenceMethod = 2 ! 1 (global force balance decrease) or 2 (residual between iterations decrease)  

  ! decreaseConv - the amount the global force imbalance decreases throughout
  ! iterations (vs. the initial one)
  ! decreaseConv = 0.5_dp
  real(dp) :: decreaseConv = 0.5 ! 3X reduction, should be enough for RAM-SCB (very hard
  ! to get better w/ coarse grid resolution)
  ! decreaseConv = 0.2   ! 5X reduction, if grid resolution is better
  ! decreaseConv = 0.1   ! 10X reduction, and more needed if state has to be a
  ! near-perfect equilibrium for theoretical calculations

  real(dp) :: blendInitial = 0.25_dp ! "Blends" solution at iterations (n+1) and n

  real(dp) :: thresh = 1.1
  real(dp) :: damp = 0.5

  real(dp) :: relax = 1.5
  integer :: nrelax = 5 ! Try to increase blending factor after nrelax iterations

  integer :: iSm = 0 !Even smoother grid respacing, for difficult equilibria
  integer :: iSm2 = 1 ! Savitzky-Golay smoothing in pressure.f90 (for RAM-SCB mostly)

  integer :: nimax = 5000  ! DO NOT CHANGE FOR RAM-SCB ! # of inner iterations in SOR technique; usually 2000 is more than enough for SOR convergence

  integer :: iAMR = 1 ! Mesh refinement in magnetic flux, so that one has equidistant magnetic flux surfaces; improves convergence a lot

  integer :: iAzimOffset = 2 ! Equidistance sought for most problematic local time
  !C iAzimOffset = 1 ! Equidistance maintained at midnight

  integer :: isSORDetailNeeded = 0 ! No details for the inner SOR iterations 
  ! isSORDetailNeeded = 1 ! Details about inner SOR iterations

  ! isEnergDetailNeeded = 0
  integer :: isEnergDetailNeeded = 1 ! Dst computation (DPS formula with thermal energy inside domain)

  integer :: isFBDetailNeeded = 0 ! Computes global force imbalance

  ! iLossCone = 2 ! More realistic, empty loss cone for RAM computations (M.
  ! Liemohn's formalism, Liemohn, 2004); 
  ! can have some errors with very deformed B-field
  integer :: iLossCone = 1 ! Filled loss cone
  integer :: iOutput = 1 ! Output for Vania's RAM - hI.dat files are written to disk

  integer :: iTilt = 0     ! Zero dipole tilt
  ! iTilt = 1   ! Case with non-zero magnetic tilt

  integer :: iHourChoice = 0 ! More than 1 equilibrium calculation (for RAM-SCB)

  character(len=200) :: prefixIn     = trim(PathScbIn)  // '/'
  character(len=200) :: prefixOut    = trim(PathScbOut) // '/'
  character(len=200) :: prefixRAMOut = trim(PathRamOut) // '/'

END MODULE ModScbMain
