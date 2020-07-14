!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamParams

  use ModRamMain, ONLY: DP

  implicit none

 !!!! PARAM Variables
  logical :: DoUseRAM = .true.
  logical :: TimedRestarts = .true.
  logical :: Optim = .false.
  logical :: verbose = .false.
  logical :: Reset = .true.
  logical :: SCBonRAMTime = .false.
  integer :: RAMTie = 300
  logical :: checkMGNP = .false.
  logical :: WriteBoundary = .false.
  logical :: WritePotential = .false.
  logical :: integral_smooth = .true.
  logical :: InitializeOnFile = .true. 
  logical :: FixedComposition = .false.

  ! Use a fixed refilling time from newtau.dat file
  logical :: UseFixedTau = .false.

  ! Standalone or Component (default = standalone)  logical :: IsComponent = .false.
  logical :: IsComponent = .false.
 
  ! SHIELDS-RC Mode?
  logical :: IsSHIELDS=.false.
 
  ! RESTART variables.
  ! Is this run a restart of a previous simulation?
  logical :: IsRestart = .false., DoSaveFinalRestart=.true.
  logical :: HardRestart = .false. ! Whether to recompute SCB outer boundary
 
  ! Include SCB?  Default is yes!
  logical :: DoUseScb = .true.
 
  ! Include induced E field?  Default is no!
  logical :: UseEfind = .false.
 
  ! Include wave-particle diffusion?  Default is no!
  logical :: DoUseWPI = .false.
 
  ! Use BAS wave-particle diffusion?  Default is no!
  logical :: DoUseBASdiff = .false.
 
  ! Do use Kp-dependent diffusion coefficients?  Default is no!
  logical :: DoUseKpDiff = .false.
 
  ! Do use EMIC wave-particle diffusion?  Default is no!
  logical :: DoUseEMIC = .false.

  ! Do use FLC scattering diffusion?  Default is no!
  logical :: DoUseFLC = .false.
  logical :: DoWriteFLCDiffCoeff = .false.

  ! Flags for initialization of energy grid and initial conditions
  logical, public :: DoInitOnly = .false.
  logical :: DoUseVAPini = .false.
 
  ! Force multi-species boundary conditions through the use of 
  ! multiple files for BCS flux values.
  logical :: DoMultiBcsFile = .false.
 
  ! Dump 2D fluxes, B field, currents, pressure anis to netcdf?  Default is no!
  logical :: DoWriteFlux=.false.
  logical :: DoWriteB=.false.
  logical :: DoWriteCur=.false.
  logical :: DoWritePres=.false.
 
  ! Dump 2D plasmaspheric density to netcdf?  Default is no!
  logical :: DoWriteDensity=.false.
 
  ! Include electrons?  Default is yes!
  logical :: electrons = .true.
 
  logical :: DoAnisoPressureGMCoupling = .false.
 
  ! Limiter Beta (scales MC Limiter from min-mod=1 to superbee=2; default=1.5)
  real(DP) :: BetaLim = 1.5
 
  ! Whether or not to use boundary files
  logical :: BoundaryFiles = .true.
  real(DP) :: ElectronFluxCap = 1e10
  real(DP) :: ProtonFluxCap = 1e8

  ! Sets a fraction of oxygen to turn into nitrogen
  real(DP) :: OfracN = 0.
 
  logical :: DoVarDt = .true.                        ! Use variable timestep.
 
  character(len=100) :: StringTest=''                ! List of tested subroutines
 
  character(len=200) :: StrRamDescription='None'     ! Descript. of simulation
 
  logical :: IsStarttimeSet=.false.  ! True if #STARTTIME used in standalone.

  ! Plasmasphere Model Parameters 
  logical :: DoUsePlasmasphere = .false.
  character(len=100) :: PlasmasphereModel = 'Carpenter' ! Name of plasmasphere model to use
  character(len=100) :: TauCalculation = 'Analytic'     ! Name of method for calculating plasmasphere refilling time

  ! Coulomb Collision Parameters
  logical :: DoUseCoulomb = .false.

  character(len=4) :: NameBoundMag = 'DIPL'
  character(len=6) :: InnerMag, OuterMag
 
  character(len=6) :: event = ''
  character(len=4) :: boundary = 'LANL'
  character(len=4) :: electric = 'VOLS'
  character(len=4) :: NameDistrib = 'MAXW'
 
  character(len=200) :: NameIndexFile = 'RamIndices.txt'
  character(len=200) :: NameOmniFile  = 'omni.txt'
  character(len=200) :: BoundaryPath  = 'IM/input_ram/'
  character(len=200) :: InitializationPath = 'IM/input_ram/'
  character(len=200) :: HissFilePath = '/projects/lanl/SHIELDS/RAM_SCB_input/PlHiss/' 
  character(len=200) :: BASFilePath = '/projects/lanl/SHIELDS/shields_codes/BAS_bavDxx/'

  character(len=7) :: densityMode = "RAIRDEN"

  ! File write logicals:
  logical :: DoSaveRamSats=.false.
 
  ! Type of file name format, default will be to use new standard
  ! once the new standard is ready.
  logical :: UseNewFmt = .true.
 
  logical :: UseSWMFFile = .false.

end Module ModRamParams
!==============================================================================
