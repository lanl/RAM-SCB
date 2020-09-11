!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamSpecies
  use ModRamMain, ONLY: DP

  implicit none

  type SpeciesType
     character(len=8) :: s_name   ! Sets the species name
     character(len=2) :: s_code   ! Sets the species unique code
     real(DP)         :: s_mass   ! Sets the species mass as a multiple of the proton mass
     real(DP)         :: s_comp   ! Sets the species composition ratio
     integer          :: s_charge ! Sets the charge state of the species
     logical          :: SCB      ! Sets whether the species is used in the pressure balance
     logical          :: WPI      ! Sets whether the species uses wave particle interactions
     logical          :: CEX      ! Sets whether the species charge exchanges
     logical          :: FLC      ! Sets whether the species uses fieldline curvature scattering
     logical          :: EMIC     ! Sets whether the species uses EMIC waves
     character(len=100) :: Initialization
     character(len=100) :: CEX_file ! File name for cross sections (or na)
     character(len=100) :: CEX_species ! Species to charge exchange against (nH nO nN) (or na)
     real(DP), allocatable :: CEX_velocities(:)
     real(DP), allocatable :: CEX_nH(:)
     real(DP), allocatable :: CEX_nO(:)
     real(DP), allocatable :: CEX_nN(:)
     real(DP) :: plasmasphereRatio
  end type SpeciesType

  integer, parameter :: nSpecies = 6
  type(SpeciesType), dimension(nSpecies) :: RAMSpecies

  contains

  subroutine DefineSpecies
     ! Define the available species for use in RAM-SCB

     implicit none

     ! Electrons
     RAMSpecies(1)%s_name = "Electron"
     RAMSpecies(1)%s_code = "_e"
     RAMSpecies(1)%s_mass = 5.4462E-4
     RAMSpecies(1)%s_comp = 1.0
     RAMSpecies(1)%s_charge = -1
     RAMSpecies(1)%SCB = .true.
     RAMSpecies(1)%WPI = .true.
     RAMSpecies(1)%CEX = .false.
     RAMSpecies(1)%FLC = .false.
     RAMSpecies(1)%EMIC = .false.
     RAMSpecies(1)%CEX_file = 'na'
     RAMSpecies(1)%CEX_species = 'na'
     RAMSpecies(1)%Initialization = 'InitializationFile'
     RAMSpecies(1)%plasmasphereRatio = 1.0

     ! Protons
     RAMSpecies(2)%s_name = "Hydrogen"
     RAMSpecies(2)%s_code = "_H"
     RAMSpecies(2)%s_mass = 1.0
     RAMSpecies(2)%s_comp = 1.0
     RAMSpecies(2)%s_charge = 1
     RAMSpecies(2)%SCB = .true.
     RAMSpecies(2)%WPI = .false.
     RAMSpecies(2)%CEX = .true.
     RAMSpecies(2)%FLC = .true.
     RAMSpecies(2)%EMIC = .true.
     RAMSpecies(2)%CEX_file = 'na'
     RAMSpecies(2)%CEX_species = 'na'
     RAMSpecies(2)%Initialization = 'InitializationFile'
     RAMSpecies(2)%plasmasphereRatio = 0.77

     ! Helium +1
     RAMSpecies(3)%s_name = "HeliumP1"
     RAMSpecies(3)%s_code = "He"
     RAMSpecies(3)%s_mass = 4.0
     RAMSpecies(3)%s_comp = 1.0
     RAMSpecies(3)%s_charge = 1
     RAMSpecies(3)%SCB = .true.
     RAMSpecies(3)%WPI = .false.
     RAMSpecies(3)%CEX = .true.
     RAMSpecies(3)%FLC = .true.
     RAMSpecies(3)%EMIC = .true.
     RAMSpecies(3)%CEX_file = 'na'
     RAMSpecies(3)%CEX_species = 'na'
     RAMSpecies(3)%Initialization = 'InitializationFile'
     RAMSpecies(3)%plasmasphereRatio = 0.2

     ! Oxygen +1
     RAMSpecies(4)%s_name = "OxygenP1"
     RAMSpecies(4)%s_code = "_O"
     RAMSpecies(4)%s_mass = 16.0
     RAMSpecies(4)%s_comp = 1.0
     RAMSpecies(4)%s_charge = 1
     RAMSpecies(4)%SCB = .true.
     RAMSpecies(4)%WPI = .false.
     RAMSpecies(4)%CEX = .true.
     RAMSpecies(4)%FLC = .true.
     RAMSpecies(4)%EMIC = .true.
     RAMSpecies(4)%CEX_file = 'na'
     RAMSpecies(4)%CEX_species = 'na'
     RAMSpecies(4)%Initialization = 'InitializationFile'
     RAMSpecies(4)%plasmasphereRatio = 0.03

     ! Nitrogen +1
     RAMSpecies(5)%s_name = "Nitrogen"
     RAMSpecies(5)%s_code = "_N"
     RAMSpecies(5)%s_mass = 14.0
     RAMSpecies(5)%s_comp = 1.0
     RAMSpecies(5)%s_charge = 1
     RAMSpecies(5)%SCB = .false.
     RAMSpecies(5)%WPI = .false.
     RAMSpecies(5)%CEX = .true.
     RAMSpecies(5)%FLC = .false.
     RAMSpecies(5)%EMIC = .false.
     RAMSpecies(5)%CEX_file = 'NitrogenCrossSections.dat'
     RAMSpecies(5)%CEX_species = 'nH'
     RAMSpecies(5)%Initialization = 'na'
     RAMSpecies(5)%plasmasphereRatio = 0.0

     ! Strontium +1
     RAMSpecies(6)%s_name = "Strontium"
     RAMSpecies(6)%s_code = "Sr"
     RAMSpecies(6)%s_mass = 87.62
     RAMSpecies(6)%s_comp = 0.0
     RAMSpecies(6)%s_charge = 1
     RAMSpecies(6)%SCB = .false.
     RAMSpecies(6)%WPI = .false.
     RAMSpecies(6)%CEX = .true.
     RAMSpecies(6)%FLC = .false.
     RAMSpecies(6)%EMIC = .false.
     RAMSpecies(6)%CEX_file = 'StrontiumCrossSections.dat'
     RAMSpecies(6)%CEX_species = 'nH nO nN'
     RAMSpecies(6)%Initialization = 'StrontiumPlusOneIons.dat'
     RAMSpecies(6)%plasmasphereRatio = 0.0

  end subroutine DefineSpecies

END MODULE ModRamSpecies
