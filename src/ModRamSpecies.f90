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
     character(len=100) :: cross_sections ! File name for cross sections (or na)
     logical          :: FLC      ! Sets whether the species uses fieldline curvature scattering
     logical          :: EMIC     ! Sets whether the species uses EMIC waves
  end type SpeciesType

  integer, parameter :: nSpecies = 5
  type(SpeciesType), dimension(nSpecies) :: RAMSpecies

  contains

  subroutine DefineSpecies
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
     RAMSpecies(1)%cross_sections = 'na'
     RAMSpecies(1)%FLC = .false.
     RAMSpecies(1)%EMIC= .false.
     
     ! Protons
     RAMSpecies(2)%s_name = "Hydrogen"
     RAMSpecies(2)%s_code = "_H"
     RAMSpecies(2)%s_mass = 1.0
     RAMSpecies(2)%s_comp = 1.0
     RAMSpecies(2)%s_charge = 1
     RAMSpecies(2)%SCB = .true.
     RAMSpecies(2)%WPI = .false.
     RAMSpecies(2)%CEX = .true.
     RAMSpecies(2)%cross_sections = 'na'
     RAMSpecies(2)%FLC = .true.
     RAMSpecies(2)%EMIC= .true.
     
     ! Helium +1
     RAMSpecies(3)%s_name = "HeliumP1"
     RAMSpecies(3)%s_code = "He"
     RAMSpecies(3)%s_mass = 4.0
     RAMSpecies(3)%s_comp = 1.0
     RAMSpecies(3)%s_charge = 1
     RAMSpecies(3)%SCB = .true.
     RAMSpecies(3)%WPI = .false.
     RAMSpecies(3)%CEX = .true.
     RAMSpecies(3)%cross_sections = 'na'
     RAMSpecies(3)%FLC = .true.
     RAMSpecies(3)%EMIC= .true.
    
     ! Oxygen +1
     RAMSpecies(4)%s_name = "OxygenP1"
     RAMSpecies(4)%s_code = "_O"
     RAMSpecies(4)%s_mass = 16.0
     RAMSpecies(4)%s_comp = 1.0
     RAMSpecies(4)%s_charge = 1
     RAMSpecies(4)%SCB = .true.
     RAMSpecies(4)%WPI = .false.
     RAMSpecies(4)%CEX = .true.
     RAMSpecies(4)%cross_sections = 'na'
     RAMSpecies(4)%FLC = .true.
     RAMSpecies(4)%EMIC= .true.

     ! Nitrogen +1
     RAMSpecies(5)%s_name = "Nitrogen"
     RAMSpecies(5)%s_code = "_N"
     RAMSpecies(5)%s_mass = 14.0
     RAMSpecies(5)%s_comp = 1.0
     RAMSpecies(5)%s_charge = 1
     RAMSpecies(5)%SCB = .false.
     RAMSpecies(5)%WPI = .false.
     RAMSpecies(5)%CEX = .true.
     RAMSpecies(5)%cross_sections = 'NitrogenCrossSections.dat'
     RAMSpecies(5)%FLC = .false.
     RAMSpecies(5)%EMIC= .false.
     
  end subroutine DefineSpecies

END MODULE ModRamSpecies
