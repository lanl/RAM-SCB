!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModRamGrids

  use nrtype, ONLY: DP

  implicit none

  character(len=100) :: NameVar = '_H _O He _e'
 
  !!! RAM Grids
  integer :: nS       = 4,   &  ! number of species to run
             NR       = 20,  &  ! grid points in radial direction 
             NT       = 25,  &  ! grid points in local time direction
             NE       = 35,  &  ! number of energy bins
             NPA      = 72,  &  ! grid points in pitch angle dimension
             NEL      = 36,  &  ! Energy bins in boundary file (for SWMF or LANL)
             NEL_prot = 36,  &  ! Energy bins in bound file for protons (SWMF or LANL)
             NBD      = 288, &  ! Number of boundary data per file = (12/hour)*(24 hours)
             NTL      = 49      ! MLT bins for SWMF or LANL grid (25 for LANL, 49 for SWMF)
  integer, parameter :: NL       = 35,  &  ! grid points in radial for PLANE.F
                        NLT      = 48      ! grid points in local time for PLANE.F
  
  
  ! For B coupling, we must request additional ghostcells.
  integer :: nRextend ! = NR+3
  
  integer :: Slen = 55, &
             ENG  = 45, &
             NCF  = 5,  &
             NCO  = 5
  
  real(DP), parameter :: RadiusMin = 1.75, & ! Min radius for RAM grid
                         RadiusMax = 6.5,  & ! Max radius for RAM Grid
                         RadOut    = 6.75, & ! SCB CANDIDATE FOR REMOVAL
                         RadMaxScb = 9.00    ! Extended RAM grid for SCB use
  
  real(DP) :: EnergyMin = 0.1
  
  integer :: Nx = 72             ! Grid points for PA in mixed diffusion solver
  integer :: Ny = 37
  integer :: NyE = 300             ! Grid points for energy    "     "

  ! Grid derivs
  real(DP) :: dR, dPhi

  !EMIC diffusion coeff. 
  integer, parameter :: NCF_emic=10, ENG_emic=41, NPA_emic=89, NLZ_emic=5

end Module ModRamGrids
