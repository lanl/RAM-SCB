!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModRamGrids

implicit none
save

!!! RAM Grids
integer :: NR       = 20,  &  ! grid points in radial direction 
           NT       = 49,  &  ! grid points in local time direction
           NE       = 35,  &  ! number of energy bins
           NS       = 4,   &  ! number of species
           NPA      = 72,  &  ! grid points in pitch angle dimension
           NEL      = 36,  &  ! Energy bins in boundary file (for SWMF or LANL)
           NEL_prot = 36,  &  ! Energy bins in bound file for protons (SWMF or LANL)
           NBD      = 288, &  ! Number of boundary data per file = (12/hour)*(24 hours)
           NTL      = 49      ! MLT bins for SWMF or LANL grid (25 for LANL, 49 for SWMF)
integer, parameter :: NL       = 35,  &  ! grid points in radial for PLANE.F
                      NLT      = 48      ! grid points in local time for PLANE.F


! For B coupling, we must request additional ghostcells.
integer :: nRextend! = NR+3

integer :: Slen = 55, &
           ENG  = 45, &
           NCF  = 5,  &
           NCO  = 5

real, parameter :: RadiusMin = 1.75, & ! Min radius for RAM grid
                   RadiusMax = 6.5, &  ! Max radius for RAM Grid
                   RadOut = 6.75       ! Extended RAM grid for SCB use

real :: EnergyMin = 0.1

integer :: Nx! = NPA             ! Grid points for PA in mixed diffusion solver
integer :: Ny = 37
integer :: NyE = 300             ! Grid points for energy    "     "

end Module ModRamGrids
