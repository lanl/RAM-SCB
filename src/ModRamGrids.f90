!==============================================================================
Module ModRamGrids
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

implicit none
save

!!! RAM Grids
integer :: NR       = 20,  &  ! grid points in radial direction 
           NT       = 25,  &  ! grid points in local time direction
           NE       = 35,  &  ! number of energy bins
           NS       = 4,   &  ! number of species
           NPA      = 72,  &  ! grid points in pitch angle dimension
           NL       = 35,  &  ! grid points in radial for PLANE.F
           NLT      = 48,  &  ! grid points in local time for PLANE.F
           NEL      = 36,  &  ! Energy bins in boundary file (for SWMF or LANL)
           NEL_prot = 36,  &  ! Energy bins in bound file for protons (SWMF or LANL)
           NBD      = 288, &  ! Number of boundary data per file = (12/hour)*(24 hours)
           NTL      = 49      ! MLT bins for SWMF or LANL grid (25 for LANL, 49 for SWMF)

! For B coupling, we must request additional ghostcells.
integer :: nRextend! = NR+3

integer :: Slen = 55, &
           ENG  = 45, &
           NCF  = 5,  &
           NCO  = 5

real, parameter :: RadiusMin = 1.75, &
                   RadiusMax = 6.5

real :: EnergyMin = 0.1

integer :: Nx! = NPA             ! Grid points for PA in mixed diffusion solver
integer :: Ny = 37
integer :: NyE = 300             ! Grid points for energy    "     "

end Module ModRamGrids
