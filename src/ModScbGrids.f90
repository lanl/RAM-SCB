!==============================================================================
Module ModScbGrids
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

use nrtype, ONLY: DP, pi_d

implicit none
save

integer, parameter :: NR = 20, NT = 25, NPA = 72
!!! SCB Grids
! Must match the computational boundary input file (*.cdf)
integer, parameter :: nthe   = 51,        &  !  nthe is the number of theta grid points
                      nthem  = nthe - 1,  &
                      nthep  = nthe + 1,  &
                      npsi   = 35,        &  !  npsi is the number of psi surfaces
                      npsim  = npsi - 1,  &
                      npsip  = npsi + 1,  &
                      nzeta  = 61,        &  !  nzeta is the number of zeta surfaces 
                      nzetam = nzeta - 1, &
                      nzetap = nzeta + 1, &
                      nq     = nthe + 2,  &
                      ny     = npsi + 2,  &
                      na     = nzeta + 2, &
                      nm     = 2 * nq,    &
                      nyp2   = ny + 2

! Various quantities for the finite difference schemes
real(DP), parameter ::  dr      = 1._dp/REAL(npsi - 1, dp),     &   ! d(rho)
             dt      = pi_d/REAL(nthe - 1, dp),      &   ! d(theta)
             dpPrime = 2*pi_d / REAL(nzeta -1, dp), &   ! d(zeta)
             rdr     = 1._dp / dr,                   &
             rdt     = 1._dp / dt,                   &
             rdp     = 1._dp / dpPrime,              &
             rdrsq   = rdr**2,                       &
             rdtsq   = rdt**2,                       &
             rdpsq   = rdp**2,                       &
             rdr2    = 0.5_dp * rdr,                 &
             rdt2    = 0.5_dp * rdt,                 &
             rdp2    = 0.5_dp * rdp,                 &
             rdr4    = 0.25_dp * rdr,                &
             rdt4    = 0.25_dp * rdt,                &
             rdp4    = 0.25_dp * rdp,                &
             rdrdp4  = 0.25_dp * rdr * rdp,          &
             rdpdt4  = 0.25_dp * rdp * rdt,          &
             rdtdr4  = 0.25_dp * rdt * rdr


! Number of LTs in RAM for coupling in RAM-SCB
integer, parameter :: nAzimRAM = NT

! Number of radial positions in RAM for coupling in RAM-SCB
integer, parameter :: nXRaw = NR-1, &
                      nYRaw = nAzimRAM

end Module ModScbGrids
