!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModScbGrids

  use nrtype, ONLY: DP, pi_d
  

  implicit none
  
  !!! SCB Grids
  ! Must match the computational boundary input file (*.cdf)
  integer :: nthe   = 71, &  !  nthe is the number of theta grid points
             npsi   = 45, &  !  npsi is the number of psi surfaces
             nzeta  = 73     !  nzeta is the number of zeta surfaces 
  
  integer :: nthem, nthep, npsim, npsip, nzetam, nzetap, nq, ny, na, nm, nyp2
  
  ! Various quantities for the finite difference schemes
  real(DP) ::  dr, dt, dpPrime, rdr, rdt, rdp, rdrsq, rdtsq, rdpsq, rdr2, rdt2, &
               rdp2, rdr4, rdt4, rdp4, rdrdp4, rdpdt4, rdtdr4
  
  
  ! Number of LTs in RAM for coupling in RAM-SCB
  integer :: nAzimRAM! = NT
  
  ! Number of radial positions in RAM for coupling in RAM-SCB
  integer :: nXRaw! = NR-1
  integer :: nXRawExt! = NR+1
  integer :: nYRaw! = nAzimRAM

end Module ModScbGrids
