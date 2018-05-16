!==============================================================================
module ModRamPl_Ne
!    This module replaces the variable declarations in rampl_ne.f
!    CAJ, EDIT: 2012-12-03: original version
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.

! *** put other modules here ***
  use ModRamGrids,     ONLY: NR, NT, NL, NLT
  use ModRamVariables, ONLY: NECR, VT, LZ, PHI, DL1, IR1, DPHI
!\
! Grid size parameters from ModRamMain
!/
!     NR  no. grids in radial direction
!     NT  no. grids in local time
!     NL  no. radial grids for plane_scb.f90
!     NLT no. local time grids for plane_scb.f90
!     NECR(NL,0:NLT)     ! plasmasphere Ne density
!     VT(NR+1,NT)        ! electric conv potential in kV
!     LZ(NR+1)           ! L shell array
!     Phi(NT)            ! Magnetic local time in radian

  implicit none
  save

!\
! New parameters for type, version, etc.
! These will make the true ModRamMain as
! modules are decentralized.
!/

! Pl_Ne version.
real, parameter :: CodeVersion = 1.0

! Kind for double precision, set here for portability.
integer, parameter :: Real8_ = selected_real_kind(14,200)
integer, parameter :: Real4_ = selected_real_kind( 6, 37)

! Include IRI 2012?  Default is yes!
logical :: UseIRI_2012 = .false.

! Use the latest NRLMSISE model?  Default is yes!
logical :: UseNRL_MSISE = .false.

! Use a non-dipole magnetic field from SCB .dat file
logical :: UseSCB_nondipole = .false.

! Use a fixed refilling time from newtau.dat file
logical :: UseFixedTau = .true.

!\
! Physical parameters
!/

real(kind=Real8_), parameter :: &
     RE = 6.371E6, &	   	! Earth's radius [m]
     RE_cm = 6.371E8, &	   	! Earth's radius [cm]
     RE_km = 6371, &	   	! Earth's radius [km]
     PI = 3.141592654, &
     FTPI = 4./3.*PI, &
     KB = 1.38065E-16           ! Boltzmann constant [cm^2 g s^-2 K^-1]

!\
! Grid parameters
!/

! Now we take these from ModRamVariables
!real(kind=Real8_), parameter :: &
!     DL1=0.25, &	       ! Grid size of L shell
!     IR1=DL1/0.25, &
!     DPhi=2.*PI/(NT-1)         ! Grid size for local time [rad]

!\
! New hI file parameters
!/
character(len=32) :: hI_file_name = 'hI_output_0113.dat'

integer, parameter :: &
     Nd = 500, &     ! no. records in file
     NLsh = 20, &    ! no. L shells in file
     Nmlt = 25, &    ! no. local times in file
     Npa = 72        ! no. PAs

end Module ModRamPl_Ne
!==============================================================================
