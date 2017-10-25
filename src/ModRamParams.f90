!==============================================================================
module ModRamParams
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModRamMain, ONLY: Real8_
  implicit none
  save

!!!! PARAM Variables
 ! Standalone or Component (default = standalone)  logical :: IsComponent = .false.
 logical :: IsComponent = .false.

 ! SHIELDS-RC Mode?
 logical :: IsSHIELDS=.false.

 ! RESTART variables.
 ! Is this run a restart of a previous simulation?
 logical :: IsRestart = .false., DoSaveFinalRestart=.true.

 ! Use plasmasphere density from Rasmussen model, coupled to MSIS and IRI
 ! For plane_scb options, see ModRamPl_Ne.f90
 logical :: DoUsePlane_SCB = .false.

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
 real(kind=Real8_) :: BetaLim = 1.5

 logical :: DoVarDt = .true.                        ! Use variable timestep.

 character(len=100) :: StringTest=''                ! List of tested subroutines

 character(len=200) :: StrRamDescription='None'     ! Descript. of simulation

 logical :: IsStarttimeSet=.false.  ! True if #STARTTIME used in standalone.

 character(len=4) :: NameBoundMag

 character(len=6) :: event = 'de_flt'
 character(len=4) :: boundary, electric, NameDistrib

 character(len=200) :: NameIndexFile = 'RamIndices.txt'
 character(len=200) :: NameOmniFile  = 'omni.txt'

 ! File write logicals:
 logical :: DoSaveRamSats=.false.

 ! Type of file name format, default will be to use new standard
 ! once the new standard is ready.
 logical :: UseNewFmt = .true.

 logical :: UseSWMFFile = .false.

end Module ModRamParams
!==============================================================================
