!==============================================================================
module ModRamMain
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModTimeConvert, ONLY: TimeType
  use nrtype,         ONLY: DP, SP
  use ModRamGrids,    ONLY: NS, NR, NT, NE, NPA

  implicit none
  save

! RAM version.
  real, parameter :: CodeVersion = 3.0

! Kind for double precision, set here for portability.
  integer, parameter :: Real8_ = DP
  integer, parameter :: Real4_ = SP

!!!!! PATHS
  ! File paths for restarts
  character(len=*), parameter :: PathRestartIn = 'IM/restartIN'
  character(len=*), parameter :: PathRestartOut= 'IM/restartOUT'

  ! File paths to fix hardcoding of path issues.
  character(len=*), parameter :: PathRAMOut  = 'IM/output/ram/'
  character(len=*), parameter :: PathSCBOut  = 'IM/output/scb/'
  character(len=*), parameter :: PathSCEOut  = 'IM/output/scb/'
  character(len=*), parameter :: PathSWMFOut = 'IM/output_swmf/'
  character(len=*), parameter :: PathRAMIn   = 'IM/input_ram/'
  character(len=*), parameter :: PathSCBIn   = 'IM/input_scb/'
  character(len=*), parameter :: PathSCEIn   = 'IM/input_scb/'
!!!!!

!!!!! ITERATIONS
  integer :: iCal     ! Tracks iteration number, kind of.
  integer :: nIter=0  ! Tracks iteration, including restart.
!!!!!

!\
! PARAMS Block.
!/
integer :: S

end Module ModRamMain
!==============================================================================
