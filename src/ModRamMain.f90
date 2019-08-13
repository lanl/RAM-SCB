!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamMain

  use ModTimeConvert, ONLY: TimeType
  use nrtype,         ONLY: DP, SP
  use ModRamGrids,    ONLY: NS, NR, NT, NE, NPA

  implicit none

! RAM version.
  real, parameter :: CodeVersion = 2.2

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
  integer :: iCal  ! Tracks iteration number, excluding restart
  integer :: nIter ! Tracks iteration, including restart.
!!!!!

!  integer :: S

end Module ModRamMain
!==============================================================================
