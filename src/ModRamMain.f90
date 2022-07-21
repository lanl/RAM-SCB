!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamMain

  use ModTimeConvert, ONLY: TimeType, time_int_to_real
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
  character(len=999) :: PathRestartOut= 'IM/restartOUT'

  ! File paths to fix hardcoding of path issues.
  character(len=*), parameter :: PathRAMOut  = 'IM/output/ram/'
  character(len=*), parameter :: PathSCBOut  = 'IM/output/scb/'
  character(len=*), parameter :: PathSCEOut  = 'IM/output/sce/'
  character(len=*), parameter :: PathSWMFOut = 'IM/output_swmf/'
  character(len=*), parameter :: PathRAMIn   = 'IM/input_ram/'
  character(len=*), parameter :: PathSCBIn   = 'IM/input_scb/'
  character(len=*), parameter :: PathSCEIn   = 'IM/input_scb/'
!!!!!

!!!!! ITERATIONS
  integer :: iCal  ! Tracks iteration number, excluding restart
  integer :: nIter ! Tracks iteration, including restart.
!!!!!

!!!!! TESTING
  integer :: nTestPassed
  integer :: nTestRun
!  integer :: S

  contains

  subroutine make_time(yr, mon, day, hr, min, sec, frs, timeObj)
    integer, intent(in) :: yr, mon, day, hr, min, sec
    real(Real8_), intent(in) :: frs
    type(TimeType), intent(inout) :: timeObj

    timeObj%iYear = yr
    timeObj%iMonth = mon
    timeObj%iDay = day
    timeObj%iHour = hr
    timeObj%iMinute = min
    timeObj%iSecond = sec
    timeObj%FracSecond = frs
    call time_int_to_real(timeObj)
  end subroutine make_time

  subroutine test_neq_abs(val1, val2, atol, result)
    real(Real8_), intent(in) :: val1, val2, atol
    logical, intent(inout) :: result
    result = merge(.TRUE., .FALSE., abs(val1 - val2).gt.atol)
  end subroutine test_neq_abs

end Module ModRamMain
!==============================================================================
