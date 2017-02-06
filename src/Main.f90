program ram_scb

!============================================================================
!    RAM_SCB, Version 2.X
!
!    This is the main program unit for RAM_SCB; for compilation & execution,
!    please see accompanying documentation.
!
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

use ModMpi
use ModReadParam
use ModRamMpi
use ModRamMain, ONLY: IsComponent, TimeMax, TimeRamElapsed, iCal

implicit none

!----------------------------------------------------------------------------
! Ensure code is set to StandAlone mode.
IsComponent = .false.

! Initialize MPI.
call MPI_INIT(ierror)
iComm = MPI_COMM_WORLD

call MPI_COMM_RANK(iComm,iProc,ierror)
call MPI_COMM_SIZE(iComm,nProc,ierror)

! Read parameter files
call read_file('PARAM.in', iComm)
call read_init('  ', iSessionIn=1, iLineIn=0)
call IM_set_parameters

! Initialize RAM_SCB
call init_ramscb

! Run RAM-SCB
if(TimeRamElapsed .lt. TimeMax) then ! No wasted cycles, please.
   MAIN: do  
      call run_ramscb
      if(TimeRamElapsed .ge. TimeMax) exit MAIN
   end do MAIN
end if

call finalize_ramscb

call MPI_Finalize(iError)

end program ram_scb

!============================================================================
subroutine CON_stop(String)
  ! "Safely" stop RAM-SCB on all nodes.
  use ModRamMpi,  ONLY: iProc, iComm, iError
  use ModRamMain, ONLY: TimeRamElapsed
  use ModMpi
  implicit none

  character(len=*), intent(in) :: String
  write(*,*)'Stopping execution! me=',iProc,' at time=',TimeRamElapsed,&
       ' with msg:'
  write(*,*)String
  call MPI_abort(iComm, 0, iError)
  stop
end subroutine CON_stop

!============================================================================

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  ! Replaces the SWMF testing routine.
  
  use ModRamMain, ONLY: StringTest

  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest   = i_sub_string(' '//StringTest,' '//String//' ')>0
  DoTestMe = DoTest !.and. i_proc()==iProcTest
contains
  !===========================================================================
  integer function i_sub_string(StringA,StringB)

    ! This is needed to avoid some SGI f90 compiler bug 
    ! (which results in a memory leak) if we use
    !
    ! index(' '//StringTest,' '//str//' ')
    !
    ! directly.

    implicit none

    character (len=*), intent(in) :: StringA, StringB

    i_sub_string=index(StringA, StringB)

  end function i_sub_string

end subroutine CON_set_do_test

!============================================================================
