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
!!!! Module Variables
use ModRamMain,      ONLY: DP, iCal, nIter
use ModRamParams,    ONLY: DoSaveFinalRestart, DoVarDt, IsComponent, &
                           verbose, reset, DoUseRam, SCBonRAMTime, RAMTie
use ModRamTiming,    ONLY: DtsFramework, DtsMax, DtsMin, DtsNext, Dts, Dt_hI, &
                           TimeRamStart, TimeRamNow, TimeMax, TimeRamElapsed, &
                           UTs, Dt_bc, DtEfi
use ModRamVariables, ONLY: species
!!!! Module Subroutines and Functions
use ModRamSpecies,   ONLY: DefineSpecies
use ModRamScbRun,    ONLY: run_ramscb
use ModRamGSL,       ONLY: GSL_Initialize
use ModRamTiming,    ONLY: init_timing, finalize_timing
use ModRamIO,        ONLY: init_output
use ModRamRestart,   ONLY: write_restart
use ModRamInit,      ONLY: ram_allocate, ram_init, init_input, ram_deallocate
use ModScbInit,      ONLY: scb_allocate, scb_init, scb_deallocate
use ModSceInit,      ONLY: sce_allocate, sce_init, sce_deallocate
use ModRamScb,       ONLY: ramscb_allocate, ramscb_deallocate

!!!! External Modules (share/Library/src)
use ModReadParam
use CON_planet,      ONLY: set_planet_defaults
use CON_axes,        ONLY: init_axes, test_axes
use ModPlanetConst,  ONLY: init_planet_const
use ModTimeConvert,  ONLY: time_real_to_int

!!!! MPI Modules (needed for coupling with SWMF)
use ModMpi
use ModRamMpi


implicit none

logical :: triggerSCB
real(DP) :: DtOutputMax, DtEndMax
!----------------------------------------------------------------------------
! Ensure code is set to StandAlone mode.
IsComponent = .false.

!!!!!!!!!! START RAM-SCB INITIALIZATION !!!!!!!!!!!\
! Initialize MPI.
call MPI_INIT(ierror)
iComm = MPI_COMM_WORLD
call MPI_COMM_RANK(iComm,iProc,ierror)
call MPI_COMM_SIZE(iComm,nProc,ierror)

! Begin tracking code efficiency.
call init_timing()

! Read and set parameters
call read_file('PARAM.in', iComm)
call read_init('  ', iSessionIn=1, iLineIn=0)
call IM_set_parameters
call DefineSpecies

! Allocate Arrays
call ram_allocate
call scb_allocate
call ramscb_allocate
call sce_allocate

! Initialize RAM_SCBE
call ram_init
call scb_init
call sce_init

! Initialize GSL
call GSL_Initialize

! Initialize planet and axes if in stand-alone mode.
! In component mode, the framework takes care of this.
! This is needed for correct coordinate transformations.
if (.not.IsComponent) then
   call init_planet_const
   call set_planet_defaults
   call init_axes(TimeRamStart % Time)
end if

! Get initial Data (checks for restart)
call init_input

! Initialize output files as necessary.
call init_output

iCal = 1
!!!!!!!!!!! END RAM-SCB INITIALIZATION !!!!!!!!!!!

!!!!!!!!!!! START RAM-SCB RUN !!!!!!!!!!!!
if (TimeRamElapsed .lt. TimeMax) then ! No wasted cycles, please.
   MAIN: do
      call run_ramscb
      if(TimeRamElapsed .ge. TimeMax) exit MAIN
   end do MAIN
end if

if ((.not.IsComponent)) then
     write(*,*) &
          '==============================================================='
     write(*,*) '                     RAM_SCB finished!'

     write(*,'(f12.2,a,i10.10,a)') &
          TimeRamElapsed, 's simulated in ', iCal-1, ' iterations.'

     if(DoSaveFinalRestart) call write_restart

     ! Stop timing Ram; write final timing report.
     call finalize_timing

     ! Deallocate arrays
     call ram_deallocate
     call scb_deallocate
     call ramscb_deallocate
     call sce_deallocate

     write(*,*) &
          '==============================================================='
end if

call MPI_Finalize(iError)

end program ram_scb
!==================================================================================================
subroutine CON_stop(String)
  ! "Safely" stop RAM-SCB on all nodes.
  use ModScbGrids,     ONLY: nthe, npsi, nzeta
  use ModScbVariables, ONLY: x,y,z
  use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamNow
  use ModRamFunctions, ONLY: RamFileName

  use ModIOUnit,       ONLY: UNITTMP_


  implicit none

  character(len=*), intent(in) :: String
  character(len=200) :: FileName

  integer :: i, j, k

  write(*,*)'Stopping execution! at time=',TimeRamElapsed,' with msg:'
  write(*,*)String

  FileName = RamFileName('MAGxyz2','dat',TimeRamNow)
  open(UNITTMP_, File=FileName)
  write(UNITTMP_,*) nthe, npsi, nzeta
  do i=1,nthe
     do j = 1, npsi
        do k=1,nzeta
           write(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k)
        enddo
    enddo
  enddo
  close(UNITTMP_)

  stop
end subroutine CON_stop

!==================================================================================================

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  ! Replaces the SWMF testing routine.
  
  use ModRamParams, ONLY: StringTest


  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest   = i_sub_string(' '//StringTest,' '//String//' ')>0
  DoTestMe = DoTest !.and. i_proc()==iProcTest
contains
  !================================================================================================
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

!==================================================================================================
