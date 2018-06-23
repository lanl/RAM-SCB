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
use ModRamVariables, ONLY: Kp, F107, dBdt, dIdt, dIbndt
use ModScbVariables, ONLY: hICalc, SORFail
use ModScbParams,    ONLY: method

!!!! Module Subroutines and Functions
use ModRamGSL,       ONLY: GSL_Initialize
use ModRamInjection, ONLY: injection
use ModRamFunctions, ONLY: ram_sum_pressure, RamFileName
use ModRamIndices,   ONLY: get_indices
use ModRamTiming,    ONLY: max_output_timestep, init_timing, finalize_timing, do_timing
use ModRamIO,        ONLY: init_output, handle_output, ram_write_pressure
use ModRamRestart,   ONLY: write_restart
use ModRamInit,      ONLY: ram_allocate, ram_init, init_input, ram_deallocate
use ModRamRun,       ONLY: ram_run
use ModRamBoundary,  ONLY: get_boundary_flux
use ModRamEField,    ONLY: get_electric_field
use ModScbInit,      ONLY: scb_allocate, scb_init, scb_deallocate
use ModScbRun,       ONLY: scb_run
use ModScbIO,        ONLY: computational_domain
use ModSceInit,      ONLY: sce_allocate, sce_init, sce_deallocate
use ModRamScb,       ONLY: ramscb_allocate, computehI, ramscb_deallocate, compute3DFlux

!!!! External Modules (share/Library/src)
use ModReadParam
use CON_planet,      ONLY: set_planet_defaults
use CON_axes,        ONLY: init_axes, test_axes
use ModPlanetConst,  ONLY: init_planet_const
use ModTimeConvert,  ONLY: time_real_to_int

!!!! MPI Modules
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

! Allocate Arrays
call ram_allocate
call scb_allocate
call ramscb_allocate
call sce_allocate

! Initialize RAM_SCB
call ram_init
call scb_init
call sce_init
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
!!!!!!!!!! SET CURRENT TIME AND TIMESTEP
      UTs = TimeRamElapsed + TimeRamStart%iHour*3600.0
      ! Set the half-timestep based on CFL, SWMF, and Maximum allowed.
      ! Use CFL Number if #VARIABLEDT is set in PARAM file.
      if (DoVarDt) then
         DtEndMax   = (TimeMax-TimeRamElapsed)/2.0
         DtOutputMax = max_output_timestep(TimeRamElapsed)
         DTs = min(DTsNext,DTsmax,DtOutputMax,DtEndMax,DTsFramework)
         if (Kp.gt.6.0 .AND. DTs.gt.5.0) DTs = 5.0
      else if(abs(mod(UTs, Dt_hI)) .le. 1e-9) then
         DTs = 5.0
         if(Kp .ge. 5.0) DTs = min(DTsMin,DTs)
         if(Kp .gt. 6.0) DTs = 1.0   !1. or 0.25
      endif
      DTsNext = DTsmax
!!!!!!!!!!

!!!!!!!!!! UPDATES AS NEEDED
      call get_indices(TimeRamNow%Time, Kp, f107)

      ! Update Boundary Flux if Dt_bc has passed
      if (abs(mod(TimeRamElapsed, Dt_bc)).le.1e-9) then
         call get_boundary_flux
      end if

      ! Update Electric Fields if DtEfi has passed
      if (abs(mod(TimeRamElapsed, DtEfi)).le.1e-9) then
         call get_electric_field
      end if
!!!!!!!!!

!!!!!!!!!! RUN RAM
      ! Broadcast current call to ram_all
      call write_prefix
      write(*,*) 'Calling ram_run for UTs, DTs,Kp = ', UTs, Dts, Kp
      ! Call RAM for each species.
      if (DoUseRAM) call ram_run
      FLUSH(6)

      ! Increment and update time
      TimeRamElapsed = TimeRamElapsed + 2.0 * DTs
      TimeRamNow % Time = TimeRamStart % Time + TimeRamElapsed
      call time_real_to_int(TimeRamNow)
!!!!!!!!!

!!!!!!!!! RUN SCB
      ! Call SCB if conditions are met
      triggerSCB = .false.
      if (SCBonRAMTime) then
         if ((mod(nIter,RAMTie).eq.0).and.(nIter.gt.1)) triggerSCB = .true.
      else
         if (abs(mod(TimeRamElapsed, Dt_hI)).le.1e-9) triggerSCB = .true.
      endif
      if (triggerSCB) then
         call write_prefix
         write(*,*) 'Running SCB model to update B-field...'

         call ram_sum_pressure
         call scb_run(nIter)
         FLUSH(6)
         if ((SORFail).and.(Reset)) then
            if (verbose) print*, 'Error in SCB calculation, attempting a full reset'
            call computational_domain
            call scb_run(0)
            if (SORFail) hICalc = .false.
         endif
         FLUSH(6)

         ! Couple SCB -> RAM
         if ((hICalc).and.(method.ne.3)) then ! Calculate full h's and I's if SCB was successful
            call computehI(nIter)
         else                                 ! If SCB wasn't successful use previous h's and I's
            dBdt = 0._dp                      ! which implies dXdt = 0
            dIdt = 0._dp
            dIbndt = 0._dp
         endif
         call compute3DFlux
         FLUSH(6)

         call write_prefix
         write(*,*) 'Finished 3D Equilibrium code.'
         DtsNext = DtsMin ! This is so we don't accidently take to big of a step after an SCB call
      end if
!!!!!!!!!

!!!!!!!!! CREATE OUTPUTS
      ! Do timing.
      call do_timing
      ! Check and write output, as necessary.
      call handle_output(TimeRamElapsed)
      FLUSH(6)
!!!!!!!!!

!!!!!!!!! FINISH TIME STEP
      iCal  = iCal + 1
      nIter = nIter + 1
      if(TimeRamElapsed .ge. TimeMax) exit MAIN
!!!!!!!!!
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
  subroutine write_prefix

    use ModRamParams, ONLY: IsComponent


    implicit none

    character(len=7) :: StringPrefix = 'IM:'

    if(.not. IsComponent) RETURN
    write(*,'(a)',ADVANCE='NO')trim(StringPrefix)
   
  end subroutine write_prefix

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
