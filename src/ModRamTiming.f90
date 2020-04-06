!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamTiming
!    A module for tracking code efficiency and other timing metrics.

  use ModRamMain,   ONLY: DP, niter, PathRamOut
  use ModRamParams, ONLY: DoSaveRamSats

  use ModTimeConvert, ONLY: TimeType

  implicit none

  type(TimeType) :: TimeRamNow, TimeRamStart, TimeRamStop, TimeRamRealStart, TimeRamFinish

  real(DP) :: TimeRestart    = 0.0, &
              TimeRamElapsed = 0.0, &
              TOld           = 0.0, &
              TimeMax        = 0.0 ! Simulation max in seconds

  integer :: MaxIter     ! Simulation max iterations

  !!!!! TEMPORAL GRIDS
  real(DP) :: DTs          = 5.0,    &  ! Variable that stores the current time step (changes during run)
              DTsNext      = 5.0,    &  ! Variable that decides the next time step (changes during run)
              DTsMin       = 0.1,    &  ! Minimum time step the code will take (doesn't actually adhear to this)
              DTsMax       = 100.0,  &  ! Max time step the code will take (configurable in PARAM)
              DTsFramework = 1000.0, &  ! Variable that dictates SWMF time steps (changes during run)
              DTOutput     = 3600.0, &  ! How often certain outputs are written
              DT_hI        = 300.0,  &  ! How often the SCB calculations are done
              DT_bc        = 300.0,  &  ! How often the boundary fluxes are updated
              DTEfi        = 300.0,  &  ! How often the electric field is updated
              DTRestart    = 3600.0, &  ! How often restart files are written (configurable in PARAM)
              DtLogFile    = 60.0,   &  ! How often the log file is written to (configurable in PARAM)
              DtWriteSat   = 60.0,   &  ! How often satellite files are written to (configurable in PARAM)
              DtW_Pressure = 300.0,  &
              DtW_hI       = 300.0,  &
              DtW_EField   = 3600.0, &
              DtW_MAGxyz   = 300.0,  &
              DtW_2DFlux   = 3600.0, &
              DtW_Losses   = 3600.0
  real(DP) :: T, UTs
  real(DP) :: Efficiency = 0.0, SysTimeStart, SysTimeNow
  real(DP) :: dtPrintTiming = 300.0

  ! Variables for writing efficiency file:
  logical :: DoWriteEffFile = .true.
  real(DP) :: dtEffFile = 300.0

  character(len=*), parameter :: NameMod='ModRamTiming'
  character(len=100) :: NameEffFile

  contains
!==================================================================================================
  subroutine init_timing()

    use ModIoUnit,  ONLY: UNITTMP_

    implicit none

    character(len=100) :: NameFile
    integer :: iError

    character(len=*), parameter :: NameSub = NameMod // '::init_timing'
    integer :: t1, clock_rate = 100, clock_max = 100000
    !-------------------------------------------------------------------------
    !Set system time for beginning of simulation.
    !call cpu_time(SysTimeStart)
    call system_clock(t1,clock_rate,clock_max)
    SysTimeStart = t1/clock_rate

    ! Initialize efficiency file.
    write(NameEffFile, '(a,a,i8.8,a)') &
         trim(PathRamOut), '/efficiency_n', nIter, '.txt'

    ! Open file:
    open(unit=UNITTMP_, file=NameEffFile, status='REPLACE', &
         action='WRITE', iostat=iError)
    if(iError .ne. 0) call CON_stop(NameSub// 'cannot open file '//NameFile)

    ! Write header:
    write(UNITTMP_, '(a,f10.1)')'Run started at t=',TimeRestart
    write(UNITTMP_, *)'SysTime  RunTime  Efficiency=Run/Sys'
    write(UNITTMP_, *)'------------------------------------'
    close(UNITTMP_)

  end subroutine init_timing

  !===========================================================================
  subroutine do_timing()

    use ModIoUnit,  ONLY: UNITTMP_

    implicit none

    real(DP) :: CpuTimeNow
    integer :: t1, clock_rate, clock_max

    clock_rate = 1000
    clock_max = 100000

    ! Update system time.
    call system_clock(t1,clock_rate,clock_max)
    CpuTimeNow = t1/clock_rate
    SysTimeNow = CpuTimeNow - SysTimeStart

    ! Update timing metrics (only efficiency so far...)
    Efficiency = (TimeRamElapsed-TimeRestart)/SysTimeNow

    if (DoWriteEffFile.and.(abs(mod(TimeRamElapsed, dtEffFile)).le.1e-9)) then
       open(UNITTMP_,FILE=NameEffFile,ACCESS='APPEND')
       write(UNITTMP_, '(2(f12.2, 1x), f12.8)') SysTimeNow, TimeRamElapsed, Efficiency
       close(UNITTMP_)
    endif

    ! Write timing report.
    if(abs(mod(TimeRamElapsed, dtPrintTiming)) .le. 1e-9) then
       write(*,'(a, f12.2, a, f12.2, a, f10.6, a)') &
            'Simulated ', TimeRamElapsed, 's in ', SysTimeNow, &
            's (',Efficiency,'x real-time)'
    end if

  end subroutine do_timing

  !===========================================================================
  subroutine finalize_timing()

    implicit none

    real(DP) :: CpuTimeNow
    integer :: t1, clock_rate, clock_max
    !-------------------------------------------------------------------------
    clock_rate = 1000
    clock_max = 100000

    ! Update system time.
    !call cpu_time(CpuTimeNow)
    !SysTimeNow = CpuTimeNow - SysTimeStart
    call system_clock(t1,clock_rate,clock_max)
    CpuTimeNow = t1/clock_rate
    SysTimeNow = CpuTimeNow - SysTimeStart

    ! Update timing metrics (only efficiency so far...)
    Efficiency = (TimeRamElapsed-TimeRestart)/SysTimeNow

    write(*,'(a, f12.2, a, f12.2, a, f12.8, a)') &
         'Simulated ', TimeRamElapsed, 's in ', SysTimeNow, &
         's (',Efficiency,'X)'

  end subroutine finalize_timing

  !===========================================================================
  function max_output_timestep(TimeIn)
    ! Return the largest timestep RAM-SCB can take before it must stop
    ! and write output.  For example, at time=0s, if satellite files are
    ! written every 5 minutes and logfiles every 1 minute, the largest time
    ! step is one minute.  Similarly, at t=30s, DtMax=30s.  This function
    ! ensures that writing output is never skipped by taking large timesteps.
    ! Because RAM uses a time-splitting approach, the answer is divided by
    ! two (because each step moves forward in time by Dt twice.)

    implicit none

    ! Arguments:
    real(DP) :: max_output_timestep
    real(DP), intent(in) :: TimeIn

    real(DP) :: DtSatTemp, hI_temp, bc_temp, ef_temp, rt_temp

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='max_output_timestep'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    hI_temp = Dt_hI
    if (Dt_hI.le.1.0) hI_temp = 9999999.9
    bc_temp = Dt_bc
    if (Dt_bc.le.1.0) bc_temp = 9999999.9
    ef_temp = DtEfi
    if (DtEfi.le.1.0) ef_temp = 9999999.9
    rt_temp = DtRestart
    if (DtRestart.le.1.0) rt_temp = 9999999.9

    ! Only include sats if we are writing them.
    DtSatTemp=999999.9
    if(DoSaveRamSats) DtSatTemp=DtWriteSat

    ! Biggest timestep we can take is the smallest difference of the amount
    ! of time that passed since the file was last written and the write
    ! frequency.  Divide by two because of time splitting.
    max_output_timestep=min( &
         DtLogfile   -mod(TimeIn, DtLogfile   ), &
         DtSatTemp   -mod(TimeIn, DtSatTemp   ), &
         hI_temp     -mod(TimeIn, hI_temp     ), &
         bc_temp     -mod(TimeIn, bc_temp     ), &
         ef_temp     -mod(TimeIn, ef_temp     ), &
         rt_temp     -mod(TimeIn, rt_temp     ), &
         DtW_Pressure-mod(TimeIn, DtW_Pressure), &
         DtW_EField  -mod(TimeIn, DtW_EField  ), &
         DtW_hI      -mod(TimeIn, DtW_hI      ), &
         DtW_MAGxyz  -mod(TimeIn, DtW_MAGxyz  ), &
         DtW_2DFlux  -mod(TimeIn, DtW_2DFlux  ), &
         DtW_Losses  -mod(TimeIn, DtW_Losses  )) / 2.0

    if(DoTestMe)then
       write(*,'(2a,f11.2,a)')NameSub,' using these values at t=',TimeIn,':'
       write(*,'(a,f8.2)')'    DtLogfile=',DtLogfile
       write(*,'(a,f8.2)')'    DtSats   =',DtSatTemp
       write(*,'(a,f8.2)')'    DtRestart=',DtRestart
       write(*,'(a,f8.2)')'    DtScb    =',Dt_hI
       write(*,'(a,f8.2)')' MaxDtOutput=', max_output_timestep
    end if
  end function max_output_timestep

  !===========================================================================
end module ModRamTiming
!=============================================================================
