module ModRamTiming
!    A module for tracking code efficiency and other timing metrics.
!    This is NOT for handling simulation time.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.

  use ModTimeConvert
  use ModRamMpi,     ONLY: iProc
  use ModRamMain,    ONLY: Real8_, PathRamOut, nIter, &
       TimeRamElapsed, TimeRestart

  implicit none
  save

  real(kind=Real8_) :: Efficiency = 0.0, SysTimeStart, SysTimeNow
  real(kind=Real8_) :: dtPrintTiming = 300.0

  ! Variables for writing efficiency file:
  logical :: DoWriteEffFile = .true.
  integer :: iUnitEffFile
  real(kind=Real8_) :: dtEffFile = 300.0

  character(len=*), parameter :: NameMod='ModRamTiming'
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  subroutine init_timing()

    use ModIoUnit,  ONLY: io_unit_new

    character(len=100) :: NameFile
    integer :: iError

    character(len=*), parameter :: NameSub = NameMod // '::init_timing'
    !-------------------------------------------------------------------------
    if(iProc .ne. 0) return

    !Set system time for beginning of simulation.
    call cpu_time(SysTimeStart)

    ! Initialize efficiency file.
    write(NameFile, '(a,a,i8.8,a)') &
         trim(PathRamOut), '/efficiency_n', nIter, '.txt'

    ! Open file:
    iUnitEffFile = io_unit_new()
    open(unit=iUnitEffFile, file=NameFile, status='REPLACE', &
         action='WRITE', iostat=iError)
    if(iError .ne. 0) call CON_stop(NameSub// 'cannot open file '//NameFile)

    ! Write header:
    write(iUnitEffFile, '(a,f10.1)')'Run started at t=',TimeRestart
    write(iUnitEffFile, *)'SysTime  RunTime  Efficiency=Run/Sys'
    write(iUnitEffFile, *)'------------------------------------'

  end subroutine init_timing

  !===========================================================================
  subroutine do_timing()
    use ModRamMpi, ONLY: iProc
    use ModRamIO,  ONLY: write_prefix
    real(kind=Real8_) :: CpuTimeNow
    !-------------------------------------------------------------------------
    if(iProc .ne. 0) return

    ! Update system time.
    call cpu_time(CpuTimeNow)
    SysTimeNow = CpuTimeNow - SysTimeStart

    ! Update timing metrics (only efficiency so far...)
    Efficiency = (TimeRamElapsed-TimeRestart)/SysTimeNow

    if(DoWriteEffFile .and. (mod(TimeRamElapsed, dtEffFile) .eq. 0) ) &
         write(iUnitEffFile, '(2(f12.2, 1x), f12.8)') &
         SysTimeNow, TimeRamElapsed, Efficiency

    ! Write timing report.
    if(mod(TimeRamElapsed, dtPrintTiming) .eq. 0) then
       call write_prefix
       write(*,'(a, f12.2, a, f12.2, a, f10.6, a)') &
            'Simulated ', TimeRamElapsed, 's in ', SysTimeNow, &
            's (',Efficiency,'x real-time)'
    end if

  end subroutine do_timing

  !===========================================================================
  subroutine finalize_timing()
    
    use ModRamIO, ONLY: write_prefix

    real(kind=Real8_) :: CpuTimeNow
    !-------------------------------------------------------------------------
    close(iUnitEffFile)

    ! Update system time.
    call cpu_time(CpuTimeNow)
    SysTimeNow = CpuTimeNow - SysTimeStart

    ! Update timing metrics (only efficiency so far...)
    Efficiency = (TimeRamElapsed-TimeRestart)/SysTimeNow

    call write_prefix
    write(*,'(a, f12.2, a, f12.2, a, f12.8, a)') &
         'Simulated ', TimeRamElapsed, 's in ', SysTimeNow, &
         's (',Efficiency,'X)'

  end subroutine finalize_timing

  !===========================================================================
end module ModRamTiming
!=============================================================================
