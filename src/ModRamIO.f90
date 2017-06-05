!============================================================================
module ModRamIO
!   A set of variables and tools for general I/O handling, file name
!   creation, etc.
!   Copyright (c) 2016, Los Alamos National Security, LLC
!   All rights reserved.
!******************************************************************************

  use ModRamMain, ONLY: Real8_

  implicit none
  save

  logical :: IsFramework = .false. 
  character(len=7) :: StringPrefix = 'IM:'

  ! File output names and units
  integer :: iUnitLog ! Logfile IO unit
  character(len=100) :: NameFileLog, NameFileDst

  ! File output frequencies:
  real(kind=Real8_) :: DtLogfile=60.0, DtWriteSat =  60.0

  ! File write logicals:
  logical :: DoSaveRamSats=.false.

  ! Type of file name format, default will be to use new standard
  ! once the new standard is ready.
  logical :: UseNewFmt = .false.

  ! String that contains the date and time for which this simulation
  ! was initialized (used for output metadata.)
  character(len=21) :: StringRunDate
  
contains
  !===========================================================================
  function max_output_timestep(TimeIn)
    ! Return the largest timestep RAM-SCB can take before it must stop
    ! and write output.  For example, at time=0s, if satellite files are
    ! written every 5 minutes and logfiles every 1 minute, the largest time
    ! step is one minute.  Similarly, at t=30s, DtMax=30s.  This function
    ! ensures that writing output is never skipped by taking large timesteps.
    ! Because RAM uses a time-splitting approach, the answer is divided by
    ! two (because each step moves forward in time by Dt twice.)
    use ModRamMain, ONLY: Real8_, Dt_hI, DtRestart

    ! Arguments:
    real(kind=Real8_) :: max_output_timestep
    real(kind=Real8_), intent(in) :: TimeIn
    
    real(kind=Real8_) :: DtSatTemp=999999.9

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='max_output_timestep'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Only include sats if we are writing them.
    if(DoSaveRamSats) DtSatTemp=DtWriteSat
    
    ! Biggest timestep we can take is the smallest difference of the amount
    ! of time that passed since the file was last written and the write
    ! frequency.  Divide by two because of time splitting.
    max_output_timestep=min( &
         DtLogfile-mod(TimeIn, DtLogfile), &
         DtSatTemp-mod(TimeIn, DtSatTemp), &
         Dt_hI    -mod(TimeIn, Dt_hI    ), &
         DtRestart-mod(TimeIn, DtRestart)  ) / 2.0

    if(DoTestMe)then
       call write_prefix
       write(*,'(2a,f11.2,a)')NameSub,' using these values at t=',TimeIn,':'
       write(*,'(a,f8.2)')'    DtLogfile=',DtLogfile
       write(*,'(a,f8.2)')'    DtSats   =',DtSatTemp
       write(*,'(a,f8.2)')'    DtRestart=',DtRestart
       write(*,'(a,f8.2)')'    DtScb    =',Dt_hI
       write(*,'(a,f8.2)')' MaxDtOutput=', max_output_timestep
    end if
  end function max_output_timestep

  !===========================================================================
  subroutine write_prefix

    use ModRamMain, ONLY: IsComponent

    if(.not. IsComponent) RETURN
    write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

  !==========================================================================
  function RamFileName(PrefixIn, SuffixIn, TimeIn)
    ! Create a file name using the new RAM-SCB I/O filename standards:
    ! FileNameOut = PrefixIn_dYYYYMMDD_tHHMMSS.SuffixIn
    ! PrefixIn and SuffixIn are short, descriptive strings (e.g. 
    ! 'pressure', 'ram_o' for prefixes and 'dat', 'cdf' for suffixes.

    use ModTimeConvert

    implicit none

    character(len=200)           :: RamFileName
    type(TimeType),   intent(in) :: TimeIn
    character(len=*), intent(in) :: PrefixIn, SuffixIn
    !------------------------------------------------------------------------
    write(RamFileName, "(a,'_d',i4.4,2i2.2,'_t',3i2.2,'.',a)") &
         trim(PrefixIn), &
         TimeIn % iYear, &
         TimeIn % iMonth, &
         TimeIn % iDay, &
         TimeIn % iHour, &
         TimeIn % iMinute, &
         TimeIn % iSecond, &
         trim(SuffixIn)

  end function RamFileName
  
  !==========================================================================
  subroutine gen_old_timetags(TimeStart, TimeNow, nDay, nFiveTotal, &
       nHourTotal, nFiveDay)
    ! In older versions of RAM, file names were built using many different
    ! time dependent strings.  This function generates those strings if the
    ! old file naming scheme is used.
    !
    ! Time(Start, Now) -- Time structures of the start time of the simulation
    !             and of the current time.  See ModTimeConvert for more info.
    ! nDay --     Number of days elapsed (starting from 0)
    !             Used to sort files into subfolders called Day01, etc.
    !             Note that if the simulation starts at 23:00UT, Day00
    !             ends 1 hour into the simulation...
    ! nFiveDay -- Number of five-minute intervals that has passed in the
    !             current day.
    ! nFiveTotal--Time elapsed normalized by 5 minutes.
    !             Starts from zero, used to order in/out files.
    ! nHourTotal--Time elapsed normalized by hours.  Starts from zero.
    !
    ! *CJ mod*: changed nFiveTotal from 5-minutes to DT_bc
    ! *CJ mod*: changed nHourTotal from 600-minutes to DtOutput
    
    use ModTimeConvert
    use ModRamMain, ONLY: Dt_bc, DtOutput, IsRestart
    
    type(TimeType), intent(in)     :: TimeStart, TimeNow
    integer, intent(out), optional :: nDay
    integer, intent(out), optional :: nFiveDay
    integer, intent(out), optional :: nFiveTotal
    integer, intent(out), optional :: nHourTotal
    !------------------------------------------------------------------------
    ! Current number of days progressed in simulation:
    if(present(nDay)) nDay = floor( (TimeNow % Time - TimeStart % Time + &
         TimeStart % iHour*3600.0) / 86400.0 )
    if(present(nFiveTotal)) then
          nFiveTotal = &
               floor( (TimeNow % Time - TimeStart % Time) / Dt_bc )
    endif
    if(present(nHourTotal)) nHourTotal = &
         floor( (TimeNow%Time - TimeStart%Time)/ DtOutput )
    if(present(nFiveDay)) nFiveDay     = &
         12 * TimeNow % iHour + TimeNow % iMinute / 5


  end subroutine gen_old_timetags

  !==========================================================================
  subroutine ncdf_check(iStatusIn, NameSubIn, NameFileIn)
    ! Check the status flag of a netcdf operation, CON_stop as necessary.

    use netcdf
    
    integer, intent(in) :: iStatusIn
    character(len=100), intent(in) :: NameSubIn
    character(len=200), intent(in), optional :: NameFileIn

    character(len=*), parameter :: NameSub = 'ncdf_check'
    !------------------------------------------------------------------------
    if(iStatusIn /= nf90_noerr) then
       write(*,*) 'ERROR WITH NETCDF:'
       write(*,*) '      Subroutine: ', trim(NameSubIn)
       if(present(NameFileIn)) &
            write(*,*) '      File: ', trim(NameFileIn)
       write(*,*) '      Error code: ', iStatusIn
       write(*,*) '      ', trim(nf90_strerror(iStatusIn))
       call CON_stop(NameSub // ': NetCDF error.')
    end if

  end subroutine ncdf_check

  !==========================================================================
  subroutine write_ncdf_globatts(iFileID)
    ! For a netcdf file that is still in define mode, write a plethora of 
    ! global attributes that all RAM-SCB output NetCDFs should have.

    use netcdf
    use ModRamMain, ONLY: StrRamDescription, TimeRamStart

    integer, intent(in) :: iFileID  ! NetCDF file ID number for open file.

    integer :: iStatus
    character(len=100) :: NameSub = 'write_ncdf_globatts'
    !------------------------------------------------------------------------
    ! Description of run:
    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'description', &
         StrRamDescription)
    call ncdf_check(iStatus, NameSub)

    ! Time of run:
    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'run_date', &
         StringRunDate)
    call ncdf_check(iStatus, NameSub)

    ! Start time of run:
    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'start_time', &
         TimeRamStart % String)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_ncdf_globatts

  !==========================================================================
  subroutine write_restart

    use ModIoUnit, ONLY: UnitTMP_
    use ModRamMain, ONLY: TimeRamStart, TimeRamElapsed, nIter, PathRestartOut, &
         f2, nR, nT, nE, nPA, PParT, PPerT
    
    implicit none
    
    integer                       :: s
    character(len=2), dimension(4):: NameSpecies = (/'e_','h_','he','o_'/)
    character(len=100)            :: NameFile

    character(len=*), parameter :: NameSub='write_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    if(DoTest)write(*,'(a,f11.2)') 'RAM-SCB: Writing restarts at t=',&
         TimeRamElapsed

    ! Write ascii portion of restart.
    NameFile=PathRestartOut//'/restart_info.txt'
    open(unit=UnitTMP_, file=trim(NameFile), status='replace')
    write(UnitTMP_, *) 'TIMING:'
    write(UnitTMP_, '(a, i4.4, 2i2.2, 1x, 3i2.2)')'Start (YYYYMMDD HHMMSS)= ', &
         TimeRamStart%iYear, TimeRamStart%iMonth, TimeRamStart%iDay, &
         TimeRamStart%iHour, TimeRamStart%iMinute,TimeRamStart%iSecond
    write(UnitTMP_, '(a,f15.4)')'Elapsed (seconds)      = ', TimeRamElapsed
    write(UnitTMP_, '(a, i15)') 'Completed Iterations   = ', nIter
    write(UnitTMP_, *)'GRID:'
    write(UnitTMP_, '(a, 4i3)') 'nR, nL, nE, nPA        = ', nR, nT, nE, nPA
    close(unitTMP_)

    ! Write binary portions of restart.
    do s=1, 4
       NameFile=PathRestartOut//'/restart_'//NameSpecies(s)//'.rst'
       if(DoTest) then
          call write_prefix
          write(*,*) 'Restart file for ', NameSpecies(s), ' = ', NameFile
       endif
       open(unit=UnitTMP_, file=trim(NameFile), status='replace', &
            form='unformatted')
       write(UnitTMP_) F2(s,:,:,:,:)
       close(UnitTMP_)
    end do

    NameFile=PathRestartOut//'/restart_ppar.rst'
    if(DoTest) then
       call write_prefix
       write(*,*) 'Restart file for parallel pressure', ' = ', NameFile
    endif
    open(unit=UnitTMP_, file=trim(NameFile), status='replace', &
         form='unformatted')
    write(UnitTMP_) PParT(:,:,:)
    close(UnitTMP_)

    NameFile=PathRestartOut//'/restart_pper.rst'
    if(DoTest) then
       call write_prefix
       write(*,*) 'Restart file for perpendicular pressure', ' = ', NameFile
    endif
    open(unit=UnitTMP_, file=trim(NameFile), status='replace', &
        form='unformatted')
    write(UnitTMP_) PPerT(:,:,:)
    close(UnitTMP_)

  end subroutine write_restart

  !==========================================================================
  subroutine read_restart
    
    use ModIoUnit, ONLY: UnitTMP_
    use ModTimeConvert, ONLY: time_int_to_real, time_real_to_int
    use ModRamMain, ONLY: TimeRamStart, TimeRamElapsed, nIter, PathRestartIn, &
         f2, nR, nT, nE, nPA, Real8_, TimeRamNow, TimeRestart, PParT, PPerT
    
    implicit none
    
    integer                        :: s, nrIn, ntIn, neIn, npaIn
    character(len=2), dimension(4) :: NameSpecies = (/ 'e_','h_','he','o_' /)
    character(len=100)             :: NameFile, StringLine

    character(len=*), parameter :: NameSub='read_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    call write_prefix
    write(*,*) 'Loading restart files.'

    ! Open ascii info file, read start time, elapsed time, and grid info.
    NameFile=PathRestartIn//'/restart_info.txt'
    open(unit=UnitTMP_, file=trim(NameFile), status='old')
    read(UnitTMP_,*)StringLine
    read(UnitTMP_, '(a25,i4.4, 2i2.2, 1x, 3i2.2)')StringLine, &
         TimeRamStart%iYear, TimeRamStart%iMonth, TimeRamStart%iDay, &
         TimeRamStart%iHour, TimeRamStart%iMinute, TimeRamStart%iSecond
    TimeRamStart%FracSecond=0.0
    read(UnitTMP_,'(a25, f15.4)') StringLine, TimeRestart
    read(UnitTMP_,'(a25, i15)') StringLine, nIter
    read(UnitTMP_, *) StringLine
    read(UnitTMP_, '(a25, 4i3)') StringLine, nrIn, ntIn, neIn, npaIn
    close(UnitTMP_)

    ! Debug:
    if(DoTest) then
       write(*,*)NameSub//': Restart info loaded.'
       write(*, *) 'TIMING:'
       write(*, '(a, i4.4, 2i2.2, 1x, 3i2.2)')'Start (YYYYMMDD HHMMSS)= ', &
            TimeRamStart%iYear, TimeRamStart%iMonth, TimeRamStart%iDay, &
            TimeRamStart%iHour, TimeRamStart%iMinute,TimeRamStart%iSecond
       write(*, '(a,f13.7)')'Elapsed (seconds)      = ', TimeRestart
       write(*, '(a, i15)') 'Completed Iterations   = ', nIter
       write(*, *)'GRID:'
       write(*, '(a, 4i3)') 'nR, nL, nE, nPA        = ', nrIn,ntIn,neIn,npaIn
    end if

    ! Update TimeRamNow
    TimeRamElapsed=TimeRestart
    call time_int_to_real(TimeRamStart)
    TimeRamNow % Time = TimeRamStart % Time + TimeRamElapsed
    call time_real_to_int(TimeRamNow)

    if(DoTest)then
       write(*,*)'Restart Timing:'
       write(*,*)'Start Time   =', TimeRamStart%String
       write(*,'(a,f13.2)')' Time Elapsed =',TimeRamElapsed
       write(*,*)'Current Time =', TimeRamNow%String
    end if

    ! Verify grid.
    if( (nR/=nrIn) .or. (nT/=ntIn) .or. (nE/=neIn) .or. (nPA/=npaIn)) &
         call CON_stop(NameSub//': Restart grid does not match current grid.')

    ! Load restart data.
    do s=1, 4
       NameFile=PathRestartIn//'/restart_'//NameSpecies(s)//'.rst'
       if(DoTest) write(*,'(4a)') 'RAM: Reading restart file for ', &
            NameSpecies(s), ' = ', NameFile
       open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
       read(UnitTMP_) F2(s,:,:,:,:)
       close(UnitTMP_)
    end do

    NameFile=PathRestartIn//'/restart_ppar.rst'
    if(DoTest) write(*,'(3a)') 'RAM: Reading restart file parallel pressure ', &
         ' = ', NameFile
    open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
    read(UnitTMP_) PParT(:,:,:)
    close(UnitTMP_)

    NameFile=PathRestartIn//'/restart_pper.rst'
    if(DoTest) write(*,'(3a)') 'RAM: Reading restart file perpendicular pressure ', &
         ' = ', NameFile
    open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
    read(UnitTMP_) PPerT(:,:,:)
    close(UnitTMP_)

  end subroutine read_restart
  !==========================================================================

end module ModRamIO
!============================================================================
