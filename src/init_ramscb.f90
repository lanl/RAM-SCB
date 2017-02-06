!******************************************************************************
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************
subroutine init_ramscb

  use ModRamMain, ONLY: iPressure, iDomain, boundary, NameBoundMag, &
       PPerH, PParH, PPerO, PParO, PPerHe, PParHe, DoInitOnly, PPerE, PParE, &
       IsComponent,iCal, TimeRamStart, TimeRamNow, RadiusMin, &
       RadiusMax, nR, nT, Lz, Phi, NTL, NEL, NEL_prot, NBD, IsStarttimeSet, &
       TimeRamElapsed, timeMax, IsRestart, electrons, kp, f107, DoVarDt, &
       electric, UseEfind, DoWriteFlux, IsSHIELDS, nRextend, GridExtend
  use Module1,  ONLY: method, iDumpRAMFlux
  use ModRamIO, ONLY: StringRunDate, read_restart, UseNewFmt, write_prefix, &
       DtLogFile
  use ModNumConst,ONLY: cTwoPi
  use ModTimeConvert
  use ModRamIndices, ONLY: init_indices, get_indices
  use CON_planet, ONLY: set_planet_defaults
  use CON_axes,   ONLY: init_axes, test_axes
  use ModPlanetConst, ONLY: init_planet_const
  use ModRamTiming, ONLY: init_timing

  use ModRamMpi

  implicit none

  character(len=8) :: StringSysDate
  character(len=10):: StringSysTime

  type(timetype) :: TimeStop

  real:: dR, dPhi
  integer:: iR, iPhi, iS

  character (len=*), parameter:: NameSub='init_ramscb'
  logical :: DoTest, DoTestMe
  !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Build time and date of when this run was started.
  call date_and_time(StringSysDate, StringSysTime)
  StringRunDate = StringSysDate // ' - ' // StringSysTime
  if(DoTest)write(*,*) 'Run started on ' // StringRunDate

  ! Check for SHIELDSRC mode.  If True, set SHIELDSRC params.
  if(IsSHIELDS)then
     call write_prefix
     write(*,*)' RAM-SCB is in SHIELDS-RC MODE.'
     DoVarDt      = .true.  ! Variable timestepping.
     boundary     = 'LANL'  ! LANL plasma BCS.
     NameBoundMag = 'DIPL'  ! Simple dipole BCS.
     electric     = 'VOLS'  ! Volland-Stern E-field.
     UseEfind     = .false. ! No induced E-field.
     UseNewFmt    = .true.  ! New output file name format.
     DoWriteFlux  = .false. ! No large flux files.
     iDumpRAMFlux = 0       ! No large flux files.
     electrons    = .true.  ! Electrons ON.
     DtLogFile    = 300.0   ! Sparser log files.
  end if

  ! Set iPressure and iDomain according to selected BC's
  select case(boundary)
  case('SWMF')
     iPressure = 5
  case('LANL')
     iPressure = 6
     NEL = 36
     NTL = 25
  case('PTM')
     iPressure = 6
     NEL = 9 !These are Joachim's energy levels
     NEL_prot = 35 ! These are the SWMF energy levels
     NTL = 25
     NBD = 130
  case('TM03')
     iPressure = 5 ! Both RAM and SCB model expanded to the same distance
  case default
     call CON_stop(NameSub//': invalid boundary='//boundary)
  end select

  select case(NameBoundMag)
  case('SWMF')
     iDomain = 20
  case('T89C')
     iDomain = 3
  case('TS04')
     iDomain = 2
  case('DIPL') ! Dipole w/o SCB calculation.
     iDomain = 4
     method = 3
  case('DIPS') ! Dipole w/  SCB calculation.
     iDomain = 4
     ! method  = 3 ! DIPS is the case with calculation (within dipole boundary)
  case default
     call CON_stop(NameSub//': invalid NameBoundMag='//NameBoundMag)
  end select

  if(IsRestart) then
     ! If Restart, read restart params and set timings appropriately.
     if (IsStarttimeSet) call CON_stop(NameSub//&
          ': Cannot use #STARTTIME command with #RESTART!')
     call read_restart
  else
     ! Otherwise, initialize time correctly.
     TimeRamNow % Time = TimeRamStart % Time
     call time_real_to_int(TimeRamNow)
     TimeRamElapsed=0.0
  end if

  iCal = 1   !AND set iCal properly.

    ! Begin tracking code efficiency.
  call init_timing()

  ! Initialize Pressures.
  PPerH  = 0
  PParH  = 0
  PPerO  = 0
  PParO  = 0
  PPerHe = 0
  PParHe = 0
  PPerE  = 0
  PParE  = 0

  ! Initialize grid.
  ! Radial distance
  dR = (RadiusMax - RadiusMin)/(nR - 1)
  do iR = 1, nR+1 ! DANGER WE SHOULD CHANGE THIS AND ALL ARRAYS USING NR+1
     Lz(iR) = RadiusMin + (iR - 1)*dR
  end do

  ! Create extended radial grid for coupling:
  do iR=1, nRextend
     GridExtend(iR) = RadiusMin + (iR-1)*dR
  end do

  ! Longitude in radians
  dPhi = cTwoPi/(nT - 1)
  do iPhi = 1, nT
     Phi(iPhi) = (iPhi - 1)*dPhi
  end do

  ! Call ram_all.f in initialization mode to 
  ! generate pitch angle and energy grid.
  DoInitOnly = .true.
  if (electrons) call ram_all(1)
  do iS=2,4
     call ram_all(iS)
  end do
  DoInitOnly = .false.

  ! Set up input indices (Kp, F10.7)
  TimeStop%Time = TimeRamStart%Time + TimeMax
  call time_real_to_int(TimeStop)
  call init_indices(TimeRamStart, TimeStop)
  call get_indices(TimeRamNow%Time, Kp, f107)

  ! Initialize planet and axes if in stand-alone mode.
  ! In component mode, the framework takes care of this.
  ! This is needed for correct coordinate transformations.
  if(.not. IsComponent) then
     call init_planet_const
     call set_planet_defaults
     call init_axes(TimeRamStart % Time)
  end if

  ! Initialize output files as necessary.
  if(iProc==0) call init_output()

contains

!============================================================================
  subroutine init_output
    ! Initialize any output that requires such action.

    use ModRamSats, ONLY: read_sat_input, init_sats
    use ModRamMain, ONLY: PathRamOut, PathScbOut, nIter
    use ModRamIO,   ONLY: iUnitLog, DoSaveRamSats, NameFileLog, NameFileDst
    use ModIoUnit,  ONLY: io_unit_new
    use Module1, ONLY : iUnitDst

    character(len=*), parameter :: NameSubSub = NameSub//'::init_output'
    !------------------------------------------------------------------------
    ! Initialize virtual satellites.
    if(DoSaveRamSats) then
       call read_sat_input
       call init_sats
    end if

    ! Initialize logfile.
    iUnitLog = io_unit_new()
    write(NameFileLog, '(a,i6.6, a)') PathRamOut//'log_n', nIter, '.log'
    open(iUnitLog, FILE=NameFileLog, STATUS='REPLACE')
    write(iUnitLog,*)'RAM-SCB Log'
    write(iUnitLog, '(a,a)') &
         'time year mo dy hr mn sc msc dstRam dstBiot pparh pperh ', & 
         'pparo ppero pparhe pperhe ppare ppere'
    
    ! Initialize Dstfile.
    iUnitDst = io_unit_new()
    write(NameFileDst, '(a,i6.6, a)') PathScbOut//'dst_computed_3deq_', nIter, '.txt'
    open(iUnitDst, FILE=NameFileDst, STATUS='REPLACE')
    write(iUnitDst,*)'SCB Dst File'
    write(iUnitDst, *) &
         'time year mo dy hr mn sc msc dstDPS dstDPSGeo ', & 
         'dstBiot dstBiotGeo'
       
  end subroutine init_output

end subroutine init_ramscb
!============================================================================
