!==============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================
subroutine IM_set_parameters

  use ModRamMain, ONLY: PathRestartIn, nIter
  use ModRamGrids, ONLY: NEL, NTL
  use ModRamTiming, ONLY: TimeRamElapsed, TimeRamStart, TimeRamRealStart, TimeRamNow, &
                          DtLogFile, DtRestart, DtsMax, TimeMax, TimeRestart, MaxIter
  use ModRamIndices, ONLY: NameOmniFile
  use ModRamParams

  use ModScbMain, ONLY: method
  use ModScbParams

  use ModRamSats, ONLY: read_sat_params  

  use ModReadParam
  use ModIOUnit, ONLY: UNITTMP_
  use ModTimeConvert, ONLY: time_real_to_int, time_int_to_real
  implicit none

  integer :: iPressure, iDomain
  character(len=4) :: sPressure, sDomain

  integer :: iDate, nrIn, ntIn, neIn, npaIn
  logical :: TempLogical
  character(len=100) :: StringLine, NameCommand, NameFile
  character(len=*), parameter  :: NameSub = 'IM_set_parameters'

  !---------------------------------------------------------------------------
 
  do
     if (.not.read_line()) EXIT
     if (.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! RAM-SCB specific Params:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case('#EVENT')
        call read_var('NameEvent', event)

     case('#USESCB')
        call read_var('DoUseScb',   DoUseScb)
        call read_var('DoWriteB',   DoWriteB)
        call read_var('DoWriteCur', DoWriteCur)
!        call read_var('DoWritePres', DoWritePres)

     case('#USEPLANE')
        call read_var('DoUsePlane_SCB', DoUsePlane_SCB)

     case('#USEWPI')
        call read_var('DoUseWPI',     DoUseWPI)
        call read_var('DoUseBASdiff', DoUseBASdiff)
        call read_var('DoUseKpDiff',  DoUseKpDiff)

     case('#OUTERBOUNDARY')
        call read_var('NameBoundPlasma',  boundary)
        call read_var('NameBoundMag',     NameBoundMag)
        call read_var('NameDistribution', NameDistrib)

     case('#MULTISPECIESBCS')
        call read_var('DoMultispeciesBcs', DoMultiBcsFile)
        call read_var('DoElectrons',       electrons)

     case('#EFIELD')
        call read_var('NameEfield', electric)
        call read_var('UseEfind',   UseEfind)

     case('#SATELLITE')
        call read_sat_params()
        call read_var('DoUseVAPini',DoUseVAPini)

     case('#LOGFILE')
        call read_var('DtWrite', DtLogfile)
        
     case('#NAMEFORMAT')
        call read_var('UseNewFormat', UseNewFmt)

     case('#SAVEFLUX')
        call read_var('DoSaveFlux',DoWriteFlux)

     case('#DUMP3DFLUX')
        call read_var('DoDump3dFlux', TempLogical)
        if(TempLogical) then
           iDumpRAMFlux=1 
        else 
           iDumpRAMFlux=0
        end if

     case('#OMNIFILE')
        call read_var('NameOmniFile', NameOmniFile)

     case('#MAXTIMESTEP')
        call read_var('MaxHalfStep', DtsMax)

     case('#VARIABLEDT')
        call read_var('DoVariableDt', DoVarDt)

     case('#RAMLIMITER')
        call read_var('LimiterBeta', BetaLim)
        if(BetaLim > 2.0) BetaLim=2.0
        if(BetaLim < 1.0) BetaLim=1.0

     case('#SCBSCHEME')
        call read_var('DecreaseConvAlphaMin', decreaseConvAlphaMin)
        call read_var('DecreaseConvAlphaMax', decreaseConvAlphaMax)
        call read_var('DecreaseConvPsiMin'  , decreaseConvPsiMin)
        call read_var('DecreaseConvPsiMax'  , decreaseConvPsiMax)
        call read_var('BlendAlpha'          , blendAlphaInit)
        call read_var('BlendPsi'            , blendPsiInit)

     ! Restart commands:
     case('#RESTART')
        IsRestart=.true.
     case('#SAVERESTART')
        call read_var('DtSaveRestart', DtRestart)
        call read_var('DoSaveFinalRestart', DoSaveFinalRestart)

     ! Activate "SHIELDS-RC" mode?
     case('#SHIELDSRC')
        IsSHIELDS=.true.  ! "SHIELDS-RC" mode (simpler runs)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Shared w/ Framework:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case('#STOP')
        if(IsComponent)then
           write(*,*)'IM WARNING: #STOP command is ignored in the framework'
        else
           call read_var('MaxIter', MaxIter)
           call read_var('TimeMax', TimeMax)
        end if

     case('#STARTTIME')
        if(IsComponent)then
           write(*,*) &
                'IM WARNING: #STARTTIME command is ignored in the framework'
        else
           IsStarttimeSet=.true.
           !read in iYear,iMonth,iDay,iHour,iMinute,iSecond,FracSecond
           call read_var('iYear'  ,   TimeRamStart % iYear)
           call read_var('iMonth' ,   TimeRamStart % iMonth)
           call read_var('iDay'   ,   TimeRamStart % iDay)
           call read_var('iHour'  ,   TimeRamStart % iHour)
           call read_var('iMinute',   TimeRamStart % iMinute)
           call read_var('iSecond',   TimeRamStart % iSecond)
           call read_var('FracSecond',TimeRamStart % FracSecond)
           ! Set integer time.
           call time_int_to_real(TimeRamStart)

        end if

     case('#DESCRIPTION')
        if(IsComponent)then
           write(*,*) &
                'IM WARNING: #DESCRIPTION should not be set ' , &
                'by IM in component mode.'
        else
           call read_var('StringDescription',StrRamDescription)
        end if

     case('#TEST')
        if(IsComponent)then
           write(*,*) &
                'IM WARNING: #TEST should not be set by IM in component mode.'
        else
           call read_var('StringTest',StringTest)
        endif

     !!! Default crash:
     case default
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand), ' !!!'
           call CON_stop(NameSub // ' Correct PARAM.in!')

     end select
  end do

  ! Check for SHIELDSRC mode.  If True, set SHIELDSRC params.
  if (IsSHIELDS)then
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
     sPressure = 'SWMF'
  case('LANL')
     iPressure = 6
     sPressure = 'LANL'
     NEL = 36
     NTL = 25
  case default
     call CON_stop(NameSub//': invalid boundary='//boundary)
  end select

  select case(NameBoundMag)
  case('SWMF')
     iDomain = 20
     sDomain = 'SWMF'
  case('DIPL') ! Dipole w/o SCB calculation.
     iDomain = 4
     method = 3
     sDomain = 'DIPL'
  case('DIPS') ! Dipole w/  SCB calculation.
     iDomain = 4
     sDomain = 'DIPS'
  case('T89C')
     iDomain = 3
     sDomain = 'T89C'
  case default
     call CON_stop(NameSub//': invalid NameBoundMag='//NameBoundMag)
  end select

  if (IsRestart) then
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
     call time_int_to_real(TimeRamStart)
     TimeRamRealStart%Time = TimeRamStart%Time + TimeRestart
     TimeRamElapsed = TimeRestart
     call time_real_to_int(TimeRamRealStart)
  else
     TimeRamElapsed = 0
     TimeRamRealStart = TimeRamStart
  end if
     TimeRamNow = TimeRamRealStart

end subroutine IM_set_parameters
