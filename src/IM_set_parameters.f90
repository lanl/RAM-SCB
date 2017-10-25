!==============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================
subroutine IM_set_parameters

  use ModRamMain,    ONLY: PathRestartIn, nIter
  use ModRamGrids,   ONLY: NEL, NTL, NR, NT, NE, NPA, EnergyMin
  use ModRamTiming,  ONLY: TimeRamElapsed, TimeRamStart, TimeRamRealStart, &
                           TimeRamNow, DtLogFile, DtRestart, DtsMax, TimeMax, &
                           TimeRestart, MaxIter, Dt_hI, Dt_bc, DtEfi, DtW_hI, &
                           DtW_Pressure, DtW_EField, DtsMin, TOld
  use ModRamParams

  use ModScbMain,  ONLY: method
  use ModScbGrids, ONLY: nthe, npsi, nzeta
  use ModScbParams

  use ModRamSats, ONLY: read_sat_params  

!--- SCE Components
  use ModRamCouple, ONLY: DoPassJr, DoIEPrecip
  use ModIE_Interface
  use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case
  use ModIonosphere
  use IE_ModIo
  use IE_ModMain
!---

  use ModReadParam
  use ModIOUnit, ONLY: UNITTMP_
  use ModTimeConvert, ONLY: time_real_to_int, time_int_to_real
  implicit none

  integer :: iPressure, iDomain, iFile, iDebugProc
  character(len=4) :: sPressure, sDomain

  integer :: iDate, nrIn, ntIn, neIn, npaIn, dummyi
  logical :: TempLogical, UseStrict = .true.
  character(len=50)  :: plot_string
  character(len=100) :: StringLine, NameCommand, RestartFile
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

!     case('#CONSTRAINTS')
!        call read_var('MinEnergy', EnergyMin)

     case('#FILESOURCE')
        call read_var('UseSWMFFile',UseSWMFFile)

     case('#RAMGRID')
        call read_var('nR',  NR)
        call read_var('nT',  NT)
        call read_var('nE',  NE)
        call read_var('nPa', NPA)

     case('#SCBGRID')
        call read_var('nTheta', nthe)
        call read_var('nPsi',   npsi)
        call read_var('nZeta',  nzeta)

     case('#INDICES_FILE')
        call read_var('NameIndicesFile', NameIndexFile)

     case('#COMPONENT_TIMESTEPS')
        call read_var('SCBTimeStep', Dt_hI)
        call read_var('BCTimeStep',  Dt_bc)
        call read_var('EFTimeStep',  DtEfi)

     case('#OUTPUT_FREQUENCY')
        call read_var('DtPressureFileWrite', DtW_Pressure)
        if (DtW_Pressure.lt.1.0) DtW_Pressure = 9999999999.9
        call read_var('DthIFileWrite',       DtW_hI)
        if (DtW_hI.lt.1.0) DtW_hI = 9999999999.9
        call read_var('DtEFieldFileWrite',   DtW_EField)
        if (DtW_EField.lt.1.0) DtW_EField = 9999999999.9

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

     case('#MINTIMESTEP')
        call read_var('MinHalfStep', DtsMin)

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
        call read_var('BlendMin'            , blendMin)
        call read_var('BlendMax'            , blendMax)

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

!--- SCE Components
     !!! IE coupled
     case("#STRICT")
        call read_var('UseStrict',UseStrict)
     case("#SAVEPLOT", "#IE_SAVEPLOT")
        call read_var('nPlotFile',nFile)
        if (nFile > MaxFile)call CON_stop(NameSub//&
             ' IE_ERROR number of ouput files is too large in #IE_SAVEPLOT:'&
             //' nFile>MaxFile')
        do iFile=1,nFile

           call read_var('StringPlot',plot_string)
           call lower_case(plot_string)

           ! Plotting frequency
           call read_var('DnSavePlot',dn_output(iFile))
           call read_var('DtSavePlot',dt_output(iFile))

           ! Plot file format
           if(index(plot_string,'idl')>0)then
              plot_form(iFile)='idl'
           elseif(index(plot_string,'tec')>0)then
              plot_form(iFile)='tec'
           else
              call CON_stop(NameSub//&
                   ' IE_ERROR format (idl,tec) missing from plot_string='&
                   //plot_string)
           end if
           if(index(plot_string,'min')>0)then
              plot_vars(iFile)='minimum'
           elseif(index(plot_string,'max')>0)then
              plot_vars(iFile)='maximum'
           elseif(index(plot_string,'aur')>0)then
              plot_vars(iFile)='aur'
           elseif(index(plot_string,'uam')>0)then
              plot_vars(iFile)='uam'
           else
              call CON_stop(NameSub//&
                   ' IE_ERROR variable definition missing in #IE_SAVEPLOT'//&
                   ' from plot_string='//plot_string)
           end if
        end do
     case("#SAVEPLOTNAME")
        call read_var('IsPlotName_e',IsPlotName_e)
     case("#SAVELOGNAME")
        call read_var('IsLogName_e',IsLogName_e)
     case("#IONOSPHERE")
        call read_var('iConductanceModel',conductance_model)
        call read_var('UseFullCurrent' ,UseFullCurrent)
        call read_var('UseFakeRegion2' ,UseFakeRegion2)
        call read_var('F10.7 Flux',f107_flux)
        call read_var('StarLightPedConductance',StarLightPedConductance)
        call read_var('PolarCapPedConductance',PolarCapPedConductance)
     case("#IM")
        call read_var('TypeImCouple',TypeImCouple)
        call lower_case(TypeImCouple)
        call read_var('FractionImJr',FractionImJr)
     case("#BOUNDARY")
        call read_var('LatBoundary',LatBoundary)
        LatBoundary = LatBoundary * cDegToRad
     case("#UA")
        call read_var('DoCoupleUaCurrent',DoCoupleUaCurrent)
        if(DoCoupleUaCurrent)then
           call read_var('LatBoundary',LatBoundary)
           LatBoundary = LatBoundary * cDegToRad
        endif
     case("#SPS")
        call read_var('UseSPS',UseSPS)
        IE_NameOfEFieldModel = "SPS"
        UseGridBasedIE = .true.
     case("#SOLVER")
        call read_var('NameSolver',        NameSolver, IsLowerCase=.true.)
     case("#KRYLOV")
        call read_var('UsePreconditioner', UsePreconditioner)
        call read_var('UseInitialGuess',   UseInitialGuess)
        call read_var('Tolerance',         Tolerance)
        call read_var('MaxIteration',      MaxIteration)
     case("#DEBUG")
        call read_var('iDebugLevel',iDebugLevel)
        call read_var('iDebugProc',iDebugProc)


     case("#BACKGROUND")

        call read_var('NameOfModelDir',IE_NameOfModelDir)
        call read_var('NameOfEFieldModel',IE_NameOfEFieldModel)
        call read_var('NameOfAuroralModel',IE_NameOfAuroralModel)
        call read_var('NameOfSolarModel',IE_NameOfSolarModel)

        if (index(IE_NameOfAuroralModel,'IHP') > 0) &
             IE_NameOfAuroralModel = 'ihp'
        if (index(IE_NameOfAuroralModel,'PEM') > 0) &
             IE_NameOfAuroralModel = 'pem'

        if (index(IE_NameOfEFieldModel,'AMIE') > 0) &
             IE_NameOfEFieldModel = 'amie'

        if (index(IE_NameOfEFieldModel,'weimer01') > 0) &
             IE_NameOfEFieldModel = 'weimer01'
        if (index(IE_NameOfEFieldModel,'Weimer01') > 0) &
             IE_NameOfEFieldModel = 'weimer01'
        if (index(IE_NameOfEFieldModel,'WEIMER01') > 0) &
             IE_NameOfEFieldModel = 'weimer01'

        if (index(IE_NameOfEFieldModel,'weimer') > 0 .and. &
             index(IE_NameOfEFieldModel,'01') == 0) &
             IE_NameOfEFieldModel = 'weimer96'
        if (index(IE_NameOfEFieldModel,'Weimer') > 0 .and. &
             index(IE_NameOfEFieldModel,'01') == 0) &
             IE_NameOfEFieldModel = 'weimer96'
        if (index(IE_NameOfEFieldModel,'WEIMER') > 0 .and. &
             index(IE_NameOfEFieldModel,'01') == 0) &
             IE_NameOfEFieldModel = 'weimer96'

        if (index(IE_NameOfEFieldModel,'weimer96') > 0) &
             IE_NameOfEFieldModel = 'weimer96'
        if (index(IE_NameOfEFieldModel,'Weimer96') > 0) &
             IE_NameOfEFieldModel = 'weimer96'
        if (index(IE_NameOfEFieldModel,'WEIMER96') > 0) &
             IE_NameOfEFieldModel = 'weimer96'

        if (index(IE_NameOfEFieldModel,'SAMIE') > 0) &
             IE_NameOfEFieldModel = 'samie'

        UseGridBasedIE = .false.

     case("#SAVELOGFILE")
        call read_var('DoSaveLogfile',DoSaveLogfile)

     case("#CONDUCTANCE")
        call read_var('OvalWidthFactor',OvalWidthFactor)
        call read_var('OvalStrengthFactor',OvalStrengthFactor)
        if (conductance_model/=4) then
           write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                ' command '//trim(NameCommand)// &
                ' can only be used with conductance model 4'
           if(UseStrict)call CON_stop('Correct PARAM.in!')
        end if
!---

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
  case('PTM')
     iPressure = 6
     sPressure = 'LANL'
     NEL = 36
     NTL = 25
  case('QDM')
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
     RestartFile=PathRestartIn//'/restart_info.txt'
     open(unit=UnitTMP_, file=trim(RestartFile), status='old')
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
  TOld = TimeRamElapsed

end subroutine IM_set_parameters
