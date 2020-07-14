!==============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================
subroutine IM_set_parameters

!!!!! Module Variables
  use ModRamMain,    ONLY: PathRestartIn, nIter, DP
  use ModRamGrids,   ONLY: nS, NEL, NTL, NR, NT, NE, NPA, NameVar
  use ModRamVariables, ONLY: composition
  use ModRamTiming,  ONLY: TimeRamElapsed, TimeRamStart, TimeRamRealStart, &
                           TimeRamNow, DtLogFile, DtRestart, DtsMax, TimeMax, &
                           TimeRestart, MaxIter, Dt_hI, Dt_bc, DtEfi, DtW_hI, &
                           DtW_Pressure, DtW_EField, DtsMin, TOld, TimeRamFinish, &
                           DtW_MAGxyz, DtW_2DFlux
  use ModScbGrids, ONLY: nthe, npsi, nzeta
  use ModRamParams
  use ModScbParams

!!!! Module Subroutines and Functions
  use ModRamIO,   ONLY: write_prefix
  use ModRamSats, ONLY: read_sat_params  

!!!! External Modules (share/Library)
  use ModReadParam
  use ModIOUnit, ONLY: UNITTMP_
  use ModTimeConvert, ONLY: time_real_to_int, time_int_to_real


  implicit none

  integer :: i, nChar, nrIn, ntIn, neIn, npaIn
  logical :: TempLogical
  logical :: StopCommand, IsStopTimeSet
  real(DP) :: TempReal
  character(len=100) :: StringLine, NameCommand, RestartFile
  character(len=*), parameter  :: NameSub = 'IM_set_parameters'
  StopCommand = .false.
  IsStopTimeSet = .false.
  !---------------------------------------------------------------------------
 
  do
     if (.not.read_line()) EXIT
     if (.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! RAM-SCB specific Params:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Misc Parameters
     case('#EVENT')
        call read_var('NameEvent', event)

     case('#VERBOSE')
        verbose = .true.

     !case('#TEST') - Shared with the framework
     !case('#DESCRIPTION') - Shared with the framework

!!!!!! Timing Parameters
     !case('#STARTTIME') - Shared with the framework
     !case('#STOPTIME') - Shared with the framework
     !case('#STOP') - Shared with the framework

     case('#COMPONENT_TIMESTEPS')
        call read_var('SCBTimeStep', Dt_hI)
        call read_var('BCTimeStep',  Dt_bc)
        call read_var('EFTimeStep',  DtEfi)

     case('#TIE_SCB_TO_RAM')
        SCBonRAMTime = .true.

     case('#MAXTIMESTEP')
        call read_var('MaxHalfStep', DtsMax)

     case('#MINTIMESTEP')
        call read_var('MinHalfStep', DtsMin)

     case('#VARIABLEDT')
        call read_var('DoVariableDt', DoVarDt)

!!!!!! RAM Parameters
     case('#USERAM')
        call read_var('DoUseRAM', DoUseRAM)

     case('#PLASMASPHERE')
        call read_var('UsePlasmasphere', DoUsePlasmasphere)
        if (DoUsePlasmasphere) then
           call read_var('PlasmasphereModel', PlasmasphereModel)
           call read_var('TauCalculation', TauCalculation)
        endif

     case('#COULOMB')
        call read_var('UseCoulomb', DoUseCoulomb)

     case('#SPECIES')
        call read_var('nS', nS)
        call read_var('NameVar', NameVar)
        call read_var('FixedComp',FixedComposition)
        if (FixedComposition) then
           allocate(composition(nS))
           do i=1,nS
              call read_var('Composition',TempReal)
              composition(i) = TempReal/100
           enddo
        endif

     case('#NITROGEN_PERCENT')
        call read_var('NitrogenPercent', TempReal)
        OfracN = TempReal/100

     case('#FLAT_INITIALIZATION')
        InitializeOnFile = .false.

     case('#CHECK_MAGNETOPAUSE')
        checkMGNP = .true.

     case('#USEWPI')
        call read_var('DoUseWPI',     TempLogical)
        if (TempLogical) DoUseWPI = .true.
        call read_var('DoUseBASdiff', TempLogical)
        if (TempLogical) DoUseBASDiff = .true.
        call read_var('DoUseKpDiff',  TempLogical)
        if (TempLogical) DoUseKpDiff = .true.
        call read_var('DoUseEMIC',  TempLogical)
        if (TempLogical) DoUseEMIC = .true.

     case('#USEFLC')
        call read_var('DoUseFLC',DoUseFLC)
        call read_var('DoWriteFLCDiffCoeff', DoWriteFLCDiffCoeff)
        
     case('#FLUX_CAP')
        call read_var('ElectronFluxCap', ElectronFluxCap)
        call read_var('ProtonFluxCap', ProtonFluxCap)

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

     case('#RAMLIMITER')
        call read_var('LimiterBeta', BetaLim)
        if(BetaLim > 2.0) BetaLim=2.0
        if(BetaLim < 1.0) BetaLim=1.0

     case('#RAMGRID')
        call read_var('nR',  NR)
        call read_var('nT',  NT)
        call read_var('nE',  NE)
        call read_var('nPa', NPA)

     case('#SMOOTH_INTEGRAL')
        integral_smooth = .false.

!!!!!! SCB Parameters
     case('#RESET')
        reset=.true.

     case('#PRESS_MODE')
        call read_var('PressMod', PressMode)

     case('#SCB_FIELD')
        call read_var('constZ',constZ)
        call read_var('constTheta', constTheta)

     case('#SCB_CONVERGENCE')
        call read_var('ConvergenceDistance',convergence)

     case('#SCBSMOOTHING')
        call read_var('PressureSmoothing', iSm2)
        call read_var('SavitzkyGolayIterations', SavGolIters)

     case('#SCB_GHOSTCELLS')
        call read_var('FixedOuterShell', psiChange)
        call read_var('FixedFootPoints', theChange)

     case('#SCBSETTINGS')
        call read_var('MinSCBIterations', MinSCBIterations)
        call read_var('SCBMethod', method)
        call read_var('OuterMethod', iOuterMethod)
        call read_var('InnerConvergenceAlpha', InConAlpha)
        call read_var('InnerConvergencePsi',   InConPsi)

     case('#SCBFLAGS')
        call read_var('Isotropic', TempLogical)
           if (TempLogical) then
              isotropy = 1
           else
              isotropy = 0
           endif
        call read_var('ReduceAnisotropy', TempLogical)
           if (TempLogical) then
              iReduceAnisotropy = 1
           else
              iReduceAnisotropy = 0
           endif
        call read_var('BetaExtrapolation', TempLogical)
           if (TempLogical) then
              iWantAlphaExtrapolation = 1
           else
              iWantAlphaExtrapolation = 0
           endif
        call read_var('AzimuthalOffset', TempLogical)
           if (TempLogical) then
              iAzimOffset = 2
           else
              iAzimOffset = 1
           endif
        call read_var('EmptyLossCone', TempLogical)
           if (TempLogical) then
              iLossCone = 2
           else
              iLossCone = 1
           endif
        call read_var('AdaptiveMesh', TempLogical)
           if (TempLogical) then
              iAMR = 1
           else
              iAMR = 0
           endif

     case('#SCBDETAILS')
        call read_var('SORDetail', TempLogical)
           if (TempLogical) then
              isSORDetailNeeded = 1
           else
              isSORDetailNeeded = 0
           endif
        call read_var('EnergyDetail', TempLogical)
           if (TempLogical) then
              isEnergDetailNeeded = 1
           else
              isEnergDetailNeeded = 0
           endif
        call read_var('ForceBalanceDetail', TempLogical)
           if (TempLogical) then
              isFBDetailNeeded = 1
           else
              isFBDetailNeeded = 0
           endif
        call read_var('PressureDetail', TempLogical)
           if (TempLogical) then
              isPressureDetailNeeded = 1
           else
              isPressureDetailNeeded = 0
           endif

     case('#SCBGRID')
        call read_var('nTheta', nthe)
        call read_var('nPsi',   npsi)
        call read_var('nZeta',  nzeta)

     case('#SCBSCHEME')
        call read_var('DecreaseConvAlphaMin', decreaseConvAlphaMin)
        call read_var('DecreaseConvAlphaMax', decreaseConvAlphaMax)
        call read_var('DecreaseConvPsiMin'  , decreaseConvPsiMin)
        call read_var('DecreaseConvPsiMax'  , decreaseConvPsiMax)
        call read_var('BlendAlpha'          , blendAlphaInit)
        call read_var('BlendPsi'            , blendPsiInit)
        call read_var('BlendMin'            , blendMin)
        call read_var('BlendMax'            , blendMax)

!!!!!! Input Parameters
     case('#TS07_DIRECTORY')
        call read_var('TS07Directory', TS07Path)

     case('#QINDENTON_FILE_PATH')
        call read_var('QinDentonPath', QinDentonPath)

     case('#BOUNDARY_FILE_PATH')
        call read_var('BoundaryPath', BoundaryPath)

     case('#INITIALIZATION_FILE_PATH')
        call read_var('InitializationPath', InitializationPath)

     case('#INDICES_FILE')
        call read_var('NameIndicesFile', NameIndexFile)

     case('#OMNIFILESOURCE')
        call read_var('UseSWMFFile',UseSWMFFile)

     case('#OMNIFILE')
        call read_var('NameOmniFile', NameOmniFile)

!!!!!! Output Parameters
     case('#OUTPUT_FREQUENCY')
        call read_var('DtPressureFileWrite', DtW_Pressure)
        if (DtW_Pressure.lt.1.0) DtW_Pressure = 9999999999.9
        call read_var('DthIFileWrite',       DtW_hI)
        if (DtW_hI.lt.1.0) DtW_hI = 9999999999.9
        call read_var('DtEFieldFileWrite',   DtW_EField)
        if (DtW_EField.lt.1.0) DtW_EField = 9999999999.9
        call read_var('DtMAGxyzWrite', DtW_MAGxyz)
        if (DtW_MAGxyz.lt.1.0) DtW_MAGxyz = 9999999999.9

     case('#WRITE_BOUNDARY')
        WriteBoundary = .true.

     case('#WRITE_POTENTIAL')
        WritePotential = .true.

     case('#SATELLITE')
        call read_sat_params()
        call read_var('DoUseVAPini',DoUseVAPini)

     case('#LOGFILE')
        call read_var('DtWrite', DtLogfile)

     case('#NAMEFORMAT')
        call read_var('UseNewFormat', UseNewFmt)

     case('#SAVEFLUX')
        call read_var('DoSaveFlux',DtW_2DFlux)
        if (DtW_2DFlux.lt.1.0) DtW_2DFlux = 99999999999.9

     case('#DUMP3DFLUX')
        call read_var('DoDump3dFlux', TempLogical)
        if(TempLogical) then
           iDumpRAMFlux=1
        else
           iDumpRAMFlux=0
        end if

!!!!!! Restart Parameters
     case('#RESTART')
        IsRestart=.true.
     case('#SAVERESTART')
        call read_var('DtSaveRestart', DtRestart)
        call read_var('DoSaveFinalRestart', DoSaveFinalRestart)
     case('#TIMEDRESTART')
        call read_var('TimedRestartFiles', TimedRestarts)
     case('#HARDRESTART')
        IsRestart=.true.
        HardRestart = .true.

!!!!!! Activate "SHIELDS-RC" mode?
     case('#SHIELDSRC')
        IsSHIELDS=.true.  ! "SHIELDS-RC" mode (simpler runs)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! Shared w/ Framework:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case('#STOP')
        if(IsComponent)then
           write(*,*)'IM WARNING: #STOP command is ignored in the framework'
        else
           StopCommand = .true.
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

     case('#STOPTIME')
        if(IsComponent)then
           write(*,*) &
                'IM WARNING: #STOPTIME command is ignored in the framework'
        else
           IsStopTimeSet=.true.
           !read in iYear,iMonth,iDay,iHour,iMinute,iSecond,FracSecond
           call read_var('iYear'  ,   TimeRamFinish % iYear)
           call read_var('iMonth' ,   TimeRamFinish % iMonth)
           call read_var('iDay'   ,   TimeRamFinish % iDay)
           call read_var('iHour'  ,   TimeRamFinish % iHour)
           call read_var('iMinute',   TimeRamFinish % iMinute)
           call read_var('iSecond',   TimeRamFinish % iSecond)
           call read_var('FracSecond',TimeRamFinish % FracSecond)
           ! Set integer time.
           call time_int_to_real(TimeRamFinish)
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

  if (SCBonRAMTime) then
     RAMTie = floor(Dt_hI)
     Dt_hI = 9999999999.9
  endif

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
     NEL = NEL
     NTL = NTL
     BoundaryFiles = .false.
     !PressMode = 'BAT'
  case('LANL')
     NEL = 36
     NTL = 25
  case('PTM')
     NEL = 36
     NTL = 25
  case('QDKP')
     NEL = 36
     NTL = 25
  case('QDBZ')
     NEL = 36
     NTL = 25
  case default
     call CON_stop(NameSub//': invalid boundary='//boundary)
  end select

  select case(NameBoundMag)
     case('SWMF')
        DoScbCalc = .true.
        InnerMag  = 'SWMF'
        OuterMag  = 'SWMF'
     case('SWML')
        DoScbCalc = .false.
        InnerMag  = 'SWMF'
        OuterMag  = 'SWMF'
     case('DIPL')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'NONE'
     case('DIPS')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'NONE'
     case('T89I')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'T89'
     case('T89D')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T89'
     case('T89L')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T89'
        NameBoundMag = 'T89D'
     case('T96I')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'T96'
     case('T96D')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T96'
     case('T96L')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T96'
        NameBoundMag = 'T96D'
     case('T02I')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'T01S'
     case('T02D')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T01S'
     case('T02L')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'T01S'
        NameBoundMag = 'T02D'
     case('T04I')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'TS05'
     case('T04D')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'TS05'
     case('T04L')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'TS05'
        NameBoundMag = 'T04D'
     case('T07I')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'TS07'
     case('T07D')
        DoScbCalc = .true.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'TS07'
     case('T07L')
        DoScbCalc = .false.
        InnerMag  = 'DIPOLE'
        OuterMag  = 'TS07'
        NameBoundMag = 'T07D'
     case('IGRF')
        DoScbCalc = .true.
        InnerMag  = 'IGRF'
        OuterMag  = 'NONE'
     case('IGRL')
        DoScbCalc = .false.
        InnerMag  = 'IGRF'
        OuterMag  = 'NONE'
        NameBoundMag = 'IGRF'
     case default
        call CON_stop(NameSub//': invalid NameBoundMag='//NameBoundMag)
  end select
  if (DoScbCalc) then
     method = 2
  else
     method = 3
  endif

end subroutine IM_set_parameters
