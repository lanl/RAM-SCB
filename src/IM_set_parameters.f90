!==============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================
subroutine IM_set_parameters

  !use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModRamMpi
  use ModRamSats, ONLY: read_sat_params
  use Module1,    ONLY: decreaseConvAlphaMin, decreaseConvAlphaMax,&
       decreaseConvPsiMin, decreaseConvPsiMax, &
       blendAlphaInit, blendPsiInit, iDumpRAMFlux
  use ModRamIO,   ONLY: DtLogfile, UseNewFmt
  use ModRamMain, ONLY: event, boundary, NameBoundMag, electric, UseEfind, &
                      TimeMax, IsComponent, MaxIter, TimeRamStart, &
                      StringTest, StrRamDescription, DoUseScb, DoWriteFlux, &
                      DoMultiBcsFile, DtsMax, BetaLim, DoVarDt, electrons, &
                      IsStarttimeSet, IsRestart, DtRestart, IsSHIELDS, &
                      DoSaveFinalRestart, NameDistrib, DoUsePlane_SCB, &
                      DoUseWPI, DoWriteB, DoWriteCur, DoWritePres, &
                      DoUseBASdiff, DoUseKpDiff, DoUseVAPini

  use ModRamIndices, ONLY: NameOmniFile
  use ModTimeConvert

  implicit none

  integer :: iDate
  logical :: TempLogical
  character (len=100)           :: NameCommand
  character (len=*), parameter  :: NameSub = 'IM_set_parameters'

  !---------------------------------------------------------------------------
 
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!! RAM-SCB specific Params:
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case('#EVENT')
        call read_var('NameEvent', event)

     case('#USESCB')
        call read_var('DoUseScb', DoUseScb)
        call read_var('DoWriteB', DoWriteB)
        call read_var('DoWriteCur', DoWriteCur)
!        call read_var('DoWritePres', DoWritePres)

     case('#USEPLANE')
        call read_var('DoUsePlane_SCB', DoUsePlane_SCB)

     case('#USEWPI')
        call read_var('DoUseWPI', DoUseWPI)
        call read_var('DoUseBASdiff', DoUseBASdiff)
        call read_var('DoUseKpDiff', DoUseKpDiff)

     case('#OUTERBOUNDARY')
        call read_var('NameBoundPlasma', boundary)
        call read_var('NameBoundMag', NameBoundMag)
        call read_var('NameDistribution', NameDistrib)

     case('#MULTISPECIESBCS')
        call read_var('DoMultispeciesBcs', DoMultiBcsFile)
        call read_var('DoElectrons', electrons)

     case('#EFIELD')
        call read_var('NameEfield', electric)
        call read_var('UseEfind', UseEfind)

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
        call read_var('BlendAlpha'       , blendAlphaInit)
        call read_var('BlendPsi'         , blendPsiInit)

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
        if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand), ' !!!'
           call CON_stop(NameSub // ' Correct PARAM.in!')
        end if

     end select
  end do

end subroutine IM_set_parameters
