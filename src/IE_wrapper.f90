! Wrapper for Ridley's ionosphere
!==============================================================================
subroutine IE_set_param(CompInfo, TypeAction)

  use ModIonosphere
  use IE_ModIo
  use IE_ModMain
  use ModRamMain, ONLY: PathSceOut
  use ModIoUnit

  implicit none

  character (len=*), parameter :: NameSub='IE_set_param'

  ! Arguments
  character, intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: iError
  character (len=10) :: NameIonoDir='ionosphere'
 
  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('READ','CHECK')
     call read_param
  case('STDOUT')
     iUnitOut=STDOUT_
     StringPrefix='IE:'
  case default
     call CON_stop(NameSub//' IE_ERROR: invalid TypeAction='//TypeAction)
  end select
  
contains

  subroutine read_param

    use ModReadParam
    use ModIE_Interface
    use ModUtilities,   ONLY: fix_dir_name, check_dir, lower_case

    ! The name of the command
    character (len=100) :: NameCommand

    ! Read parameters
    logical :: DoEcho=.false., UseStrict=.true., IsUninitialized=.true.

    ! Plot file parameters
    integer :: iFile, i, iError, iDebugProc
    character (len=50) :: plot_string

    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('CHECK')
       if(IsUninitialized)call set_defaults
       IsUninitialized=.false.

       ! We should check and correct parameters here
       write(*,*) NameSub,': CHECK iSession =',i_session_read()

       RETURN
    case('READ')
       write(*,*) NameSub,': READ iSession =',i_session_read(),&
            ' iLine=',i_line_read(),' nLine =',n_line_read()

       if(IsUninitialized)call set_defaults
       IsUninitialized=.false.
    end select

    ! Read input data from text via ModReadParam
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#STRICT")
          call read_var('UseStrict',UseStrict)
       case("#IONODIR")
          call read_var("NameIonoDir",NameIonoDir)
          call fix_dir_name(NameIonoDIr)
          call check_dir(NameIonoDir)
       case("#SAVEPLOT", "#IE_SAVEPLOT")
          call read_var('nPlotFile',nFile)
          if (nFile > MaxFile)call CON_stop(NameSub//&
               ' IE_ERROR number of ouput files is too large in #IE_SAVEPLOT:'&
               //' nFile>MaxFile')
          if (nFile>0) call check_dir(NameIonoDir)
          do iFile=1,nFile

             call read_var('StringPlot',plot_string)
             call lower_case(plot_string)

             ! Check to see if the ionosphere directory exists...
             call check_dir(NameIonoDir)

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
          if (iDebugProc >= 0) then
             iDebugLevel = -1
          endif

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
          if(DoSaveLogfile)then
             call check_dir(NameIonoDir)
          endif

       case("#CONDUCTANCE")
          call read_var('OvalWidthFactor',OvalWidthFactor)
          call read_var('OvalStrengthFactor',OvalStrengthFactor)
          if (conductance_model/=4) then
             write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                  ' command '//trim(NameCommand)// &
                  ' can only be used with conductance model 4'
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if

       case default
          write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
               ' invalid command '//trim(NameCommand)
          if(UseStrict)call CON_stop('Correct PARAM.in!')
       end select
    end do

  end subroutine read_param
  !===========================================================================
  subroutine set_defaults

    conductance_model       = 5
    UseFullCurrent          = .false.
    UseFakeRegion2          = .false.
    StarLightPedConductance = 0.25
    PolarCapPedConductance  = 0.25
    f107_flux               = 150.0

  end subroutine set_defaults

end subroutine IE_set_param
!==============================================================================
subroutine IE_finalize

  use IE_ModMain,     ONLY: Time_Array, time_simulation, nSolve, DoSaveLogFile
  use IE_ModIo,       ONLY: nFile, unitlog
  use CON_physics,    ONLY: get_time
  use ModTimeConvert, ONLY: time_real_to_int
  use ModKind,        ONLY: Real8_
  implicit none

  !INPUT PARAMETERS:
!  real(Real8_),     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_finalize'

  integer :: iFile

  !---------------------------------------------------------------------------

  if(nSolve>0)then
     do iFile=1,nFile
         call ionosphere_write_output(iFile, 1)
         call ionosphere_write_output(iFile, 2)
     end do
  end if

  if(DoSaveLogfile)then
     close(unitlog)
  end if

end subroutine IE_finalize

!==============================================================================

subroutine IE_save_restart(tSimulation)

  use IE_ModMain, ONLY: nSolve
  use ModKind

  implicit none

  !INPUT PARAMETERS:
  real(Real8_),     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_save_restart'

  RETURN
  ! IE does not need a restart file at least not in the framework
  ! ionosphere_write_restart_file is still not parallel

  call ionosphere_write_restart_file(nSolve)

end subroutine IE_save_restart

!==============================================================================

subroutine IE_run(tSimulation,tSimulationLimit)

  use IE_ModMain
  use CON_physics,  ONLY: get_time, get_axes, time_real_to_int
  use ModKind
  use ModRamTiming, ONLY: TimeRamStart
  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real(real8_), intent(inout) :: tSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real(Real8_), intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  real(Real8_) :: tStart
  integer      :: nStep

  character(len=*), parameter :: NameSub='IE_run'

  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,': tSimulation,tSimulationLimit=',&
       tSimulation,tSimulationLimit

  ! Store the current time
  time_simulation = tSimulation
  
  ! Since IE is not a time dependent component, it may advance to the 
  ! next coupling time in a time accurate run
  if(time_accurate)tSimulation = tSimulationLimit

  if(DoTest)write(*,*)NameSub,': IsNewInput=',IsNewInput

  ! Do not solve if there is no new input from GM or UA
  if(.not.IsNewInput) RETURN

  ! Check if we can have a reasonable magnetic field already
  call get_time(nStepOut=nStep)

  if(DoTest)write(*,*)NameSub,': nStep = ',nStep

  ! After the solve this input can be considered old
  IsNewInput = .false.

  ! Obtain the position of the magnetix axis
  call get_axes(time_simulation,MagAxisTiltGsmOut = ThetaTilt)

!  call get_time(tStartOut=tStart)
  tStart = TimeRamStart%time

  call time_real_to_int(tStart + time_simulation, Time_Array)

  nSolve = nSolve + 1

  if(DoTest)write(*,*) 'solve'

  ! Solve for the ionosphere potential
  call IE_solve

  if(DoTest)write(*,*) 'done with solve'

  ! Save solution (plot files) into file if required
  call IE_output

  if(DoTest)write(*,*) 'done with output'

  call IE_gather

  if(DoTest)write(*,*) 'gather done'

! Save logfile if required -- it is saved in run_ramscbe now
!  call IE_save_logfile
  
  if(DoTest)write(*,*) 'done with IE_run'

end subroutine IE_run

!=================================================================
subroutine IE_get_for_ps(Buffer_IIV, iSize, jSize, nVar)

  use ModKind
  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_ps'

  integer, intent(in)           :: iSize, jSize, nVar
  real(Real8_), intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  !NOTE: The Buffer variables must be collected to i_proc0(IE_) before return.

  write(*,*) NameSub,' -- called but not yet implemented.'

end subroutine IE_get_for_ps
!==============================================================================

subroutine IE_setnMlts(iComponent, nMLTsIn, iError)

  implicit none

  integer, intent(in)  :: iComponent, nMLTsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setnMlts'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setnMlts

!==============================================================================

subroutine IE_setnLats(iComponent, nLatsIn, iError)

  implicit none

  integer, intent(in)  :: iComponent, nLatsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setnLats'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setnLats

!==============================================================================

subroutine IE_setgrid(iComponent, MLTsIn, LatsIn, iError)

  integer, intent(in) :: iComponent
  real, intent(in) :: MLTsIn,LatsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setgrid'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setgrid

!==============================================================================

