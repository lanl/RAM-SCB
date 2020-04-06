MODULE ModRamScbRun

  implicit none

  contains
!==================================================================================================
  subroutine run_ramscb
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, iCal, nIter
    use ModRamParams,    ONLY: DoSaveFinalRestart, DoVarDt, IsComponent, &
                               verbose, reset, DoUseRam, SCBonRAMTime, RAMTie
    use ModRamTiming,    ONLY: DtsFramework, DtsMax, DtsMin, DtsNext, Dts, Dt_hI, &
                               TimeRamStart, TimeRamNow, TimeMax, TimeRamElapsed, &
                               UTs, Dt_bc, DtEfi
    use ModRamVariables, ONLY: Kp, F107, AE, dBdt, dIdt, dIbndt
    use ModScbVariables, ONLY: hICalc, SORFail
    use ModScbParams,    ONLY: method
    
    !!!! Module Subroutines and Functions
    use ModRamGSL,       ONLY: GSL_Initialize
    use ModRamFunctions, ONLY: ram_sum_pressure, RamFileName
    use ModRamIndices,   ONLY: get_indices
    use ModRamTiming,    ONLY: max_output_timestep, init_timing, finalize_timing, do_timing
    use ModRamIO,        ONLY: init_output, handle_output, ram_write_pressure, write_prefix
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
    
    implicit none
    
    logical :: triggerSCB
    real(DP) :: DtOutputMax, DtEndMax
!----------------------------------------------------------------------------
!!!!!!!! SET CURRENT TIME AND TIMESTEP
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
!!!!!!!!

!!!!!!!! UPDATES AS NEEDED
    call get_indices(TimeRamNow%Time, Kp, f107, AE)

    ! Update Boundary Flux if Dt_bc has passed
    if (abs(mod(TimeRamElapsed, Dt_bc)).le.1e-9) then
       call get_boundary_flux
    end if

    ! Update Electric Fields if DtEfi has passed
    if (abs(mod(TimeRamElapsed, DtEfi)).le.1e-9) then
       call get_electric_field
    end if
!!!!!!!

!!!!!!!! RUN RAM
    ! Broadcast current call to ram_all
    call write_prefix
    write(*,'(a, f7.2,1x,f6.2,2x,f3.1)') &
         'Calling ram_run for UTh, DTs, Kp = ', UTs/3600., Dts, Kp
    ! Call RAM for each species.
    if (DoUseRAM) call ram_run

    ! Increment and update time
    TimeRamElapsed = TimeRamElapsed + 2.0 * DTs
    TimeRamNow % Time = TimeRamStart % Time + TimeRamElapsed
    call time_real_to_int(TimeRamNow)
!!!!!!!

!!!!!!! RUN SCB
    ! Call SCB if conditions are met
    triggerSCB = .false.
    if (SCBonRAMTime) then
       if ((mod(nIter,RAMTie).eq.0).and.(nIter.gt.1)) triggerSCB = .true.
    else
       if (abs(mod(TimeRamElapsed, Dt_hI)).le.1e-9) triggerSCB = .true.
    endif
    if (triggerSCB) then
       write(*,*) ''
       call write_prefix
       write(*,'(a,F10.2)') 'Running SCB model to update B-field at time: ', TimeRamElapsed

       call ram_sum_pressure
       call scb_run(nIter)

       if ((SORFail).and.(Reset)) then
          if (verbose) write(*,*) 'Error in SCB calculation, attempting a full reset'
          call computational_domain
          call scb_run(0)
          if (SORFail) hICalc = .false.
       endif


       ! Couple SCB -> RAM
       if ((hICalc)) then ! Calculate full h's and I's if SCB was successful
          call computehI(nIter)
       else                                 ! If SCB wasn't successful use previous h's and I's
          dBdt = 0._dp                      ! which implies dXdt = 0
          dIdt = 0._dp
          dIbndt = 0._dp
       endif
       call compute3DFlux

       call write_prefix
       write(*,*) 'Finished 3D Equilibrium code.'
       write(*,*) ''
       DtsNext = DtsMin ! This is so we don't accidently take to big of a step after an SCB call
    end if
!!!!!!!

!!!!!!! CREATE OUTPUTS
    ! Do timing.
    call do_timing
    ! Check and write output, as necessary.
    call handle_output(TimeRamElapsed)
!!!!!!!

!!!!!!! FINISH TIME STEP
    iCal  = iCal + 1
    nIter = nIter + 1

    return

  end subroutine run_ramscb
!==================================================================================================
END MODULE ModRamScbRun
