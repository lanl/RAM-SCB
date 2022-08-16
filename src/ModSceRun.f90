!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModSceRun

  implicit none

  contains
!==================================================================================================
  subroutine sce_run

    use ModRamTiming,    ONLY: TimeRamElapsed
    use ModSceVariables, ONLY: ThetaTilt

    use CON_axes, ONLY: get_axes
 
    implicit none
    !----------------------------------------------------------------------------

    call get_axes(TimeRamElapsed ,MagAxisTiltGsmOut = ThetaTilt)  
    call sce_solve
  
  end subroutine sce_run
!==================================================================================================
  subroutine sce_solve

    use ModScbVariables, ONLY: pjconst
    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables
    use ModSceIono, ONLY: ionosphere_conductance,  ionosphere_solver, ionosphere_fluxes
    use ModRamSce, ONLY: calculate_precip_flux_jr
    use ModRamVariables, ONLY: species
    use ModRamGrids, ONLY: nS
    use nrtype, ONLY: DP, cDegtoRad, cRadtoDeg
    implicit none

    real(DP) :: rIonosphere
    integer :: IS, S
    !--------------------------------------------------------------------------
    do IS=1,nS
       if (species(IS)%s_name .eq. 'Electron')then ! only the electron precipitation is used for SCE calculation
          S=IS
          exit
       end if
    end do

    SinThetaTilt = sin(ThetaTilt)
    CosThetaTilt = cos(ThetaTilt)

    ! calculate the electron energy flux and average_energy 
    ! precipitating into the ionosphere altitude
    ! these grids are from the ionosphere. In the following routine, 
    ! the variables should be mapped down the iono. grids.
    rIonosphere = 1.0+0.0172

    call calculate_precip_flux_jr(S, 2*Iono_nTheta-1, Iono_nPsi, rIonosphere, &
                                  Energy_FluxIono, Ave_eIono, Num_FluxIono, &
                                  Dis_Energy_FluxIono, Dis_Ave_eIono, &
                                  JrIono, HighLatBoundaryIM, Diff_FluxIono) ! only for electrons

    JrIono              = JrIono*1.0e-6               ! convert mA/m^2 to A/m^2 
    Energy_FluxIono     = Energy_FluxIono*1.6e-9      ! convert keV/(cm^2s) to ergs/(cm^2s)
    Ave_eIono           = Ave_eIono                   ! keV
    Dis_Energy_FluxIono = Dis_Energy_FluxIono*1.6e-9  ! convert keV/cm^2/s to ergs/cm2/s
    Dis_Ave_eIono       = Dis_Ave_eIono               ! keV
    HighLatBoundary     = HighLatBoundaryIM*cRadToDeg

    iono_north_im_jr        = JrIono(1:Iono_nTheta,:)
    iono_north_im_eflux     = Energy_FluxIono(1:Iono_nTheta,:)
    iono_north_im_avee      = Ave_eIono(1:Iono_nTheta,:)
    iono_north_im_dis_eflux = Dis_Energy_FluxIono(1:Iono_nTheta,:)
    iono_north_im_dis_avee  = Dis_Ave_eIono(1:Iono_nTheta,:)
    iono_north_im_eflux_diff = Diff_FluxIono(1:Iono_nTheta,:,:)  ! differential flux !/cm^2/s/sr/keV
    
    iono_north_im_jr(:,Iono_nPsi)        = iono_north_im_jr(:,1)
    iono_north_im_eflux(:,Iono_nPsi)     = iono_north_im_eflux(:,1)
    iono_north_im_avee(:,Iono_nPsi)      = iono_north_im_avee(:,1)
    iono_north_im_dis_eflux(:,Iono_nPsi) = iono_north_im_dis_eflux(:,1)
    iono_north_im_dis_avee(:,Iono_nPsi)  = iono_north_im_dis_avee(:,1)
    iono_north_im_eflux_diff(:,Iono_nPsi,:) = iono_north_im_eflux_diff(:,1,:)
    
    IONO_NORTH_JR  = iono_north_im_jr
    call ionosphere_fluxes(Conductance_Model)
    call ionosphere_conductance(IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
                                IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, IONO_NORTH_SigmaPsPs, &
                                IONO_NORTH_dSigmaThTh_dTheta, IONO_NORTH_dSigmaThPs_dTheta, &
                                IONO_NORTH_dSigmaPsPs_dTheta, IONO_NORTH_dSigmaThTh_dPsi, &
                                IONO_NORTH_dSigmaThPs_dPsi, IONO_NORTH_dSigmaPsPs_dPsi, &
                                IONO_NORTH_EFlux, IONO_NORTH_Ave_E, IONO_NORTH_EFlux_Diff, &
                                IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_nTheta, IONO_nPsi, &
                                dTheta_North, dPsi_North, IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z,&
                                Conductance_Model)
    call ionosphere_solver(IONO_NORTH_JR, IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                           IONO_NORTH_SigmaPsPs, IONO_NORTH_dSigmaThTh_dTheta, &
                           IONO_NORTH_dSigmaThPs_dTheta, IONO_NORTH_dSigmaPsPs_dTheta, &
                           IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
                           IONO_NORTH_dSigmaPsPs_dPsi, IONO_NORTH_Theta, IONO_NORTH_Psi, &
                           dTheta_North, dPsi_North, IONO_NORTH_PHI)
    write(*,*) "Northern Cross Polar Cap Potential=",cpcp," kV"

    return
  end subroutine sce_solve
!==================================================================================================
END MODULE ModSceRun

