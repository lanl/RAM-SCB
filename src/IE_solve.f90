subroutine IE_solve

  use IE_ModMain
  use IE_ModIO, ONLY: write_prefix
  use ModIonosphere
  use ModNumConst
!  use ModMpi
  use ModRamCouple, ONLY:  energy_fluxIono, ave_eIono, num_fluxIono, JrIono, &
                           calculate_precip_flux_Jr, DoPassJr, DoIEPrecip,&
                           dis_energy_fluxIono, dis_ave_eIono

  implicit none
  character(len=*), parameter :: NameSub = 'IE_solve'
  real(real8_)    :: CurrentSum
  integer :: iBlock,j
  integer :: nSize, iError
  real(real8_), dimension(2*Iono_nTheta-1, Iono_nPsi, 5) :: Buffer_IIV

  logical DoTest, DoTestMe
  real(real8_)    :: rIonosphere, HighLatBoundaryIm(Iono_nPsi)
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  SinThetaTilt = sin(ThetaTilt)
  CosThetaTilt = cos(ThetaTilt)

  if (DoTest) &
       write(*,'(a,i4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3)')&
       " "//NameSub//" at ",Time_Array

  if (TypeIMCouple == 'north')then

     if (.not. allocated(energy_fluxIono)) then
        allocate(         JrIono(2*Iono_nTheta-1, Iono_nPsi), &
             energy_fluxIono    (2*Iono_nTheta-1, Iono_nPsi), &
             ave_eIono          (2*Iono_nTheta-1, Iono_nPsi), &
             num_fluxIono       (2*Iono_nTheta-1, Iono_nPsi), &
             dis_energy_fluxIono(2*Iono_nTheta-1, Iono_nPsi), &
             dis_ave_eIono      (2*Iono_nTheta-1, Iono_nPsi))
        
     end if
 
   ! calculate the electron energy flux and average_energy 
   ! precipitating into the ionosphere altitude
   ! these grids are from the ionosphere. In the following routine, 
   ! the variables should be mapped down the iono. grids.
     rIonosphere = 1.0+0.0172
     call calculate_precip_flux_jr(1, 2*Iono_nTheta-1, Iono_nPsi, rIonosphere, &
          Energy_FluxIono, Ave_eIono, Num_FluxIono, Dis_Energy_FluxIono, Dis_Ave_eIono, &
          JrIono, HighLatBoundaryIm, DoIEPrecip, DoPassJr) ! only for electrons

     IsNewInput = .true.
     
     Buffer_IIV(:,:,1) = JrIono*1.0e-6            ! convert mA/m^2 to A/m^2 
     Buffer_IIV(:,:,2) = Energy_FluxIono * 1.6e-9 ! convert keV/(cm^2s) to ergs/(cm^2s)
     Buffer_IIV(:,:,3) = Ave_eIono                ! keV
     Buffer_IIV(:,:,4) = Dis_Energy_FluxIono*1.6e-9  ! convert keV/cm^2/s to  ergs/cm2/s
     Buffer_IIV(:,:,5) = Dis_Ave_eIono            ! keV
     HighLatBoundary = HighLatBoundaryIm*cRadToDeg

     iono_north_im_jr = Buffer_IIV(1:Iono_nTheta,:,1)
     iono_north_im_eflux = Buffer_IIV(1:Iono_nTheta,:,2) !ergs/cm^2s
     iono_north_im_avee  = Buffer_IIV(1:Iono_nTheta,:,3) !keV
     iono_north_im_dis_eflux = Buffer_IIV(1:Iono_nTheta,:,4) !ergs/cm^2s
     iono_north_im_dis_avee  = Buffer_IIV(1:Iono_nTheta,:,5) !keV
     iono_north_im_jr(:,Iono_nPsi) = iono_north_im_jr(:,1)
     iono_north_im_eflux(:,Iono_nPsi) = iono_north_im_eflux(:,1)
     iono_north_im_avee(:,Iono_nPsi) = iono_north_im_avee(:,1)
     iono_north_im_dis_eflux(:,Iono_nPsi) = iono_north_im_dis_eflux(:,1)
     iono_north_im_dis_avee(:,Iono_nPsi) = iono_north_im_dis_avee(:,1)

     write(*,*)'iono_north_im_eflux:',maxval(iono_north_im_eflux)
     write(*,*)'iono_north_im_dis_eflux:',maxval(iono_north_im_dis_eflux)
     write(*,*)'iono_north_im_jr:',maxval(iono_north_im_jr)

     iono_south_im_jr    = Buffer_IIV(Iono_nTheta:2*Iono_nTheta-1,:,1)
     iono_south_im_eflux = Buffer_IIV(Iono_nTheta:2*Iono_nTheta-1,:,2)
     iono_south_im_avee  = Buffer_IIV(Iono_nTheta:2*Iono_nTheta-1,:,3)
     iono_south_im_dis_eflux = Buffer_IIV(Iono_nTheta:2*Iono_nTheta-1,:,4)
     iono_south_im_dis_avee  = Buffer_IIV(Iono_nTheta:2*Iono_nTheta-1,:,5)
     iono_south_im_jr(:,Iono_nPsi) = iono_south_im_jr(:,1)
     iono_south_im_eflux(:,Iono_nPsi) = iono_south_im_eflux(:,1)
     iono_south_im_avee(:,Iono_nPsi) = iono_south_im_avee(:,1)
     iono_south_im_dis_eflux(:,Iono_nPsi) = iono_south_im_dis_eflux(:,1)
     iono_south_im_dis_avee(:,Iono_nPsi) = iono_south_im_dis_avee(:,1)

     write(*,*)'iono_south_im_eflux:',maxval(iono_south_im_eflux)
     write(*,*)'iono_south_im_dis_eflux:',maxval(iono_south_im_dis_eflux)
     write(*,*)'iono_south_im_jr:',maxval(iono_south_im_jr)
  end if

  do iBlock = 1, 2

     if(DoTest)write(*,*) 'iblock',iblock

     select case(iBlock)
     case(1) ! Northern hemisphere

        north = .true.

!        CurrentSum = sum(abs(IONO_NORTH_JR))
!        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_NORTH_JR))=', CurrentSum
!        if(CurrentSum < 1e-6)CYCLE

!        ! Add the IM currents before the conductances are calculated
!        IONO_NORTH_JR = IONO_NORTH_JR + FractionImJr*iono_north_im_jr
        IONO_NORTH_JR  = iono_north_im_jr
        call FACs_to_fluxes(conductance_model, iBlock)
        call ionosphere_conductance(IONO_NORTH_Sigma0,               &
             IONO_NORTH_SigmaH, IONO_NORTH_SigmaP,    &
             IONO_NORTH_SigmaThTh,                    &
             IONO_NORTH_SigmaThPs,                    &
             IONO_NORTH_SigmaPsPs,                    &
             IONO_NORTH_dSigmaThTh_dTheta,            &
             IONO_NORTH_dSigmaThPs_dTheta,            &
             IONO_NORTH_dSigmaPsPs_dTheta,            &
             IONO_NORTH_dSigmaThTh_dPsi,              &
             IONO_NORTH_dSigmaThPs_dPsi,              &
             IONO_NORTH_dSigmaPsPs_dPsi,              &
             IONO_NORTH_EFlux, IONO_NORTH_Ave_E,      &
             IONO_NORTH_Theta, IONO_NORTH_Psi,        &
             IONO_nTheta, IONO_nPsi,                  &
             dTheta_North, dPsi_North,                &
             conductance_model, f107_flux)

        ! Add in ionospheric currents after the conductances
        ! are calculated
!        IONO_NORTH_JR = IONO_NORTH_JR - Iono_North_Tgcm_Jr
        call ionosphere_solver(iBlock, &
             IONO_NORTH_JR,     &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs,   &
             IONO_NORTH_SigmaPsPs,                         &
             IONO_NORTH_dSigmaThTh_dTheta,                 &
             IONO_NORTH_dSigmaThPs_dTheta, &
             IONO_NORTH_dSigmaPsPs_dTheta, &
             IONO_NORTH_dSigmaThTh_dPsi,  &
             IONO_NORTH_dSigmaThPs_dPsi, &
             IONO_NORTH_dSigmaPsPs_dPsi, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             dTheta_North, dPsi_North, &
!             HighLatBoundaryIm, &
             IONO_NORTH_PHI)
        if(DoTest)then
           call write_prefix; 
           write(*,*) "Northern Cross Polar Cap Potential=",cpcp_north," kV"
        end if

        ! Calculate Currents and Boundary Conditions for North

        call ionosphere_currents(iBlock, &
             IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz,&
             IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
             IONO_NORTH_ETh,IONO_NORTH_EPs, &
             IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
             IONO_NORTH_PHI, &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
             IONO_NORTH_SigmaPsPs, &
             IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             dTheta_North, dPsi_North)
        !add joule heating for north (JouleHeating = sigmaP * E^2)
        ! Yiqun
        call ionosphere_jouleheating_ionflux(iBlock, &
             IONO_NORTH_ETh, IONO_NORTH_EPs, &
             IONO_NORTH_SigmaP, &
             IONO_NORTH_Joule,  &
             IONO_NORTH_IonNumFlux)

     case(2) ! Southern hemisphere

        north = .false.

!        CurrentSum = sum(abs(IONO_SOUTH_JR))
!        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_SOUTH_JR))=', CurrentSum
!        if(CurrentSum < 1e-6)CYCLE
!
        ! Add the IM currents before the conductances are calculated
!        IONO_SOUTH_JR = IONO_SOUTH_JR + FractionImJr*iono_south_im_jr
        IONO_SOUTH_JR = iono_south_im_jr
        call FACs_to_fluxes(conductance_model, iBlock)
        call ionosphere_conductance(IONO_SOUTH_Sigma0,               &
             IONO_SOUTH_SigmaH, &
             IONO_SOUTH_SigmaP, &
             IONO_SOUTH_SigmaThTh, &
             IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_dSigmaThTh_dTheta, &
             IONO_SOUTH_dSigmaThPs_dTheta, &
             IONO_SOUTH_dSigmaPsPs_dTheta, &
             IONO_SOUTH_dSigmaThTh_dPsi, &
             IONO_SOUTH_dSigmaThPs_dPsi, &
             IONO_SOUTH_dSigmaPsPs_dPsi, &
             IONO_SOUTH_EFlux, IONO_SOUTH_Ave_E,  &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             IONO_nTheta, IONO_nPsi,                  &
             dTheta_South, dPsi_South,                &
             conductance_model, f107_flux)
        ! Add in ionospheric currents after the conductances
        ! are calculated
        IONO_SOUTH_JR = IONO_SOUTH_JR - Iono_South_Tgcm_Jr

        call ionosphere_solver(iBlock, &
             IONO_SOUTH_JR, &
             IONO_SOUTH_SigmaThTh, &
             IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_dSigmaThTh_dTheta, &
             IONO_SOUTH_dSigmaThPs_dTheta, &
             IONO_SOUTH_dSigmaPsPs_dTheta, &
             IONO_SOUTH_dSigmaThTh_dPsi, &
             IONO_SOUTH_dSigmaThPs_dPsi, &
             IONO_SOUTH_dSigmaPsPs_dPsi, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             dTheta_South, dPsi_South, &
!             HighLatBoundaryIm, &
             IONO_SOUTH_PHI)

        if(DoTest)then
           call write_prefix; 
           write(*,*) "Southern Cross Polar Cap Potential=",cpcp_south," kV"
        end if

        ! Calculate Currents and Boundary Conditions for South

        call ionosphere_currents(iBlock, &
             IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz,&
             IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
             IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
             IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
             IONO_SOUTH_PHI, &
             IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             dTheta_South, dPsi_South)
        call ionosphere_jouleheating_ionflux(iBlock, &
             IONO_SOUTH_ETh, IONO_SOUTH_EPs, &
             IONO_SOUTH_SigmaP, &
             IONO_SOUTH_Joule,  &
             IONO_SOUTH_IonNumFlux)

     end select

  end do

end subroutine IE_solve
