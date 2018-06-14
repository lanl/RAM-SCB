MODULE ModSceInit

  implicit none

  contains
!==================================================================================================
  subroutine sce_allocate

    use ModSceGrids, ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables

    use nrtype, ONLY: DP
    implicit none

    integer :: a, b

    a = Iono_nTheta
    b = Iono_nPsi

    allocate(JrIono(2*a-1,b), energy_fluxIono(2*a-1,b), ave_eIono(2*a-1,b), &
             num_fluxIono(2*a-1,b), dis_energy_fluxIono(2*a-1,b), &
             dis_ave_eIono(2*a-1,b))
    JrIono = 0.0_dp; energy_fluxIono = 0.0_dp; ave_eIono = 0.0_dp
    num_fluxIono = 0.0_dp; dis_energy_fluxIono = 0.0_dp; dis_ave_eIono = 0.0_dp

    allocate(dTheta_North(a), dPsi_North(b), HighLatBoundary(b))
    dTheta_North = 0.0_dp; dPsi_North = 0.0_dp; HighLatBoundary = 0.0_dp

    allocate(PhiIono_Weimer(a,b), iHighBnd(b))
    PhiIono_Weimer = 0.0_dp; iHighBnd = 0;

    allocate(PhiOld_CB(a,b,2))
    PhiOld_CB = 0.0_dp

    allocate(IONO_NORTH_PHI(a,b), IONO_NORTH_X(a,b), IONO_NORTH_Y(a,b), &
             IONO_NORTH_Z(a,b), IONO_NORTH_Theta(a,b), IONO_NORTH_Psi(a,b), &
             IONO_NORTH_Ex(a,b), IONO_NORTH_Ey(a,b), IONO_NORTH_Ez(a,b), & 
             IONO_NORTH_ETh(a,b), IONO_NORTH_EPs(a,b), IONO_NORTH_Ux(a,b), &
             IONO_NORTH_Uy(a,b), IONO_NORTH_Uz(a,b), IONO_NORTH_UTh(a,b), &
             IONO_NORTH_UPs(a,b), IONO_NORTH_EFlux(a,b), IONO_NORTH_Ave_E(a,b), &              
             IONO_NORTH_Sigma0(a,b), IONO_NORTH_SigmaH(a,b), IONO_NORTH_SigmaP(a,b), &
             IONO_NORTH_SigmaThTh(a,b), IONO_NORTH_SigmaThPs(a,b), &
             IONO_NORTH_SigmaPsPs(a,b), IONO_NORTH_dSigmaThTh_dTheta(a,b), &
             IONO_NORTH_dSigmaThPs_dTheta(a,b), IONO_NORTH_dSigmaPsPs_dTheta(a,b), &
             IONO_NORTH_dSigmaThTh_dPsi(a,b), IONO_NORTH_dSigmaThPs_dPsi(a,b), &
             IONO_NORTH_dSigmaPsPs_dPsi(a,b), SAVE_NORTH_SigmaH(a,b), &
             SAVE_NORTH_SigmaP(a,b), IONO_NORTH_Joule(a,b))
    IONO_NORTH_PHI = 0.0_dp; IONO_NORTH_X = 0.0_dp; IONO_NORTH_Y = 0.0_dp
    IONO_NORTH_Z = 0.0_dp; IONO_NORTH_Theta = 0.0_dp; IONO_NORTH_Psi = 0.0_dp
    IONO_NORTH_Ex = 0.0_dp; IONO_NORTH_Ey = 0.0_dp; IONO_NORTH_Ez = 0.0_dp
    IONO_NORTH_ETh = 0.0_dp; IONO_NORTH_EPs = 0.0_dp; IONO_NORTH_Ux = 0.0_dp
    IONO_NORTH_Uy = 0.0_dp; IONO_NORTH_Uz = 0.0_dp; IONO_NORTH_UTh = 0.0_dp
    IONO_NORTH_UPs = 0.0_dp; IONO_NORTH_EFlux = 0.0_dp; IONO_NORTH_Ave_E = 0.0_dp
    IONO_NORTH_Sigma0 = 0.0_dp; IONO_NORTH_SigmaH = 0.0_dp; IONO_NORTH_SigmaP = 0.0_dp
    IONO_NORTH_SigmaThTh = 0.0_dp; IONO_NORTH_SigmaThPs = 0.0_dp
    IONO_NORTH_SigmaPsPs = 0.0_dp; IONO_NORTH_dSigmaThTh_dTheta = 0.0_dp
    IONO_NORTH_dSigmaThPs_dTheta = 0.0_dp; IONO_NORTH_dSigmaPsPs_dTheta = 0.0_dp
    IONO_NORTH_dSigmaThTh_dPsi = 0.0_dp; IONO_NORTH_dSigmaThPs_dPsi = 0.0_dp
    IONO_NORTH_dSigmaPsPs_dPsi = 0.0_dp; SAVE_NORTH_SigmaH = 0.0_dp
    SAVE_NORTH_SigmaP = 0.0_dp; IONO_NORTH_Joule(a,b) = 0.0_dp
  
    allocate(IONO_NORTH_JR(a,b), IONO_NORTH_JTh(a,b), IONO_NORTH_JPs(a,b), &
             IONO_NORTH_Jx(a,b), IONO_NORTH_Jy(a,b), IONO_NORTH_Jz(a,b),   &
             IONO_NORTH_TGCM_JR(a,b), IONO_NORTH_AMIE_JR(a,b), &
             IONO_NORTH_Fake_JR(a,b), IONO_NORTH_IonNumFlux(a,b))
    IONO_NORTH_JR = 0.0_dp; IONO_NORTH_JTh = 0.0_dp; IONO_NORTH_JPs = 0.0_dp
    IONO_NORTH_Jx = 0.0_dp; IONO_NORTH_Jy = 0.0_dp; IONO_NORTH_Jz = 0.0_dp
    IONO_NORTH_TGCM_JR = 0.0_dp; IONO_NORTH_AMIE_JR = 0.0_dp
    IONO_NORTH_Fake_JR = 0.0_dp; IONO_NORTH_IonNumFlux = 0.0_dp
  
    allocate(iono_north_im_jr(a,b), iono_north_im_avee(a,b), iono_north_im_eflux(a,b), &
             iono_north_im_dis_avee(a,b), iono_north_im_dis_eflux(a,b))
    iono_north_im_jr = 0.0_dp; iono_north_im_avee = 0.0_dp; iono_north_im_eflux = 0.0_dp
    iono_north_im_dis_avee = 0.0_dp; iono_north_im_dis_eflux = 0.0_dp

    allocate(C_A(a,b), C_B(a,b), C_C(a,b), C_D(a,b), C_E(a,b))
    C_A = 0.0_dp; C_B = 0.0_dp; C_C = 0.0_dp; C_D = 0.0_dp; C_E = 0.0_dp

    return
  end subroutine sce_allocate

!==================================================================================================
  subroutine sce_deallocate

    use ModSceVariables

    implicit none

    deallocate(JrIono, energy_fluxIono, ave_eIono, &
               num_fluxIono, dis_energy_fluxIono, &
               dis_ave_eIono)

    deallocate(dTheta_North, dPsi_North, HighLatBoundary)

    deallocate(PhiIono_Weimer, iHighBnd)

    deallocate(PhiOld_CB)

    deallocate(IONO_NORTH_PHI, IONO_NORTH_X, IONO_NORTH_Y, &
               IONO_NORTH_Z, IONO_NORTH_Theta, IONO_NORTH_Psi, &
               IONO_NORTH_Ex, IONO_NORTH_Ey, IONO_NORTH_Ez, & 
               IONO_NORTH_ETh, IONO_NORTH_EPs, IONO_NORTH_Ux, &
               IONO_NORTH_Uy, IONO_NORTH_Uz, IONO_NORTH_UTh, &
               IONO_NORTH_UPs, IONO_NORTH_EFlux, IONO_NORTH_Ave_E, &              
               IONO_NORTH_Sigma0, IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
               IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
               IONO_NORTH_SigmaPsPs, IONO_NORTH_dSigmaThTh_dTheta, &
               IONO_NORTH_dSigmaThPs_dTheta, IONO_NORTH_dSigmaPsPs_dTheta, &
               IONO_NORTH_dSigmaThTh_dPsi, IONO_NORTH_dSigmaThPs_dPsi, &
               IONO_NORTH_dSigmaPsPs_dPsi, SAVE_NORTH_SigmaH, &
               SAVE_NORTH_SigmaP, IONO_NORTH_Joule)
 
    deallocate(IONO_NORTH_JR, IONO_NORTH_JTh, IONO_NORTH_JPs, &
               IONO_NORTH_Jx, IONO_NORTH_Jy, IONO_NORTH_Jz,   &
               IONO_NORTH_TGCM_JR, IONO_NORTH_AMIE_JR, &
               IONO_NORTH_Fake_JR)
 
    deallocate(iono_north_im_jr, iono_north_im_avee, iono_north_im_eflux, &
               iono_north_im_dis_avee, iono_north_im_dis_eflux)

    deallocate(C_A, C_B, C_C, C_D, C_E)

    return
  end subroutine sce_deallocate

!==================================================================================================
  subroutine sce_init

    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables, ONLY: Radius, IONO_Radius, IONO_Height, IONO_NORTH_Theta, &
                               IONO_NORTH_Psi, IONO_NORTH_X, IONO_NORTH_Y, &
                               IONO_NORTH_Z, dPsi_North, dTheta_North, &
                               Iono_Theta_0

    use nrtype, ONLY: DP, pio2_d, twopi_d, pi_d

    implicit none

    integer :: i,j
    real(DP) :: dTheta_l, dPsi_l
    !------------------------------------------------------------------------
    dTheta_l = pio2_d/(IONO_nTheta-1)
    dPsi_l   = twopi_d/(IONO_nPsi-1)

    IONO_Radius = 6378.0e+3
    IONO_Height = 110000.00


    Radius = IONO_Radius + IONO_Height

    do j = 1, IONO_nPsi
       IONO_NORTH_Theta(1,j) = IONO_Theta_0
       IONO_NORTH_Psi(1,j) = real(j-1)*dPsi_l
       IONO_NORTH_X(1,j) = sin(IONO_NORTH_Theta(1,j)) * cos(IONO_NORTH_Psi(1,j))
       IONO_NORTH_Y(1,j) = sin(IONO_NORTH_Theta(1,j)) * sin(IONO_NORTH_Psi(1,j))
       IONO_NORTH_Z(1,j) = cos(IONO_NORTH_Theta(1,j))

       do i = 2, IONO_nTheta
          IONO_NORTH_Theta(i,j) = pio2_d*real(i-1)/real(IONO_nTheta-1)
          IONO_NORTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
          IONO_NORTH_X(i,j) = sin(IONO_NORTH_Theta(i,j)) * cos(IONO_NORTH_Psi(i,j))
          IONO_NORTH_Y(i,j) = sin(IONO_NORTH_Theta(i,j)) * sin(IONO_NORTH_Psi(i,j))
          IONO_NORTH_Z(i,j) = cos(IONO_NORTH_Theta(i,j))
       end do
    end do

    !do j = 1, IONO_nPsi
    !   do i = 1, IONO_nTheta
    !      IONO_SOUTH_Theta(i,j) = pi_d - IONO_NORTH_Theta(IONO_nTheta-(i-1),j)
    !      IONO_SOUTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
    !      IONO_SOUTH_X(i,j) = sin(IONO_SOUTH_Theta(i,j)) * cos(IONO_SOUTH_Psi(i,j))
    !      IONO_SOUTH_Y(i,j) = sin(IONO_SOUTH_Theta(i,j)) * sin(IONO_SOUTH_Psi(i,j))
    !      IONO_SOUTH_Z(i,j) = cos(IONO_SOUTH_Theta(i,j))
    !   end do
    !end do

    !
    ! dPsi is going to be defined as 2*dPsi for all of the points:
    !
    do j = 2, IONO_nPsi-1
       dPsi_North(j) = IONO_NORTH_Psi(1,j+1)-IONO_NORTH_Psi(1,j-1)
    enddo
    dPsi_North(1)    = IONO_NORTH_Psi(1,2) - (IONO_NORTH_Psi(1,IONO_nPsi-1) - twopi_d)
    dPsi_North(IONO_nPsi) = dPsi_North(1)

    !
    ! dTheta is going to be defined as 2*dTheta for all but the top and bottom
    !   points.  As these locations, we switch to 1st order
    !
    do i = 2, IONO_nTheta-1
       dTheta_North(i) = IONO_NORTH_Theta(i+1,1) - IONO_NORTH_Theta(i-1,1)
    enddo
    dTheta_North(1)           = IONO_NORTH_Theta(2,1) - IONO_NORTH_Theta(1,1)
    dTheta_North(IONO_nTheta) = IONO_NORTH_Theta(IONO_nTheta,1) - IONO_NORTH_Theta(IONO_nTheta-1,1)

    !
    ! dPsi is going to be defined as 2*dPsi for all of the points:
    !
    !do j = 2, IONO_nPsi-1
    !   dPsi_South(j) = IONO_SOUTH_Psi(1,j+1)-IONO_SOUTH_Psi(1,j-1)
    !enddo
    !dPsi_South(1)         = IONO_SOUTH_Psi(1,2) - (IONO_SOUTH_Psi(1,IONO_nPsi-1) - twopi_d)
    !dPsi_South(IONO_nPsi) = dPsi_South(1)

    !
    ! dTheta is going to be defined as 2*dTheta for all but the top and bottom
    !   points.  As these locations, we switch to 1st order
    !
    !do i = 2, IONO_nTheta-1
    !   dTheta_South(i) = IONO_SOUTH_Theta(i+1,1) - IONO_SOUTH_Theta(i-1,1)
    !enddo
    !dTheta_South(1)           = IONO_SOUTH_Theta(2,1) - IONO_SOUTH_Theta(1,1)
    !dTheta_South(IONO_nTheta) = IONO_SOUTH_Theta(IONO_nTheta,1) - IONO_SOUTH_Theta(IONO_nTheta-1,1)

    return

  end subroutine sce_init

!==================================================================================================
END MODULE ModSceInit
