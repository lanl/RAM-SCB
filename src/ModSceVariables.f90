!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================
Module ModSceVariables

  use nrtype, ONLY: DP, cDegToRad

  implicit none

  logical :: north
  logical :: UseWeimer, DoPrecond, UseInitialGuess

  real(DP) :: SinThetaTilt, CosThetaTilt, ThetaTilt

  ! From ModRamCouple
  real(DP), allocatable :: JrIono(:,:), energy_fluxIono(:,:), ave_eIono(:,:), &
                           num_fluxIono(:,:), dis_energy_fluxIono(:,:), dis_ave_eIono(:,:)

  ! From IE_Solve
  real(DP), allocatable :: HighLatBoundary(:)

  ! From IE_ModMain
  integer :: MaxIteration = 600
  real(DP) :: Hall_to_Ped_Ratio
  real(DP) :: PolarCapPedConductance = 0.25, &
              StarLightPedConductance = 0.50, &
              LatBoundary = 30.0 * cDegToRad
  real(DP) :: Tolerance = 1.e-2

  ! From ModIonosphere
  real(DP), parameter :: IONO_TOLER = 5.0e-05,     &
                         IONO_MU = 1.256637e-06,   &
                         IONO_Theta_0 = 0.0001,    &
                         IONO_Min_EFlux = 0.1e-16, & ! W/m2
                         IONO_Min_Ave_E = 0.5,     & ! keV
                         Polar_Rain = 0.1e-2         ! W/m2

  integer, parameter :: IONO_Model_No_Hall = 1,            &
                        IONO_Model_With_Hall = 2,          &
                        IONO_Model_With_Simple_Aurora = 3, &
                        IONO_Model_With_Complex_Aurora = 4

  integer  :: nThetaUsed, nThetaSolver, nX
  real(DP) :: Radius, cpcp_north, IONO_Radius, IONO_Height

  real(DP), allocatable :: PhiIono_Weimer(:,:)
  integer, allocatable  :: iHighBnd(:)

  real(DP), allocatable :: PhiOld_CB(:,:,:)

  real(DP), allocatable :: IONO_NORTH_PHI(:,:), IONO_NORTH_X(:,:), IONO_NORTH_Y(:,:), &
                           IONO_NORTH_Z(:,:), IONO_NORTH_Theta(:,:), IONO_NORTH_Psi(:,:), &
                           IONO_NORTH_Ex(:,:), IONO_NORTH_Ey(:,:), IONO_NORTH_Ez(:,:), &
                           IONO_NORTH_ETh(:,:), IONO_NORTH_EPs(:,:), IONO_NORTH_Ux(:,:), &
                           IONO_NORTH_Uy(:,:), IONO_NORTH_Uz(:,:), IONO_NORTH_UTh(:,:), &
                           IONO_NORTH_UPs(:,:), IONO_NORTH_EFlux(:,:), IONO_NORTH_Ave_E(:,:), &
                           IONO_NORTH_Sigma0(:,:), IONO_NORTH_SigmaH(:,:), IONO_NORTH_SigmaP(:,:), &
                           IONO_NORTH_SigmaThTh(:,:), IONO_NORTH_SigmaThPs(:,:), & 
                           IONO_NORTH_SigmaPsPs(:,:), IONO_NORTH_dSigmaThTh_dTheta(:,:), &
                           IONO_NORTH_dSigmaThPs_dTheta(:,:), IONO_NORTH_dSigmaPsPs_dTheta(:,:), &
                           IONO_NORTH_dSigmaThTh_dPsi(:,:), IONO_NORTH_dSigmaThPs_dPsi(:,:), &
                           IONO_NORTH_dSigmaPsPs_dPsi(:,:), SAVE_NORTH_SigmaH(:,:), &
                           SAVE_NORTH_SigmaP(:,:), IONO_NORTH_Joule(:,:)

  real(DP), allocatable :: IONO_NORTH_JR(:,:), IONO_NORTH_JTh(:,:), IONO_NORTH_JPs(:,:), &
                           IONO_NORTH_Jx(:,:), IONO_NORTH_Jy(:,:), IONO_NORTH_Jz(:,:),   &
                           IONO_NORTH_TGCM_JR(:,:), IONO_NORTH_AMIE_JR(:,:), &
                           IONO_NORTH_Fake_JR(:,:), IONO_NORTH_IonNumFlux(:,:)

  real(DP), allocatable :: iono_north_im_jr(:,:), iono_north_im_avee(:,:), iono_north_im_eflux(:,:), &
                           iono_north_im_dis_avee(:,:), iono_north_im_dis_eflux(:,:)

  real(DP), allocatable :: dTheta_North(:), dPsi_North(:)

  real(DP), allocatable :: C_A(:,:), C_B(:,:), C_C(:,:), C_D(:,:), C_E(:,:)

End Module ModSceVariables
