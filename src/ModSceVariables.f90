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
       num_fluxIono(:,:), dis_energy_fluxIono(:,:), dis_ave_eIono(:,:),&
       diff_fluxIono(:,:,:)

  ! From IE_Solve
  real(DP), allocatable :: HighLatBoundary(:), HighLatBoundaryIM(:)

  ! From IE_ModMain
  integer :: MaxIteration = 600
  real(DP) :: StarLightPedConductance = 0.50, LatBoundary = 30.0 * cDegToRad
  real(DP) :: Tolerance = 1.e-2
  integer  :: Conductance_Model

  ! From ModIonosphere
  real(DP), parameter :: IONO_TOLER = 5.0e-05,     &
       IONO_theta_0 = 0.001, &
       IONO_Min_EFlux = 0.1e-16, & ! W/m2
                         IONO_Min_Ave_E = 0.5

  integer  :: nThetaUsed, nThetaSolver, nX
  real(DP) :: Radius, cpcp, IONO_Radius, IONO_Height

  real(DP), allocatable :: PhiIono_Weimer(:,:)
  integer, allocatable  :: iHighBnd(:)

  real(DP), allocatable :: PhiOld_CB(:,:)

  real(DP), allocatable :: IONO_NORTH_PHI(:,:), IONO_NORTH_X(:,:), IONO_NORTH_Y(:,:), &
                           IONO_NORTH_Z(:,:), IONO_NORTH_Theta(:,:), IONO_NORTH_Psi(:,:), &
                           IONO_NORTH_EFlux(:,:), IONO_NORTH_Ave_E(:,:), &
                           IONO_NORTH_EFlux_diff(:,:,:), &
                           IONO_NORTH_Sigma0(:,:), IONO_NORTH_SigmaH(:,:), IONO_NORTH_SigmaP(:,:), &
                           IONO_NORTH_SigmaThTh(:,:), IONO_NORTH_SigmaThPs(:,:), & 
                           IONO_NORTH_SigmaPsPs(:,:), IONO_NORTH_dSigmaThTh_dTheta(:,:), &
                           IONO_NORTH_dSigmaThPs_dTheta(:,:), IONO_NORTH_dSigmaPsPs_dTheta(:,:), &
                           IONO_NORTH_dSigmaThTh_dPsi(:,:), IONO_NORTH_dSigmaThPs_dPsi(:,:), &
                           IONO_NORTH_dSigmaPsPs_dPsi(:,:), SAVE_NORTH_SigmaH(:,:), &
                           SAVE_NORTH_SigmaP(:,:), IONO_NORTH_Joule(:,:)

  real(DP), allocatable :: IONO_NORTH_JR(:,:), IONO_NORTH_IonNumFlux(:,:)

  real(DP), allocatable :: iono_north_im_jr(:,:), iono_north_im_avee(:,:), iono_north_im_eflux(:,:), &
       iono_north_im_dis_avee(:,:), iono_north_im_dis_eflux(:,:), &
       iono_north_im_eflux_diff(:,:,:)

  real(DP), allocatable :: dTheta_North(:), dPsi_North(:)

  real(DP), allocatable :: C_A(:,:), C_B(:,:), C_C(:,:), C_D(:,:), C_E(:,:)
  !\
  ! GLOW related paremater
  !/
  integer, parameter :: nzGLOW = 112

  integer :: nSolve = 0
  character(len=10):: NameSolver='bicgstab' ! Name of krylov solver
  logical :: DoSaveLogfile = .False.
  logical :: UsePreconditioner = .True. 
  logical :: DoUseFullSpec = .False.
  logical :: DoSaveGLOWConductivity = .False.
  
End Module ModSceVariables
