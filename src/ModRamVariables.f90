!============================================================================
!    Copyright (c) 2016, Los Alamos National Laboratory
!    All rights reserved.
!============================================================================

Module ModRamVariables

  use ModRamSpecies, ONLY: SpeciesType
  use ModRamGrids, ONLY: dR, dPhi
  use nrtype, ONLY: DP

  implicit none

  ! Species variables
  type(SpeciesType), allocatable :: species(:)
  real(DP), allocatable :: composition(:)

  integer, allocatable :: outsideMGNP(:,:)

  ! ModRamPlasmasphere Variables
  real(DP), ALLOCATABLE :: flux_volume(:,:), uL(:,:), uT(:,:), tau(:,:), &
                           smLon(:,:), smLat(:,:), xRAM(:,:,:), yRAM(:,:,:), zRAM(:,:,:)

! UNKNOWN VARIABLES
  real(DP), ALLOCATABLE :: XNE(:,:)
!

  real :: IG(3), RZ(3), rsn, IAPO(7)   ! Ap-index, Rs-index

! Main RAM Variables (Pressures, Fluxes, and hI variables)
  real(DP), ALLOCATABLE :: F2(:,:,:,:,:), FLUX(:,:,:,:,:), PPerH(:,:), PParH(:,:), &
                           PPerO(:,:), PParO(:,:), PPerHe(:,:), PParHe(:,:), &
                           PPerE(:,:), PParE(:,:), PAllSum(:,:), PParSum(:,:), &
                           PPerT(:,:,:), PParT(:,:,:), dIdt(:,:,:), dIbndt(:,:,:), &
                           HDNS(:,:,:), FNHS(:,:,:), FNIS(:,:,:), BOUNHS(:,:,:), &
                           BOUNIS(:,:,:), BNES(:,:), dBdt(:,:), dHdt(:,:,:), &
                           ODNS(:,:,:), NDNS(:,:,:)
  real(DP), ALLOCATABLE :: NAllSum(:,:), DensO(:,:), DensH(:,:), DensHe(:,:), HPAllSum(:,:), &
                           OPAllSum(:,:), HePAllSum(:,:), ePAllSum(:,:), HNAllSum(:,:), &
                           ONAllSum(:,:), HeNAllSum(:,:)

! ModRamInit variables
  real(DP), ALLOCATABLE :: RMAS(:), V(:,:), VBND(:,:), GREL(:,:), GRBND(:,:), &
                           FACGR(:,:), EPP(:,:), ERNH(:,:), UPA(:), WE(:), DE(:), &
                           EKEV(:), EBND(:), PHI(:), LT(:), MLT(:), MU(:), DMU(:), &
                           WMU(:), PAbn(:), LZ(:), RLZ(:), AMLA(:), BE(:,:), &
                           GridExtend(:), ZRPabn(:,:,:), FFACTOR(:,:,:,:), PA(:)
  real(DP) :: PHIOFS, IR1, DL1, MDR, IP1, CONF1, CONF2, RFACTOR

! ModRamWPI variables
  real(DP), ALLOCATABLE :: WALOS1(:,:), WALOS2(:,:), WALOS3(:,:), fpofc(:), &
                           NDVVJ(:,:,:,:), NDAAJ(:,:,:,:), ENOR(:), ECHOR(:), &
                           BDAAR(:,:,:,:), CDAAR(:,:,:,:),&
                           Daa_emic_h(:,:,:,:), Daa_emic_he(:,:,:,:), &
                           EKEV_emic(:), fp2c_emic(:), Ihs_emic(:,:,:), Ihes_emic(:,:,:)
  integer, parameter :: NKpDiff = 5, &
                        NR_Dxx  = 20, &
                        NT_Dxx  = 25, &
                        NE_Dxx  = 45, &
                        NPA_Dxx = 72  
  integer, dimension(5) :: Kp_Chorus = (/0,1,2,3,4/) 
  real(DP) :: CDAAR_Chorus(NR_Dxx,NT_Dxx,NE_Dxx,NPA_Dxx,NKpDiff)
  real(DP) :: RCHOR_Dxx(NR_Dxx), TCHOR_Dxx(NT_Dxx), ECHOR_Dxx(NE_Dxx), &
              PACHOR_Dxx(NPA_Dxx)

! ModRamLoss variables
  real(DP), ALLOCATABLE :: ATLOS(:,:,:), CHARGE(:,:,:,:,:), &
                           COULE(:,:,:), COULI(:,:,:), ATA(:,:,:), GTA(:,:,:), &
                           GTAE(:,:,:), GTAI(:,:,:), CEDR(:,:,:), CIDR(:,:,:)

! ModRamIndices variables
  character(len=4)   :: NameIndexSource = 'file'
  integer :: nRawKp, nRawF107, nRawAE
  integer, parameter :: kptime(8) = (/1, 4, 7, 10, 13, 16, 19, 22/)
  real(DP) :: KP, F107
  integer  :: AE
  real(DP), allocatable :: timeKp(:),timeF107(:),timeAE(:), rawKp(:),rawF107(:)
  integer,  allocatable :: rawAE(:)

! ModRamEField variables
  real(DP), ALLOCATABLE :: VT(:,:), EIR(:,:), EIP(:,:), VTOL(:,:), VTN(:,:)
  real(DP) :: TOLV

! ModRamBoundary variables
  real(DP), ALLOCATABLE :: FGEOS(:,:,:,:)
  logical :: IsInitialized = .false.
  ! Start time offset in seconds.
  real(DP) :: timeOffset
  ! Date of file:
  character(len=8) :: StringFileDate
  ! Values that are read and stored from file:
  real(DP), allocatable :: flux_SIII(:,:,:,:), fluxLast_SII(:,:,:), eGrid_SI(:,:), &
                           avgSats_SI(:,:), tGrid_SI(:,:), lGrid_SI(:,:)

! ModRamDrift variables
  !real(DP), ALLOCATABLE :: P1(:), VR(:), P2(:,:), EDOT(:,:), MUDOT(:,:), &
  !                         CDriftR(:,:,:,:), CDriftP(:,:,:,:), sgnDriftR(:,:,:,:), &
  !                         CDriftE(:,:,:,:), CDriftMu(:,:,:,:)
  real(DP), ALLOCATABLE :: DTDriftR(:), DTDriftP(:), DTDriftE(:), DTDriftMu(:)
  real(DP) :: FracCFL = 0.8

! ModRamRun variables
  real(DP), ALLOCATABLE :: SETRC(:), ELORC(:), LSDR(:), LSCHA(:), LSATM(:), LSCOE(:), &
                           LSCSC(:), LSWAE(:), XNN(:,:), XND(:,:), LNCN(:,:), LNCD(:,:), &
                           LECN(:,:), LECD(:,:), ENERN(:,:), ENERD(:,:), ATEW(:,:,:,:), &
                           ATAW(:,:,:,:), ATAC(:,:,:,:), ATEC(:,:,:,:), ATMC(:,:,:,:), &
                           NECR(:,:), ESUM(:), NSUM(:),&
                           ATAW_emic_h(:,:,:,:), ATAW_emic_he(:,:,:,:), &
                           r_curvEq(:,:), zeta1Eq(:,:), zeta2Eq(:,:), FLC_coef(:,:,:,:,:)

End Module ModRamVariables
