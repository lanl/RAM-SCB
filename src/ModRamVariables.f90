Module ModRamVariables

  use ModRamMain,  ONLY: Real8_
  use ModRamGrids, ONLY: NS, NE, NT, NR, NPA, Slen, ENG, NCF, NL, NRextend

  implicit none

! UNKNOWN VARIABLES
  real(kind=Real8_) :: XNE(NR,NT)
!

!!!!! MAIN RAM VARIABLES (Pressures, Fluxes, and hI variables)
  real(kind=Real8_), dimension(NS,NR,NT,NE,NPA) :: F2, FLUX
  real(kind=Real8_), dimension(NR,NT) :: PPerH, PParH, PPerO, PParO, PPerHe, &
                                         PParHe, PPerE, PParE, PAllSum, PParSum
  real(kind=Real8_), dimension(NS,NR,NT) :: PPerT, PParT
  real(kind=Real8_), dimension(NR+1,NT,NPA) :: dIdt, dIbndt, HDNS, FNHS, FNIS, &
                                               BOUNHS, BOUNIS
  real(kind=Real8_), dimension(NR+1,NT) :: BNES, dBdt
!!!!!

! ModRamInit variables
  real(kind=Real8_), dimension(NS) :: RMAS
  real(kind=Real8_), dimension(NS,NE) :: V, VBND, GREL, GRBND, FACGR, EPP, ERNH
  real(kind=Real8_), dimension(NR) :: UPA
  real(kind=Real8_), dimension(NE) :: WE, DE, EKEV, EBND
  real(kind=Real8_), dimension(NT) :: PHI, LT, MLT
  real(kind=Real8_), dimension(NPA) :: MU, DMU, WMU, PAbn
  real(kind=Real8_), dimension(NR+1) :: LZ, RLZ
  real(kind=Real8_), dimension(Slen) :: AMLA
  real(kind=Real8_), dimension(NR+1, Slen) :: BE
  real(kind=Real8_), dimension(NRExtend) :: GridExtend
  real(kind=Real8_), dimension(NR,NPA,Slen) :: ZRPabn
  real(kind=Real8_), dimension(NS,NR,NE,NPA) :: FFACTOR
  real(kind=Real8_) :: PHIOFS, IR1, DL1, MDR, dPhi, IP1, CONF1, CONF2, RFACTOR

! ModRamWPI variables
  integer, parameter :: NKpDiff = 5, &
                        NR_Dxx  = 20, &
                        NT_Dxx  = 25, &
                        NE_Dxx  = 45, &
                        NPA_Dxx = 72  
  integer, dimension(5) :: Kp_Chorus = (/0,1,2,3,4/) 
  real(kind=Real8_) :: CDAAR_Chorus(NR_Dxx,NT_Dxx,NE_Dxx,NPA_Dxx,NKpDiff)
  real(kind=Real8_), dimension(NR,NE)          :: WALOS1, WALOS2, WALOS3
  real(kind=Real8_), dimension(NCF)            :: fpofc
  real(kind=Real8_), dimension(NR,ENG,NPA,NCF) :: NDVVJ, NDAAJ
  real(kind=Real8_), dimension(ENG)            :: ENOR, ECHOR
  real(kind=Real8_), dimension(NR,NT,ENG,NPA)  :: BDAAR
  real(kind=Real8_), dimension(NR,NT,NE,NPA)   :: CDAAR
  real(kind=Real8_), dimension(NR_Dxx)  :: RCHOR_Dxx
  real(kind=Real8_), dimension(NT_Dxx)  :: TCHOR_Dxx
  real(kind=Real8_), dimension(NE_Dxx)  :: ECHOR_Dxx
  real(kind=Real8_), dimension(NPA_Dxx) :: PACHOR_Dxx

! ModRamLoss variables
  real(kind=Real8_), dimension(NR,NT,NE,NPA) :: ACHAR
  real(kind=Real8_), dimension(NR,NE)        :: ATLOS

! ModRamIndices variables
  character(len=4)   :: NameIndexSource = 'file'
  character(len=200) :: NameIndexFile = 'RamIndices.txt'
  character(len=200) :: NameOmniFile  = 'omni.txt'
  integer :: nRawKp, nRawF107
  integer, parameter :: kptime(8) = (/1, 4, 7, 10, 13, 16, 19, 22/)
  real(kind=Real8_) :: KP, F107
  real(kind=Real8_), allocatable :: timeKp(:),timeF107(:),rawKp(:),rawF107(:)

! ModRamEField variables
  real(kind=Real8_), dimension(NR+1,NT) :: VT, EIR, EIP, VTOL, VTN
  real(kind=Real8_) :: TOLV

! ModRamBoundary variables
  real(kind=Real8_), dimension(NS,NT,NE,NPA) :: FGEOS
  logical :: IsInitialized = .false.
  ! Start time offset in seconds.
  real(kind=Real8_) :: timeOffset
  ! Date of file:
  character(len=8) :: StringFileDate
  ! Values that are read and stored from file:
  real(kind=Real8_), allocatable :: flux_SIII(:,:,:,:), fluxLast_SII(:,:,:), eGrid_SI(:,:), avgSats_SI(:,:)
  real(kind=Real8_) :: lGrid_SI(2,0:25)

! ModRamDrift variables
  real(kind=Real8_), dimension(NR)     :: P1, VR
  real(kind=Real8_), dimension(NR,NE)  :: P2, EDOT
  real(kind=Real8_), dimension(NR,NPA) :: MUDOT
  real(kind=Real8_) :: DTDriftR, DTDriftP, DTDriftE, DTDriftMu
  real(kind=Real8_) :: FracCFL = 0.8

! ModRamRun variables
  real(kind=Real8_), dimension(NS) :: SETRC, ELORC, LSDR, LSCHA, LSATM, LSCOE, &
                                      LSCSC, LSWAE
  real(kind=Real8_), dimension(NS,NR) :: XNN, XND, LNCN, LNCD, LECN, LECD, &
                                         ENERN, ENERD
  real(kind=Real8_), dimension(NR,NT,NE,NPA) :: ATEW, ATAW, ATAC, ATEC, ATMC, &
                                                ATAW_emic
  real(kind=Real8_) :: NECR(NL,0:48)

End Module ModRamVariables
