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
  real(kind=Real8_), dimension(NR,NPA,Slen) :: ZRPabn
  real(kind=Real8_), dimension(NS,NR,NE,NPA) :: FFACTOR
  real(kind=Real8_) :: PHIOFS, IR1, DL1, MDR, dPhi, IP1, CONF1, CONF2, RFACTOR

! ModRamWPI variables
  integer, parameter :: NKpDiff = 5, &
                        Kp_Chorus = (/0,1,2,3,4/), &
                        NR_Dxx  = 20, &
                        NT_Dxx  = 25, &
                        NE_Dxx  = 45, &
                        NPA_Dxx = 72   
  real(kind=Real8_) :: CDAAR_Chorus(NR_Dxx,NT_Dxx,NE_Dxx,NKpDiff)
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
  real(kind=Rael8_) :: KP, F107
  real(kind=Real8_), allocatable :: timeKp(:),timeF107(:),rawKp(:),rawF107(:)

! ModRamEField variables
  real(kind=Real8_), dimension(NR+1,NT) :: VT, EIR, EIP, VTOL, VTN
  real(kind=Real8_) :: TOLV

! ModRamBoundary variables
  real(kind=Real8_), dimension(NS,NT,NE,NPA) :: FGEOS

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



real(kind=Real8_) :: DT, &
                     UTs, &

! Pressures
real(kind=Real8_), dimension(NR,NT)             :: &
    PAllSum,           &
    NAllSum,           &
    DensO,             &
    DensH,             &
    DensHe,            &
    HPAllSum,          &
    OPAllSum,          &
    HePAllSum,         &
    ePAllSum,          &
    HNAllSum,          &
    ONAllSum,          &
    HeNAllSum,         &
    PparSum


!\
! CAR Block.
!/

real(kind=Real8_), dimension(nRextend)          :: &
    GridExtend
! EKEV(k): center of energy channel 'k'.  WE(k): Width of energy channel 'k'.
! EBND(k): upper limit of E-channel 'k'.  DE(k): EKEV(k+1) - EKEV(k)

!\
! CFAC Block.
!/

real(kind=Real8_), dimension(NE,NPA)            :: &
    CEDR,              &
    CIDR


real(kind=Real8_)                               :: &
    A1C,               &

!\
! CCOUL Block.
!/

real(kind=Real8_), dimension(NE,NPA)            :: &
    COULE,             &
    COULI,             &
    GTAE,              &
    ATA,               &
    BTA,               &
    GTA,               &
    GTAI

!\
! CINIT Block.
!/

real(kind=Real8_), dimension(NS,NR)             :: &
    ENERN,             &
    ENERD,             &
    LECN,              &
    LECD
real                                            :: &
    IG(3),             &
    RZ(3),             &
    rsn,               &
    IAPO(7)   ! Ap-index, Rs-index

!\
! CANIS Block.
!/

real(kind=Real8_), dimension(NR+1)              :: &
    BFC

!\
! CWPI Block.
!/


real(kind=Real8_) :: &
    fpcho(NCO), &
    CDEER(NR,NT,NE,NPA), &

! Variables for interpolating to regular grid
integer, dimension(2,Nx)                        :: x_ind
integer, dimension(2,NyE)                       :: y_ind

integer, dimension(2,NPA)                       :: PA_ind
integer, dimension(2,NE)                        :: E_ind

real(kind=Real8_), dimension(Nx)                :: x1,Jcbn1
real(kind=Real8_), dimension(NyE)               :: y1,Jcbn2

real(kind=Real8_), dimension(2,2,NyE,Nx)        :: xycoeff
real(kind=Real8_), dimension(2,2,NE,NPA)        :: EPAcoeff

DATA fpcho/1.5,2.5,5.,7.5,10./    ! fpe/fce for chorus waves

!\
! RDIFF2 Block.
! WLO Block.
!/
     ! Not used variables.

!\
! SCB Block
!/

real(kind=Real8_), dimension(npsi)              :: &
    xEqMidNew,         &
    radEqMidNew,       &
    g,                 &
real(kind=Real8_), dimension(152,202,48)        :: &
    alpha_cyl,         &
    beta_cyl
real(kind=Real8_)                               :: &
    mjac,              &
    aisw2,             &
    xmx,               &
    omega,             &
    delta,             &
    eps,               &
    xmin,              &
    xmaj,              &
    epsb,              &
    aguess,            &
    xplmin,            &
    xplmax,            &
    bcfactor,          &
    bzimf,             &
    psimax,            &
    psimin,            &
    psimm,             &
    decreaseConvAlpha, &
    decreaseConvPsi,   &
    start_time,        &
    end_time,          &
    pressurequot,      &
    errorAlphaPrev,    &
    errorPsiPrev
integer                                         :: &
    tsygcorrect,       &
    nmaximum,          &
    iCorrectedPressure,&
    iPressureChoice,   &
    lconv,             &
    inest,             &
    nitry,             &
    iUseSavedData,     &
    iteration,         &
    isw1,              &
    itout,             &
    iderr,             &
    itooff,            &
    jorgn,             &
    icurm,             &
    nthe2,             &
    isym,              &
    nisave,            &
    isw2,              &
    kMax,              &
    numProc,           &
    rank,              &
    iConvGlobal,       &
    iSWMF,             &
    iAlphaMove,        &
    iPsiMove,          &
    iST3,              &
    iDay,              &
    iHour,             &
    iHourAbs,          &
    iMin,              &
    iHourBegin,        &
    iCallRoutine,      &
    k9,                &
    k21,               &
    iCompDomain
