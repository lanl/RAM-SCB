!==============================================================================
module ModRamMain
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModTimeConvert, ONLY: TimeType
  use ModNumConst, ONLY: cPi

  implicit none
  save

!\
! New parameters for type, version, etc.
! These will make the true ModRamMain as
! modules are decentralized.
!/

! RAM version.
real, parameter :: CodeVersion = 2.0

! Kind for double precision, set here for portability.
integer, parameter :: Real8_ = selected_real_kind(14,200)
integer, parameter :: Real4_ = selected_real_kind( 6, 37)

! Standalone or Component (default = standalone)
logical :: IsComponent = .false.

! SHIELDS-RC Mode?
logical :: IsSHIELDS=.false.

! RESTART variables.
! Is this run a restart of a previous simulation?
logical :: IsRestart = .false., DoSaveFinalRestart=.true.
real(kind=Real8_) :: DtRestart=3600.0, TimeRestart=0.0
character(len=*), parameter :: PathRestartIn = 'IM/restartIN'
character(len=*), parameter :: PathRestartOut= 'IM/restartOUT'

! Use plasmasphere density from Rasmussen model, coupled to MSIS and IRI
! For plane_scb options, see ModRamPl_Ne.f90
logical :: DoUsePlane_SCB = .false.

! Include SCB?  Default is yes!
logical :: DoUseScb = .true.

! Include induced E field?  Default is no!
logical :: UseEfind = .false.

! Include wave-particle diffusion?  Default is no!
logical :: DoUseWPI = .false.

! Use BAS wave-particle diffusion?  Default is no!
logical :: DoUseBASdiff = .false.

! Do use Kp-dependent diffusion coefficients?  Default is no!
logical :: DoUseKpDiff = .false.

! Flags for initialization of energy grid and initial conditions
logical, public :: DoInitOnly = .false.
logical :: DoUseVAPini = .false.

! Force multi-species boundary conditions through the use of 
! multiple files for BCS flux values.
logical :: DoMultiBcsFile = .false.

! Dump 2D fluxes, B field, currents, pressure anis to netcdf?  Default is no!
logical :: DoWriteFlux=.false.
logical :: DoWriteB=.false.
logical :: DoWriteCur=.false.
logical :: DoWritePres=.false.

! Dump 2D plasmaspheric density to netcdf?  Default is no!
logical :: DoWriteDensity=.false.

! Include electrons?  Default is yes!
logical :: electrons = .true.

! IO units for persistent files:
integer :: iUnit7

!\
! New parameters for Stand Alone mode.
!/
integer :: TimeMax                                 ! Simulation max in seconds.
integer :: MaxIter                                 ! Simulation max iterations
logical :: DoVarDt = .true.                        ! Use variable timestep.
character(len=100) :: StringTest=''                ! List of tested subroutines
character(len=200) :: StrRamDescription='None'     ! Descript. of simulation

!\
! Time object for convienient creation of time strings.
!/
type(TimeType) :: TimeRamNow
type(TimeType) :: TimeRamStart
logical :: IsStarttimeSet=.false.  ! True if #STARTTIME used in standalone.

!\
! Parameters (to be set dynamically later.)
! Originally set in ram_ion1.h
!/

integer, parameter :: &
     NR  = 20, &   ! grid points in radial direction 
     NT  = 25, &   ! grid points in local time direction
     NE  = 35, &   ! number of energy bins
     NS  = 4,  &   ! number of species
     NPA = 72, &   ! grid points in pitch angle dimension
     NL  = 35, &   ! grid points in radial for PLANE.F
     NLT = 48      ! grid points in local time for PLANE.F
integer :: &       ! --Below are now dynamically set.--
     NEL = 36, &   ! Energy bins in boundary file (for SWMF or LANL)
     NEL_prot = 36, &   ! Energy bins in bound file for protons (SWMF or LANL)
     NBD = 288, &  ! Number of boundary data per file = (12/hour)*(24 hours)
     NTL = 49      ! MLT bins for SWMF or LANL grid (25 for LANL, 49 for SWMF)

! For B coupling, we must request additional ghostcells.
integer, parameter :: nRextend = NR+3

integer, parameter :: Slen = 55, ENG = 45, NCF = 5, NCO = 5

real(kind=Real8_), parameter :: RadiusMin = 1.75, RadiusMax = 6.5


!\
! Variables from Main.f and ram_ion.f that are better in Modules.
!/

real(Real4_) :: FLUX(4,NR,NT,NE,NPA)=0.0

character(len=4) :: NameBoundMag ! Originally 'boundMag' in Main.f
integer          :: iCal         ! Tracks iteration number, kind of.
integer          :: nIter=0      ! Tracks iteration, including restart.
integer          :: iPressure, iDomain

! Timing variables:
real(kind=Real8_):: UTs, DTs = 5.0, DTsNext = 5.0, FracCFL=0.8, DtsMin=1.0, &
     DtsMax = 100.0, TimeRamElapsed = 0., DtsFramework=1000.0
real(kind=Real8_),parameter:: DtEfi=300.0, Dt_bc=300.0, Dt_hI=300.0
! This is the file output interval, fixed at 1 hour (3600. seconds)
! Uncomment the next line to couple DtOutput with Dt_bc
real(kind=Real8_):: DtOutput = 3600.
!real(kind=Real8_):: DtOutput = 12.*Dt_bc
! Explanations:
! FracCFL-- the percent of the CFL timestep that we will use.
! If the CFL timestep limit is 1.0s, the timestep taken would be 0.8s.
! Dts, DtsNext, DtsMax, and DtsFramework -- These control the half-timestep
! taken by RAM.  Dts is the half-timestep that will be taken and is chosen
! from DtsNext (set by the CFL limit), DtsMax (maximum timestep allowed),
! and DtsFramework (set by the SWMF to synchronize couplings).  Because it
! is a half timestep (due to the time-splitting scheme), a full time step
! is 2x Dts.
! Dt_hI is the couple frequency between RAM and SCB, default is 300s
! DtEfi is the couple frequency between RAM and E field, default is 300s
! Dt_bc is the time step to update the boundary conditions, default is 300s

! Pressures:
real(kind=Real8_), dimension(NR,NT) :: &
     PPerH, PParH, PPerO, PParO, PPerHe, PParHe, PPerE, PParE, PAllSum, &
     NAllSum, DensO, DensH, DensHe, HPAllSum, OPAllSum, HePAllSum, ePAllSum, &
     HNAllSum, ONAllSum, HeNAllSum, PparSum

! Limiter Beta (scales MC Limiter from min-mod=1 to superbee=2; default=1.5)
real(kind=Real8_) :: BetaLim = 1.5

! File paths to fix hardcoding of path issues.
character(len=*), parameter :: PathRamOut = 'IM/output/ram/'
character(len=*), parameter :: PathScbOut = 'IM/output/scb/'
character(len=*), parameter :: PathSwmfOut= 'IM/output_swmf/'
character(len=*), parameter :: PathRamIn  = 'IM/input_ram/'
character(len=*), parameter :: PathScbIn  = 'IM/input_scb/'

!\
! CFUN Block.
!/

real(kind=Real8_) ::        &
     F2(4,NR,NT,NE,NPA),    &
     FGEOS(4, NT, NE, NPA), &
     BNES(NR+1, NT),        &
     dBdt(NR+1,NT),         &
     dIdt(NR+1,NT,NPA),     &
     dIbndt(NR+1,NT,NPA)

real(kind=Real8_) :: &
     VT(NR+1, NT), &
     T,            &
     EIR(NR+1, NT),       & 
     EIP(NR+1, NT)


logical :: DoAnisoPressureGMCoupling = .false.

!\
! PARAMS Block.
!/

real(kind=Real8_) :: DT

integer :: S

character(len=3) :: ST0, ST1
character(len=2) :: ST2
character(len=19):: ST7


!\
! CONST Block.
!/

real(kind=Real8_) ::   &
     ME   = 7.9E15,    &   ! Magnetic moment of the earth [T*m3]
     RE   = 6.371E6,   &   ! Earth's radius [m]
     HMIN = 2E5,       &   ! Altitude of dense atmosphere [m]
     MP   = 1.673E-27, &   ! Mass of H+ [kg]
     Q    = 1.602E-19, &   ! Elementary charge [C]
     CS   = 2.998E8,   &   ! Speed of light [m/s]
     PI   = cPi            ! pi used to be 3.141592654
! These variables removed from latest version.
!     FLUXFACT             ! convert distribution; see SetFluxfact
!     VTO  = 5.E4, &	   ! Voll-Stern quiescent conv potential [V]
!     COLS = 8.547, &	   ! L*, or co-latitude of polar cap boundary


!\
! CAR Block.
!/

real(kind=Real8_) :: &
     RMAS(NS), &
     V(NE),    &
     VBND(NE)
! EKEV(k): center of energy channel 'k'.  WE(k): Width of energy channel 'k'.
! EBND(k): upper limit of E-channel 'k'.  DE(k): EKEV(k+1) - EKEV(k)
real(kind=Real8_)                        :: IR1, DL1, DR, dPhi, IP1
real(kind=Real8_), dimension(NR+1)       :: LZ, Z
real(kind=Real8_), dimension(nRextend)   :: GridExtend
integer,           dimension(NR)         :: UPA
real(kind=Real8_), dimension(NT)         :: PHI, LT, MLT
real(kind=Real8_), dimension(NE)         :: WE, DE, EKEV, EBND 
real(kind=Real8_), dimension(NPA)        :: MU, DMU, WMU
real(kind=Real8_), dimension(NR+1, Slen) :: BE

   ! Mass number of e-, H+, He+ and O+
   ! should be moved to constants.
real(kind=Real8_), dimension(4) :: M1  = (/5.4462E-4,1.,4.,16./)


!\
! CFAC Block.
!/

real(kind=Real8_) :: CONF1, CONF2

real(kind=Real8_) :: &
     CEDR(NE, NPA), CIDR(NE, NPA), &
     GREL(NE), GRBND(NE), FACGR(NE), &
     FFACTOR(NR, NE, NPA)


!\
! CE Block.
!/

real(kind=Real8_) :: ACHAR(NR,NT,NE,NPA), ATLOS(NR,NE)
real(kind=Real8_) :: HDNS(NR+1,NT,NPA)


!\
! CDR Block.
!/

real(kind=Real8_), dimension(NR+1,NT,NPA) :: FNHS, FNIS, BOUNHS, BOUNIS
real(kind=Real8_), dimension(NR,NE)       :: P2, EDOT

real(kind=Real8_) :: A1C, VR(NR), P1(NR), MUDOT(NR,NPA), PHIOFS

integer, dimension(NR,NT) :: ISGN


!\
! CCOUL Block.
!/

real(kind=Real8_), dimension(NE, NPA) :: &
     COULE, COULI, GTAE, ATA, BTA, GTA, GTAI


!\
! CINIT Block.
!/

real(kind=Real8_), dimension(NS,NR):: ENERN, ENERD, LECN, LECD
real(kind=Real8_) :: FACTOR, KP=0, F107, NECR(NL,0:NLT), XNE(NR,NT) 
real :: IG(3), RZ(3), rsn, IAPO(7)   ! Ap-index, Rs-index

     ! Unused variables.
real(kind=Real8_) :: ECOF(NR), WCD(NR,NT), YIND(NR,NE)


!\
! CANIS Block.
!/

real(kind=Real8_) :: &
     ATEW(NR,NT,NE,NPA), ATAW(NR,NT,NE,NPA), &
     PPERT(NS,NR,NT), PPART(NS,NR,NT), &
     EPP(NE), ERNH(NE)

real(kind=Real8_), dimension(NR+1) :: BFC


!\
! CWPI Block.
!/

integer, parameter :: &  ! grid points for chorus diffusion coefficients:
     NR_Dxx  = 20, &   ! in radial direction
     NT_Dxx  = 25, &   ! in local time direction
     NE_Dxx  = 45, &   ! number of energy bins
     NPA_Dxx = 72      ! number of pitch angles

real(kind=Real8_) :: PAbn(NPA), AMLA(Slen)
real(kind=Real8_) :: ZRpabn(NR,NPA,Slen),fpofc(NCF),fpcho(NCO)
real(kind=Real8_) :: NDVVJ(NR,ENG,NPA,NCF), NDAAJ(NR,ENG,NPA,NCF), ENOR(ENG)
real(kind=Real8_) :: WALOS1(NR,NE), WALOS2(NR,NE), WALOS3(NR,NE)
real(kind=Real8_) :: CDAAR(NR,NT,NE,NPA),CDEER(NR,NT,NE,NPA),ECHOR(ENG),BDAAR(NR,NT,ENG,NPA)
real(kind=Real8_) :: ATAW_emic(NR,NT,NE,NPA),ATAC(NR,NT,NE,NPA), &
     ATEC(NR,NT,NE,NPA),ATMC(NR,NT,NE,NPA)

! Arrays to store bounce-averaged BAS diffusion coefficients.

real(kind=Real8_), dimension(NR_Dxx) :: RCHOR_Dxx
real(kind=Real8_), dimension(NT_Dxx) :: TCHOR_Dxx
real(kind=Real8_), dimension(NE_Dxx) :: ECHOR_Dxx
real(kind=Real8_), dimension(NPA_Dxx) :: PACHOR_Dxx

! KP-dependent diffusion coefficients
integer,parameter :: NKpDiff = 5
integer,parameter,dimension(NKpDiff) :: Kp_chorus = (/ 0,1,2,3,4 /)
real(kind=Real8_) :: CDAAR_chorus(NR_Dxx,NT_Dxx,NE_Dxx,NPA_Dxx,NKpDiff)

! Variables for interpolating to regular grid
integer, parameter :: Nx = NPA  ! Grid points for PA in mixed diffusion solver
integer, parameter :: Ny = 300  ! Grid points for energy    "     "

integer, dimension(2,Nx) :: x_ind
integer, dimension(2,Ny) :: y_ind

integer, dimension(2,NPA) :: PA_ind
integer, dimension(2,NE) :: E_ind

real(kind=Real8_), dimension(Nx) :: x1,Jcbn1
real(kind=Real8_), dimension(Ny) :: y1,Jcbn2

real(kind=Real8_), dimension(2,2,Ny,Nx) :: xycoeff
real(kind=Real8_), dimension(2,2,NE,NPA) :: EPAcoeff

DATA fpcho/1.5,2.5,5.,7.5,10./    ! fpe/fce for chorus waves


!\
! CELOS Block.
!/

real(kind=Real8_), dimension(NS) :: &
     LSDR,  &
     LSCHA, &
     LSATM, &
     LSCOE, &
     LSCSC, &
     LSWAE, &
     ELORC, &
     SETRC


!\
! RDIFF2 Block.
! WLO Block.
!/
     ! Not used variables.

!\
! RESTART block
!/

character(len=6) :: event = 'de_flt'
character(len=4) :: boundary, electric, NameDistrib
real(kind=Real8_) :: hrSim

end Module ModRamMain
!==============================================================================
