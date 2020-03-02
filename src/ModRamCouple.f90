!==============================================================================
module ModRamCouple
!    This module contains variables and routines for coupling RAM-SCB to
!    the SWMF: BATS-R-US and Ridley_Serial/RIM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModRamMain, ONLY: Real8_, PathRamOut
  use ModRamGrids, ONLY: NR, NE, NT, NPA, nRextend
  use ModRamTiming, ONLY: TimeRamElapsed, TimeRamNow
  use ModRamVariables, ONLY: KP, F107, EKEV, MU, PHI, LZ, GridExtend

  use ModRamParams

  implicit none

  public :: set_type_mhd
  public :: generate_flux
  public :: sort_mhd_lines
  public :: generate_field

  ! Can be single, multiS, or multiF to indicate MHD run.
  character(len=6), public :: TypeMhd='single'

  ! Indices of total and species density and temperature.
  ! These are NEVER assumed to increase robustness.
  integer, public :: &
       TotalRho_ = -1, TotalPres_= -1, &
       RhoH_ = -1, PresH_ = -1, &
       RhoHe_= -1, PresHe_= -1, &
       RhoO_ = -1, PresO_ = -1, &
       RhoSw_= -1, PresSw_= -1, &
       Ppar_ = -1, Bmin_  = -1, &
       Bx_=-1, By_=-1, Bz_=-1,  &
       Ux_=-1, Uy_=-1, Uz_=-1

  ! Coupling P to SWMF:
  ! Container for MHD density and pressure at outer boundary.
  ! Indices are (dens:pres,Local Time,species[all, h, he, O])
  real(kind=Real8_), public, allocatable :: MhdDensPres_VII(:,:,:)
  real(kind=Real8_), public, allocatable :: FluxBats_IIS(:,:,:)
  real(kind=Real8_), public, allocatable :: FluxBats_anis(:,:,:,:)

  ! Variables for E-field coupling with MHD (E=UxB):
  ! Only need equatorial values.
  real(kind=Real8_), public, allocatable :: ETotal_DII(:,:,:)
  
  ! Variables for B-field coupling with MHD:
  ! Tracing info is stored in MhdLines_IIV(nLines, nPointsMax, nVarLine).
  ! Lines progress from 1 to nLines such that a line passing through radius
  ! index iR and MLT index iT is stored in iLine = 2*((nRextend)*(iT-1) + iR)-1
  ! (for trace from equatorial plane to northern hemisphere) and in
  ! iLine = 2*((nRextend)*(iT-1) + iR) (to southern hemisphere).  Each line has
  ! up to nPointsMax; iEnd(iLine) has index of where to stop.  As MhdLines_IIV
  ! is filled, points are added between the MHD inner boundary and R=1RE to
  ! give all lines nPointsMax points.
  integer, public, parameter   :: nPointsMax = 500
  integer, public :: nRadSWMF, nRadSwmfVar, nRadSwmfVarInner, &
       nLonSWMF, nLinesSWMF, nPoints
  real(kind=Real8_), public, allocatable :: iEnd(:)
  real(kind=Real8_), public, allocatable :: MhdLines_IIV(:,:,:)
  real(kind=Real8_), public, allocatable :: xSWMF(:,:,:), ySWMF(:,:,:), &
       zSWMF(:,:,:), LatSWMF(:,:), LonSWMF(:,:)
  real(kind=Real8_), public, allocatable :: xEqSWMF(:,:), yEqSWMF(:,:), &
       pEqSWMF(:,:), nEqSWMF(:,:)
  real(kind=Real8_), public, allocatable :: uEqSWMF_DII(:,:,:), bEqSWMF_DII(:,:,:)

  ! Updated MHD-> SCB coupling variables: Blines has the X, Y, and Z coordinates along each
  ! magnetic field line.  Example: the X-coords of the line passing through RAM cell iR, iT
  ! is stored in Blines_DIII(1, iR, iT, :).  Tracing happens from southern (index 1) to
  ! northern hemisphere in the final dimension.
  ! Corresponding logical variable denotes if line is close (.true.) or open (.false.).
  real(kind=Real8_), public, allocatable :: Blines_DIII(:,:,:,:)
  logical,           public, allocatable :: IsClosed_II(:,:)
  
  ! Variables for coupling to SWMF-IE component:
  integer, public :: nIePhi=0, nIeTheta=0
  real(kind=Real8_), public, allocatable :: SwmfIonoPot_II(:,:)
  real(kind=Real8_), public, allocatable :: SwmfPot_II(:,:)

  real(kind=real8_), allocatable :: JrIono(:,:), energy_fluxIono(:,:), ave_eIono(:,:), &
                                    num_fluxIono(:,:), dis_energy_fluxIono(:,:), dis_ave_eIono(:,:)

contains

!==============================================================================
  subroutine RAMCouple_Allocate

    use ModRamGrids, ONLY: nS, nT, nE, nRExtend, nR

    implicit none

    ALLOCATE(MhdDensPres_VII(3,nT,0:nS), &
             FluxBats_IIS(nE, nT, nS), &
             FluxBats_anis(nE,nPa,nT,nS), &
             iEnd(2*(nRExtend)*nT), &
             xEqSWMF(nRExtend,nT-1), &
             yEqSWMF(nRExtend,nT-1), &
             pEqSWMF(nRExtend,nT-1), &
             nEqSWMF(nRExtend,nT-1), &
             SwmfPot_II(nR+1, nT), &
             uEqSWMF_DII(3,nRExtend,nT-1), &
             bEqSWMF_DII(3,nRExtend,nT-1), &
             ETotal_DII(2,nRExtend-2,nT), &
             IsClosed_II(nRextend,nT))

    SWMFPot_II = 0.
    FluxBats_anis = 0.
    FluxBats_IIS = 0.
    xEqSWMF = 0.
    yEqSWMF = 0.
    pEqSWMF = 0.
    nEqSWMF = 0.
    uEqSWMF_DII = 0.
    bEqSWMF_DII = 0.
    IsClosed_II  = .true.
    
  end subroutine RAMCouple_Allocate

!==============================================================================
  subroutine RAMCouple_Deallocate

    implicit none
    
    ! Items that may not have been allocated:
    if(allocated(MHDLines_IIV)) deallocate(MHDLines_IIV)

    DEALLOCATE(MhdDensPres_VII, FluxBats_IIS, iEnd, &
         FluxBats_anis, xEqSWMF, yEqSWMF, pEqSWMF, nEqSWMF, SwmfPot_II, &
         uEqSWMF_DII, bEqSWMF_DII, ETotal_DII)

  end subroutine RAMCouple_Deallocate

  !===========================================================================
  subroutine set_type_mhd(NameVarIn, nVarIn)
    ! Parse NameVar to determine if the info from the MHD code 
    ! is single/ multispecies or multifluid.  Check for H, He, 
    ! and O+ in simulation if multi.  If those cannot be found, 
    ! revert to single species.
    

    implicit none

    character(len=*), intent(in) :: NameVarIn
    integer, intent(in)          :: nVarIn

    character(len=200) :: NameVar
    integer :: i, nChar

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'set_type_mhd'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    NameVar = trim(NameVarIn)

    ! Parse NameVar to locate variables of interest.
    do i=1, nVarIn
       nChar=index(NameVar, ' ')
       if (nChar .eq. -1) nChar = len(trim(NameVar))
       select case ( trim(NameVar(1:nChar)) )
       case('rho')
          TotalRho_ = i
       case('Rho')
          TotalRho_ = i
       case('p')
          TotalPres_ = i
       case('P')
          TotalPres_ = i
       case('RhoH')
          RhoH_ = i
       case('HpRho')
          RhoH_ = i
       case('RhoHe')
          RhoHe_ = i
       case('RhoO')
          RhoO_ = i
       case ('OpRho')
          RhoO_ = i
       case('SwRho')
          RhoSw_ = i
       case('HpP')
          PresH_ = i
       case('HeP')
          PresHe_ = i
       case('OpP')
          PresO_ = i
       case('SwP')
          PresSw_ = i
       case('Ppar')
          Ppar_ = i
       case('ppar')
          Ppar_ = i
       case('Ux')
          Ux_=i
       case('Uy')
          Uy_=i
       case('Uz')
          Uz_=i  
       case('Bx')
          Bx_=i
       case('By')
          By_=i
       case('Bz')
          Bz_=i
       end select
       NameVar = trim( NameVar(nChar + 1:len(NameVar)) )
    end do

    if ( all((/ RhoH_,  RhoHe_, RhoO_  /) .ne. -1) ) TypeMhd = 'multiS'
    if ( all((/ PresH_, PresHe_,PresO_ /) .ne. -1) ) TypeMhd = 'multiF'
    if ( all((/ PresSw_,PresH_ ,PresO_ /) .ne. -1) ) TypeMhd = 'mult3F'
    if ( Ppar_ .ne. -1) TypeMhd = 'anisoP'

    if(DoTestMe) then
       write(*,*) 'RAM_SCB: NameVarIn = ', NameVarIn
       write(*,*) 'RAM_SCB: TypeMhd = ', TypeMhd
       write(*,*) 'RAM_SCB: MHD Indices = ', &
            'total = ', TotalRho_, TotalPres_, &
            'H+    = ', RhoH_, PresH_, &
            'He+   = ', RhoHe_,PresHe_, &
            'O+    = ', RhoO_, PresO_, &
            'Sw H+ = ', RhoSw_,PresSw_, &
            'Bx, By, Bz = ', Bx_, By_, Bz_, &
            'Ux, Uy, Uz = ', Ux_, Uy_, Uz_
    end if

  end subroutine set_type_mhd

  !===========================================================================
  subroutine sort_mhd_lines(nVarIn, nPointIn, BufferLine_VI)
    ! Sort MHD field lines & associated state variables into manageable array.
    ! Extend field lines in preparation for calculation of Euler potential
    ! surfaces.  Create equatorial variables.
   
    use ModRamConst, ONLY: b0dip
    use ModRamFunctions,ONLY: get_dipole_trace

    implicit none

    integer,           intent(in) :: nVarIn, nPointIn
    real(kind=Real8_), intent(in) :: BufferLine_VI(nVarIn, nPointIn)

    integer           :: i, j, iRad, iLon, iLine, iPointBuff, iPointLine
    real(kind=Real8_) :: xRam, yRam, Corot=7.272205E-5!rads/sec

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='sort_mhd_lines'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Clean slate:
    iPointBuff = 1
    iPointLine = 1
    MhdLines_IIV = 0.0

    if(DoTestMe) write(*,*)'IM:'//NameSub//'nLinesSWMF = ', nLinesSWMF

    ! Loop through lines, find the number of points in each line (stored in
    ! iEnd), determine missing lines, etc.  Sort from big 2D array to
    ! 3D array.  2D array has all tracing info for all field lines on one
    ! dimension; 3D array has info for each tracing segregated.
    NewLine: do iLine=1, nLinesSWMF
       ! Skip missing lines (i.e., inside MHD boundary) by checking the
       ! line number reported by the MHD tracer:
       if (iLine<BufferLine_VI(1, iPointBuff)) then
          iEnd(iLine) = -1
          cycle NewLine
       end if

       ! Grab all points on current line:
       do while(BufferLine_VI(1,iPointBuff)<iLine)
          ! Do not exceed buffer:
          if (iPointBuff == nPointIn) exit
          iPointBuff = iPointBuff + 1
       enddo
       do while(BufferLine_VI(1,iPointBuff)==iLine)
          ! Extract values into 3D array:
          MhdLines_IIV(iLine, iPointLine, :) = BufferLine_VI(:, iPointBuff)
          ! Do not exceed buffer:
          if (iPointBuff == nPointIn) exit
          ! Do not exceeed max points
          if (iPointLine == nPointsMax) exit
          ! Continue along line: 
          iPointBuff = iPointBuff + 1
          iPointLine = iPointLine + 1
       end do
       
       ! Check for very long lines or open lines.
       if (iPointLine>nPointsMax-5) then
          write(*,*) NameSub, ' ERROR: nPointsMax exceeded.'
          ! This line is too long.  Why?
          if (abs(MhdLines_IIV(iLine, iPointLine-1, 5)) > 5.0) write(*,*) &
               'Open field line?  Z=',MhdLines_IIV(iLine,iPointLine-1,5)
          !call CON_stop(NameSub//' - nPointsMax exceeded.')
          iPointLine = 0
       endif

       ! Save length of current line, reset.
       iEnd(iLine) = iPointLine - 1
       iPointLine = 1
    end do NewLine

    MhdLines_IIV(:,:,3:5) = MhdLines_IIV(:,:,3:5)/6378100.0 !XYZ in R_E.

    ! Now, extend lines from MHD inner boundary to ionosphere.
    ! Fill in missing lines w/ dipole lines.
    nPoints = int(maxval(iEnd)+5) ! Find the longest line in array + extra points
    iLon=1; iRad=1
    do iLine=1, nLinesSWMF
       ! If line is missing...
       if(iEnd(iLine)<=0)then
          xRam = GridExtend(iRad)*cos(phi(iLon))
          yRam = GridExtend(iRad)*sin(phi(iLon))
          if (sqrt(xRam**2 + yRam**2) > 3.0) then
             iEnd(iLine) = 0
          else 
             MhdLines_IIV(iLine,1,3) = xRam
             MhdLines_IIV(iLine,1,4) = yRam
             MhdLines_IIV(iLine,1,5) = 0.0
             call get_dipole_trace((/xRam, yRam, 0.0_Real8_/), nPoints-1, &
                                    MhdLines_IIV(iLine,2:nPoints,3), &
                                    MhdLines_IIV(iLine,2:nPoints,4), &
                                    MhdLines_IIV(iLine,2:nPoints,5), &
                                    MhdLines_IIV(iLine,2:nPoints,Bx_), &
                                    MhdLines_IIV(iLine,2:nPoints,By_), &
                                    MhdLines_IIV(iLine,2:nPoints,Bz_))
             ! Southern hemisphere?  Use simple symmetry.
             if (mod(iLine, 2) == 0) then 
                MhdLines_IIV(iLine,:,5)   = -1.0*MhdLines_IIV(iLine,:,5)
                MhdLines_IIV(iLine,:,Bx_) = -1.0*MhdLines_IIV(iLine,:,Bx_)
                MhdLines_IIV(iLine,:,By_) = -1.0*MhdLines_IIV(iLine,:,By_)
                !MhdLines_IIV(iLine,:,Bz_) = -1.0*MhdLines_IIV(iLine,:,Bz_)
             endif
             ! Fill equatorial values with simple conditions:
             ! Magnetic field: use dipole field, SI units are Tesla.
             MhdLines_IIV(iLine,1,Bx_:By_) = 0.0
             !MhdLines_IIV(iLine,1,Bz_)     = 7.84E15 / GridExtend(iRad)**3
             MhdLines_IIV(iLine,1,Bz_)     = b0dip/10**9 / GridExtend(iRad)**3
             ! Flow velocity: corotation
             MhdLines_IIV(iLine,1,Ux_) = corot*GridExtend(iRad)*sin(phi(iLon))
             MhdLines_IIV(iLine,1,Uy_) = corot*GridExtend(iRad)*cos(phi(iLon))
             MhdLines_IIV(iLine,1,Uz_) = 0.0
             iEnd(iLine) = nPoints
          endif
       else
          i = int(iEnd(iLine))
          !write(*,*)'Extending!'
          !write(*,*)'Adding ', nPoints-i, 'points.'
          !write(*,'(a,1f9.6)') '...from ', sqrt(MhdLines_IIV(iLine,i,3)**2 &
          !                                     +MHDLines_IIV(iLine,i,4)**2 &
          !                                     +MHDLines_IIV(iLine,i,5)**2)
          call get_dipole_trace((/MhdLines_IIV(iLine,i,3), &
                                 MhdLines_IIV(iLine,i,4),MhdLines_IIV(iLine,i,5)/), nPoints-i, &
                                 MhdLines_IIV(iLine, i+1:nPoints, 3), &
                                 MhdLines_IIV(iLine, i+1:nPoints, 4), &
                                 MhdLines_IIV(iLine, i+1:nPoints, 5), &
                                 MhdLines_IIV(iLine,i+1:nPoints,Bx_), &
                                 MhdLines_IIV(iLine,i+1:nPoints,By_), &
                                 MhdLines_IIV(iLine,i+1:nPoints,Bz_))
          if (mod(iLine, 2) == 0) then
             MhdLines_IIV(iLine,i+1:nPoints,5)   = -1.0*MhdLines_IIV(iLine,i+1:nPoints,5)
             MhdLines_IIV(iLine,i+1:nPoints,Bx_) = -1.0*MhdLines_IIV(iLine,i+1:nPoints,Bx_)
             MhdLines_IIV(iLine,i+1:nPoints,By_) = -1.0*MhdLines_IIV(iLine,i+1:nPoints,By_)
             !MhdLines_IIV(iLine,i+1:nPoints,Bz_) = -1.0*MhdLines_IIV(iLine,i+1:nPoints,Bz_)
          endif
          !write(*,'(a,1f9.6)') 'Got to:',sqrt(MhdLines_IIV(iLine,nPoints,3)**2 &
          !                                   +MHDLines_IIV(iLine,nPoints,4)**2 &
          !                                   +MHDLines_IIV(iLine,nPoints,5)**2)
          iEnd(iLine) = nPoints
       end if
       ! Keep track of our lon&rad.
       if(mod(iLine, 2) == 0 )iRad=iRad+1
       if(iRad > nRextend) then
          iLon = iLon+1
          iRad = 1
       end if
    end do

    ! Save equatorial values, working backwards towards missing values.
    ! Note that the repeated longitude (0 and 2*pi) is skipped the 2nd time.
    ! This section is for debug file writing and is not currently leveraged by
    ! SCB.  In the future, this section should be checked as indices do not
    ! line up correctly.
    
    ! Save equatorial values.  EqSWMF values are written to SCB output for
    ! later reference but not used in calculation.
    do j=1, nT-1
       do i=1, nRextend
          ! Find corresponding trace.
          iLine = 2*((nRextend)*(j-1) + i)-1 
          ! Collect scalar equatorial values (iLine, 1st trace point, variable):
          xEqSWMF(i,j) = MhdLines_IIV(iLine,1,3)
          yEqSWMF(i,j) = MhdLines_IIV(iLine,1,4)
          pEqSWMF(i,j) = MhdLines_IIV(iLine,1,TotalPres_)*1.0E9    !Pa->nPa
          nEqSWMF(i,j) = MhdLines_IIV(iLine,1,TotalRho_ )*1.67E-23 !amu->cm3
          ! Collect vector equatorial values (always sequential in list):
          bEqSWMF_DII(:,i,j) = MhdLines_IIV(iLine,1,Bx_:Bz_)
          uEqSWMF_DII(:,i,j) = MhdLines_IIV(iLine,1,Ux_:Uz_)
       end do
    end do

    ! Calculate total UxB electric field in equatorial plane.
    !iRad=nRextend-2 ! Do not include ghost cells from coupling.
    !ETotal_DII(1,:,1:nT-1) = uEqSWMF_DII(2,2:iRad,:)*bEqSWMF_DII(3,2:iRad,:) &
    !     -uEqSWMF_DII(3,2:iRad,:)*bEqSWMF_DII(2,2:iRad,:) ! X-component
    !
    !ETotal_DII(2,:,1:nT-1) = uEqSWMF_DII(1,2:iRad,:)*bEqSWMF_DII(3,2:iRad,:) &
    !     -uEqSWMF_DII(3,2:iRad,:)*bEqSWMF_DII(1,2:iRad,:) ! Y-component
    
  end subroutine sort_mhd_lines
  !===========================================================================
  subroutine generate_flux
    ! Convert MHD density and pressure into flux.
    ! Methodology depends on MHD: single, multispecies, or multifluid.
    use ModRamGrids,     ONLY: nS
    use ModRamVariables, ONLY: species
    use ModIoUnit, ONLY: UnitTmp_
    use ModConst,  ONLY: cProtonMass, cElectronCharge, cPi

    implicit none

    real, parameter :: cMass2Num = 1.0E6 * cProtonMass

    real(kind=Real8_) :: MassConv
    real(kind=Real8_) :: ratioHeH, ratioOH, fracH, fracHe, fracO
    real(kind=Real8_) :: factor, eCent, MhdPpar(4), MhdPper(4)
    real(kind=Real8_) :: factor1, kappa, gamma1, gamma2, gamma3
    integer :: iS, iT, iE, iPa
    character(len=100) :: NameFile

    ! Test Variables:
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'generate_flux'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Convert mass density and pressure into number density and temp (ev).
    ! Separate species depending on type of simulation.    
    select case(TypeMhd)
    case('single')
       if (DoTestMe) then
          write(*,*) 'Fractions used: ', species(:)%s_name
          write(*,*) species(:)%s_comp
       end if

       do iT=1, nT
          ! Convert kg/m^-3 to total #/cm^-3 for given species content.
          MassConv = 0.0
          do iS = 1, nS
             MassConv = MassConv + species(iS)%s_mass*species(iS)%s_comp
          enddo
          MhdDensPres_VII(1,iT,0) = MhdDensPres_VII(1,iT,0) / &
               ( cMass2Num * MassConv )
          ! Convert P to temp in eV.
          MhdDensPres_VII(2,iT,0) = MhdDensPres_VII(2,iT,0) / &
               ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,iT,0) )
       end do

       do iS = 1, nS
          ! Set other species number density.
          MhdDensPres_VII(1,:,iS) = species(iS)%s_comp * MhdDensPres_VII(1,:,0)
          ! All other temps the same in single fluid.
          MhdDensPres_VII(2,:,iS) = MhdDensPres_VII(2,:,0)
       enddo
    case('anisoP')
       call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
    case('multiS')
       call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
    case('multiF')
       call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
    case('mult3F')
       call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
    case default
       call CON_stop(NameSub//' TypeMhd not recognized: '//TypeMhd)
    end select

    kappa = 3.
    gamma1 = 6.0    !gamma(kappa+1)  
    gamma2 = 1.3293 !gamma(kappa-0.5)  
    gamma3 = 0.886  ! gamma(1.5)
    factor1 = gamma1/gamma2/gamma3*sqrt(2*cPi)/(2*kappa-3)**1.5

    ! Convert from moments to fluxes:
    do iE=1, nE
       eCent = 1000.0 * Ekev(iE) ! Center of energy bin in eV.
       factor=4.0E6*Ekev(iE)     ! Unit conversion factor * eCent in KeV.
       do iT=1, nT
          select case(TypeMhd)
          case('single')
             if (NameDistrib .eq. 'MAXW') then
                ! Assuming a maxwellian distribution about MHD temperature, divide
                ! densities into energy bins and calculate fluxes.
                do iS = 1, nS
                   FluxBats_IIS(iE,iT,iS) = 1./sqrt(species(iS)%s_mass) &
                                          *exp(-1.0*eCent/MhdDensPres_VII(2,iT,iS)) &
                                          *factor*MhdDensPres_VII(1,iT,iS) &
                                          *(MhdDensPres_VII(2,iT,iS)/1000.0)**(-1.5)
                   if (MhdDensPres_VII(1,iT,0) < 0) then
                      FluxBats_IIS(iE,iT,iS) = 0.0
                   endif
                enddo

             elseif(NameDistrib .eq.'KAPA')then
                ! assuming a kappa distribution
                do iS = 1, nS
                   FluxBats_IIS(iE,iT,iS) = 1./sqrt(species(iS)%s_mass) &
                                          *(1+1.0*eCent/((kappa-1.5)*MhdDensPres_VII(2,iT,iS)))**(-(kappa+1)) &
                                          *factor1*factor*MhdDensPres_VII(1,iT,iS) &
                                          *(MhdDensPres_VII(2,iT,iS)/1000.0)**(-1.5)
                   if (MhdDensPres_VII(1,iT,0) < 0) then
                      FluxBats_IIS(iE,iT,iS) = 0.0
                   endif
                enddo               
             else
                call CON_stop(NameSub//' NameDistribution not recognized: '//NameDistrib)
             end if ! NameDistrib if
          case('anisoP')
             call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
          case('multiS')
             call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
          case('multiF')
             call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
          case('mult3F')
             call CON_stop(NameSub//' TypeMhd not currently supported: '//TypeMhd)
          end select ! TypeMhd select
       end do ! nT Loop
    end do ! nE Loop
    
    ! Remove ridiculously small values.
    where(FluxBats_IIS .lt. 1.0E-30) FluxBats_IIS = 0.0

    if(DoTestMe)call write_FluxGM
    
    ! With the flux array stored in the module, coupling is complete.

  end subroutine generate_flux

  !===========================================================================
  subroutine write_FluxGM
    ! Quickly write out flux from GM to files for stand-alone RAM.
    use ModIoUnit,      ONLY: UnitTmp_

    implicit none
    
    integer :: iS, iT, nFile, iPa
    character(len=100) :: NameFile
    character(len=35)  :: StringFormat

    !------------------------------------------------------------------------
    nFile = nint(TimeRamElapsed/300.0)
    write(StringFormat,'(a,i2,a)') '(f8.2, 1x, f7.2, ', nE, '(1x, 1pE11.4))'
    ! Write one file for each species.
    do iS=1,4
       write(NameFile,'(a,i4.4,a,i1.1,a)')  &
            PathRamOut//"sep_05_",nFile,"_",iS,".swf"
       open(UnitTmp_, FILE=NameFile, STATUS='replace')
       
       write(UnitTmp_, *) &
            'Flux file for TimeRam, iSpecies = ', TimeRamElapsed, iS
       do iT=1, nT
          write(UnitTmp_,StringFormat) 1.0, iT-1.0, FluxBats_IIS(:,iT,iS)
       end do
       write(UnitTmp_,*)'Necessary Footer.'
       close(UnitTmp_)
    end do
  end subroutine write_FluxGM

  !===========================================================================
  function divMaxwellian(eMin, eMax, Dens, Temp)
    ! Enter start and stop of energy window (in eV), total density of plasma,
    ! and total plasma temperature (again, in eV) to receive new density based 
    ! on energy window for a Maxwellian particle distribution.
    
    use ModConst,   ONLY: cPi

    implicit none

    ! Arguments:
    real(kind=Real8_) :: divMaxwellian
    real(kind=Real8_), intent(in) :: Temp, Dens
    real(kind=Real8_), intent(in) :: eMin, eMax

    real(kind=Real8_) :: dE, factor1, factor2, energy
    integer           :: i
    !------------------------------------------------------------------------
    ! Set starting and delta energy.
    energy = eMin
    dE = (eMax - eMin) / 1000.0
 
    ! Do Riemann sum of Maxwellian distribution for
    ! given energy window.
    divMaxwellian = 0.0
    do i=1, 1000
       ! Handle exponent part of Maxwellian.
       ! Prevent underflows in the exponent.
       factor1 = (-1 * energy) / Temp
       if (factor1 .lt. -750.0) then
          factor1 = 0.0
       else 
          factor1 = exp(factor1)
       end if

       ! Calculate non-exponent portion of Maxwellian.
       factor2 = 2.0 * sqrt( energy / (cPi * temp**3) )

       ! Add this portion to the density sum.
       divMaxwellian = divMaxwellian + dens * factor1 * factor2 * dE

       energy = energy + dE
       
    end do

    return
  end function divMaxwellian

!==================================================================================================
  subroutine generate_field(xpsiin, xpsiout)
    !!!! Module Variables
    USE ModRamVariables, ONLY: Kp, MLT, outsideMGNP
    use ModRamTiming,    ONLY: TimeRamNow
    use ModRamGrids,     ONLY: nRExtend, nT, nR
    USE ModScbParams,    ONLY: Symmetric, constZ, constTheta, method
    USE ModScbGrids,     ONLY: npsi, nthe, nzeta
    USE ModScbVariables, ONLY: psi, r0Start, bf, nThetaEquator, rhoVal, zetaVal, &
                               x, y, z, thetaVal, zetaVal, left, right, chiVal, &
                               xzero3, psiVal, alphaVal, psiin, psiout, psitot, SORFail
    !!!! Module Subroutines/Functions
    use ModRamGSL,       ONLY: GSL_Interpolation_1D
    use ModRamFunctions, ONLY: RamFileName
    use ModScbEuler,     ONLY: mappsi, psifunctions, maptheta
    !!!! Share Modules
    USE ModIoUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d, twopi_d


    implicit none
    REAL(DP), intent(out) :: xpsiin, xpsiout

    character(len=200) :: FileName
    integer  :: i, j, k, ii, ir, il(1), nSWMF, GSLerr
    real(DP) :: r0, rr, b0, rf, xf, yf, zf, t0, tf, xpl, xpsitot, psis, rc(1)

    real(DP), dimension(0:1000) :: distance, xTemp, yTemp, zTemp

    integer,  allocatable :: iOuter(:)
    REAL(DP), allocatable :: aTemp(:), pTemp(:), cTemp(:)
    REAL(DP), allocatable, dimension(:,:,:) :: &
        x1(:,:,:), y1(:,:,:), z1(:,:,:), bx1(:,:,:), by1(:,:,:), bz1(:,:,:), &
        x2(:,:,:), y2(:,:,:), z2(:,:,:), bx2(:,:,:), by2(:,:,:), bz2(:,:,:), &
        x3(:,:,:), y3(:,:,:), z3(:,:,:), bx3(:,:,:), by3(:,:,:), bz3(:,:,:)

    ! Test Variables:
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'generate_field'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    nSWMF = 2*nPoints-1
    ALLOCATE(iOuter(nT))
    ALLOCATE(aTemp(0:nT+1), pTemp(nRExtend), cTemp(nSWMF))
    ALLOCATE(x1(nRExtend,nT,nthe), y1(nRExtend,nT,nthe), z1(nRExtend,nT,nthe))
    ALLOCATE(x2(npsi,nT,nthe),     y2(npsi,nT,nthe),     z2(npsi,nT,nthe))
    ALLOCATE(x3(npsi,nzeta,nthe),  y3(npsi,nzeta,nthe),  z3(npsi,nzeta,nthe))
    ALLOCATE(bx1(nRExtend,nT,nthe), by1(nRExtend,nT,nthe), bz1(nRExtend,nT,nthe))
    ALLOCATE(bx2(npsi,nT,nthe),     by2(npsi,nT,nthe),     bz2(npsi,nT,nthe))
    ALLOCATE(bx3(npsi,nzeta,nthe),  by3(npsi,nzeta,nthe),  bz3(npsi,nzeta,nthe))
    x1 = 0._dp; x2 = 0._dp; x3 = 0._dp; bx1 = 0._dp; bx2 = 0._dp; bx3 = 0._dp
    y1 = 0._dp; y2 = 0._dp; y3 = 0._dp; by1 = 0._dp; by2 = 0._dp; by3 = 0._dp
    z1 = 0._dp; z2 = 0._dp; z3 = 0._dp; bz1 = 0._dp; bz2 = 0._dp; bz3 = 0._dp


    BLines_DIII(4,:,:,:) = 10**9*BLines_DIII(4,:,:,:)
    BLines_DIII(5,:,:,:) = 10**9*BLines_DIII(5,:,:,:)
    BLines_DIII(6,:,:,:) = 10**9*BLines_DIII(6,:,:,:)

    if (DoTestMe) then
       FileName = RamFileName('SWMFField_0','dat',TimeRamNow)
       open(UNITTMP_, File=FileName)
       write(UNITTMP_,*) nRExtend, nT-1, nSWMF
       do i = 1,nSWMF
        do j = 1, nRExtend
         do k = 1,nT-1
          write(UNITTMP_,*) BLines_DIII(1,j,k,i),BLines_DIII(2,j,k,i),BLines_DIII(3,j,k,i), &
                            BLines_DIII(4,j,k,i),BLines_DIII(5,j,k,i),BLines_DIII(6,j,k,i)
         enddo
        enddo
       enddo
       close(UNITTMP_)
    endif

    ! Fix IsClosed_II so that we only take "well behaved" field lines
    iOuter = nRExtend + 1
    do j = 1, nT-1
       do i = nRExtend, 2, -1
          r0 = sqrt(maxval(BLines_DIII(1,i,j,:)**2 &
                         + BLines_DIII(2,i,j,:)**2 &
                         + BLines_DIII(3,i,j,:)**2))
          rr = sqrt(maxval(BLines_DIII(1,i-1,j,:)**2 &
                         + BLines_DIII(2,i-1,j,:)**2 &
                         + BLines_DIII(3,i-1,j,:)**2))
          ctemp = atan2(BLines_DIII(2,i,j,:),BLines_DIII(1,i,j,:))
          il = maxloc(ctemp)
          ii = il(1)
          il = minloc(ctemp)
          ir = il(1)
          t0 = atan2(BLines_DIII(1,i,j,ii)*BLines_DIII(2,i,j,ir) - BLines_DIII(1,i,j,ir)*BLines_DIII(2,i,j,ii), &
                     BLines_DIII(1,i,j,ii)*BLines_DIII(1,i,j,ir) + BLines_DIII(2,i,j,ir)*BLines_DIII(2,i,j,ii))
          if (maxval(abs(BLines_DIII(3,i,j,:))) > 5.0) then
             IsClosed_II(i,j) = .false.
          elseif (r0-rr > 1.0) then
             IsClosed_II(i,j) = .false.
          elseif (r0 > 9.0) then
             IsClosed_II(i,j) = .false.
          elseif (abs(t0) > 0.35) then
             IsClosed_II(i,j) = .false.
          endif
          if (.not.IsClosed_II(i,j)) then
             !write(*,*) 'Open Line', i, j, r0, rr, maxval(abs(BLines_DIII(3,i,j,:)))
             iOuter(j) = i
          endif
       enddo
    enddo
    iOuter = iOuter - 1

    nSWMF_to_nthe: do
       if (GSLerr > 0) exit nSWMF_to_nthe
       do j = 1, nT-1
          do i = 1, iOuter(j)
             distance(1) = 0._dp ! Distance along field line in Re
             DO k = 2, nSWMF
                distance(k) = distance(k-1) + SQRT((BLines_DIII(1,i,j,k)-BLines_DIII(1,i,j,k-1))**2 &
                                                  +(BLines_DIII(2,i,j,k)-BLines_DIII(2,i,j,k-1))**2 &
                                                  +(BLines_DIII(3,i,j,k)-BLines_DIII(3,i,j,k-1))**2)
             END DO
             distance(1:nSWMF) = distance(1:nSWMF) / distance(nSWMF)
             cTemp(1:nSWMF) = distance(1:nSWMF) * pi_d ! Distance along field line in radians (0-pi)
             ! Convert x,y,z points
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(1,i,j,1:nSWMF), &
                                       chiVal(1:nthe), x1(i,j,1:nthe), GSLErr)
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(2,i,j,1:nSWMF), &
                                       chiVal(1:nthe), y1(i,j,1:nthe), GSLErr)
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(3,i,j,1:nSWMF), &
                                       chiVal(1:nthe), z1(i,j,1:nthe), GSLErr)
             if (GSLerr > 0) exit nSWMF_to_nthe
             ! Convert bx,by,bz
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(4,i,j,1:nSWMF), &
                                       chiVal(1:nthe), bx1(i,j,1:nthe), GSLErr)
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(5,i,j,1:nSWMF), &
                                       chiVal(1:nthe), by1(i,j,1:nthe), GSLErr)
             call GSL_Interpolation_1D(ctemp(1:nSWMF), BLines_DIII(6,i,j,1:nSWMF), &
                                       chiVal(1:nthe), bz1(i,j,1:nthe), GSLErr)
             if (GSLerr > 0) exit nSWMF_to_nthe
          enddo
       enddo
       exit nSWMF_to_nthe
    enddo nSWMF_to_nthe

    if (DoTestMe) then
       FileName = RamFileName('SWMFField_1','dat',TimeRamNow)
       open(UNITTMP_, File=FileName)
       write(UNITTMP_,*) nRExtend, nT-1, nthe
       do i = 1, nthe
        do j = 1, nRExtend
         do k = 1, nT-1
          write(UNITTMP_,*) x1(j,k,i), y1(j,k,i), z1(j,k,i), &
                            bx1(j,k,i),by1(j,k,i),bz1(j,k,i)
         enddo
        enddo
       enddo
       close(UNITTMP_)
    endif

    nRExtend_to_npsi: do
       if (GSLerr > 0) exit nRExtend_to_npsi
       do j = 1, nT-1
          ii = iOuter(j)
          do k = 1, nthe
             distance(1) = 0._dp
             do i = 2, ii
                distance(i) = distance(i-1) + SQRT((x1(i,j,k)-x1(i-1,j,k))**2 &
                                                  +(y1(i,j,k)-y1(i-1,j,k))**2 &
                                                  +(z1(i,j,k)-z1(i-1,j,k))**2)
             enddo
             distance(1:ii) = distance(1:ii)/distance(ii)
             pTemp(1:ii) = distance(1:ii)
             ! Convert x,y,z points
             call GSL_Interpolation_1D(pTemp(1:ii),    x1(1:ii,j,k), &
                                       rhoVal(1:npsi), x2(1:npsi,j,k), GSLErr)
             call GSL_Interpolation_1D(pTemp(1:ii),    y1(1:ii,j,k), &
                                       rhoVal(1:npsi), y2(1:npsi,j,k), GSLErr)
             call GSL_Interpolation_1D(pTemp(1:ii),    z1(1:ii,j,k), &
                                       rhoVal(1:npsi), z2(1:npsi,j,k), GSLErr)
             if (GSLerr > 0) exit nRExtend_to_npsi
             ! Convert bx,by,bz
             call GSL_Interpolation_1D(pTemp(1:ii),    bx1(1:ii,j,k), &
                                       rhoVal(1:npsi), bx2(1:npsi,j,k), GSLErr)
             call GSL_Interpolation_1D(pTemp(1:ii),    by1(1:ii,j,k), &
                                       rhoVal(1:npsi), by2(1:npsi,j,k), GSLErr)
             call GSL_Interpolation_1D(pTemp(1:ii),    bz1(1:ii,j,k), &
                                       rhoVal(1:npsi), bz2(1:npsi,j,k), GSLErr)
             if (GSLerr > 0) exit nRExtend_to_npsi
          enddo
       enddo
       exit nRExtend_to_npsi
    enddo nRExtend_to_npsi

    if (DoTestMe) then
       FileName = RamFileName('SWMFField_2','dat',TimeRamNow)
       open(UNITTMP_, File=FileName)
       write(UNITTMP_,*) npsi, nT-1, nthe
       do i = 1, nthe
        do j = 1, npsi
         do k = 1, nT-1
          write(UNITTMP_,*) x2(j,k,i), y2(j,k,i), z2(j,k,i), &
                            bx2(j,k,i),by2(j,k,i),bz2(j,k,i)
         enddo
        enddo
       enddo
       close(UNITTMP_)
    endif

    nT_to_nzeta: do
       if (GSLerr > 0) exit nT_to_nzeta
       do j = 1,nT-1
          aTemp(j) = atan2(y2(1,j,nThetaEquator),x2(1,j,nThetaEquator))
       enddo
       where (aTemp(5:nT-3) < 0.0_dp) aTemp(5:nT-3) = aTemp(5:nT-3) + twopi_d
       aTemp(nT-2) = aTemp(nT-2) + twopi_d
       aTemp(nT-1) = aTemp(nT-1) + twopi_d
       aTemp(0) = aTemp(nT-1) - twopi_d
       aTemp(nT) = aTemp(1)   + twopi_d
       aTemp(nT+1) = aTemp(2) + twopi_d
       do i = 1, npsi
          do k = 1, nthe
             ! Convert x,y,z points
             xTemp(1:nT-1) = x2(i,1:nT-1,k)
             xTemp(0)      = xTemp(nT-1)
             xTemp(nT)     = xTemp(1)
             xTemp(nT+1)   = xTemp(2)
             yTemp(1:nT-1) = y2(i,1:nT-1,k)
             yTemp(0)      = yTemp(nT-1)
             yTemp(nT)     = yTemp(1)
             yTemp(nT+1)   = yTemp(2)
             zTemp(1:nT-1) = z2(i,1:nT-1,k)
             zTemp(0)      = zTemp(nT-1)
             zTemp(nT)     = zTemp(1)
             zTemp(nT+1)   = zTemp(2)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    xTemp(0:nT+1), &
                                       zetaVal(2:nzeta), x3(i,2:nzeta,k), GSLerr)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    yTemp(0:nT+1), &
                                       zetaVal(2:nzeta), y3(i,2:nzeta,k), GSLerr)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    zTemp(0:nT+1), &
                                       zetaVal(2:nzeta), z3(i,2:nzeta,k), GSLerr)
             if (GSLerr > 0) exit nT_to_nzeta
             ! Convert bx,by,bz
             xTemp(1:nT-1) = bx2(i,1:nT-1,k)
             xTemp(0)      = xTemp(nT-1)
             xTemp(nT)     = xTemp(1)
             xTemp(nT+1)   = xTemp(2)
             yTemp(1:nT-1) = by2(i,1:nT-1,k)
             yTemp(0)      = yTemp(nT-1)
             yTemp(nT)     = yTemp(1)
             yTemp(nT+1)   = yTemp(2)
             zTemp(1:nT-1) = bz2(i,1:nT-1,k)
             zTemp(0)      = zTemp(nT-1)
             zTemp(nT)     = zTemp(1)
             zTemp(nT+1)   = zTemp(2)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    xTemp(0:nT+1), &
                                       zetaVal(2:nzeta), bx3(i,2:nzeta,k), GSLerr)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    yTemp(0:nT+1), &
                                       zetaVal(2:nzeta), by3(i,2:nzeta,k), GSLerr)
             CALL GSL_Interpolation_1D(aTemp(0:nT+1),    zTemp(0:nT+1), &
                                       zetaVal(2:nzeta), bz3(i,2:nzeta,k), GSLerr)
             if (GSLerr > 0) exit nT_to_nzeta
          enddo
       enddo
       exit nT_to_nzeta
    enddo nT_to_nzeta
    x3(:,1,:)  = x3(:,nzeta,:)
    y3(:,1,:)  = y3(:,nzeta,:)
    z3(:,1,:)  = z3(:,nzeta,:)
    bx3(:,1,:) = bx3(:,nzeta,:)
    by3(:,1,:) = by3(:,nzeta,:)
    bz3(:,1,:) = bz3(:,nzeta,:)

    if (DoTestMe) then
       FileName = RamFileName('SWMFField_3','dat',TimeRamNow)
       open(UNITTMP_, File=FileName)
       write(UNITTMP_,*) npsi, nzeta-1, nthe
       do i = 1, nthe
        do j = 1, npsi
         do k = 1, nzeta-1
          write(UNITTMP_,*) x3(j,k,i), y3(j,k,i), z3(j,k,i), &
                            bx3(j,k,i),by3(j,k,i),bz3(j,k,i)
         enddo
        enddo
       enddo
       close(UNITTMP_)
    endif

    !!!!!!
    do k = 1, nthe
       x(k,1:npsi,1:nzeta) = x3(1:npsi,1:nzeta,k)
       y(k,1:npsi,1:nzeta) = y3(1:npsi,1:nzeta,k)
       z(k,1:npsi,1:nzeta) = z3(1:npsi,1:nzeta,k)
    enddo
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)

    do j = 1, nzeta+1
       do i = 1, npsi
          psi(:,i,j) = sqrt(x(1,i,j)**2+y(1,i,j)**2+z(1,i,j)**2) &
                       /(1._dp - z(1,i,j)**2/(x(1,i,j)**2+y(1,i,j)**2+z(1,i,j)**2))
       enddo
    enddo
    rc = maxval(psi(1,1,:))
    xpsiin = rc(1)
    rc = minval(psi(1,npsi,:))
    xpsiout = rc(1)
    xpsitot = xpsiout - xpsiin
    do i = 1, npsi
       psis = REAL(i-1,DP)/REAL(npsi-1,DP)
       xpl = xpsiin*(xpsiout/xpsiin)**(psis)
       psival(i) = -xzero3/xpl
    enddo
    psi = -xzero3/psi
    CALL mappsi
    CALL psiFunctions
    CALL maptheta
    !!!!!!

    if (DoTestMe) then
       FileName = RamFileName('SWMFField_4','dat',TimeRamNow)
       open(UNITTMP_, File=FileName)
       write(UNITTMP_,*) npsi, nzeta-1, nthe
       do i = 1, nthe
        do j = 1, npsi
         do k = 1, nzeta-1
          write(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k)
         enddo
        enddo
       enddo
       close(UNITTMP_)
    endif

    if (GSLerr > 0) SORFail = .true.
    DEALLOCATE(iOuter)
    DEALLOCATE(aTemp, pTemp, cTemp)
    DEALLOCATE(x1,  y1,  z1,  x2,  y2,  z2,  x3,  y3,  z3)
    DEALLOCATE(bx1, by1, bz1, bx2, by2, bz2, bx3, by3, bz3)

    return
  end subroutine generate_field
!==================================================================================================
END MODULE ModRamCouple
