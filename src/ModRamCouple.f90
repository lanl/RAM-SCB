!==============================================================================
module ModRamCouple
!    This module contains variables and routines for coupling RAM-SCB to
!    the SWMF: BATS-R-US and Ridley_Serial/RIM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModRamMain, ONLY: Real8_, PathRamOut
  use ModRamConst, ONLY: M1
  use ModRamGrids, ONLY: NR, NE, NT, NPA, NRExtend
  use ModRamTiming, ONLY: TimeRamElapsed, TimeRamNow
  use ModRamVariables, ONLY: KP, F107, EKEV, MU, PHI, LZ

  use ModRamParams

  implicit none
  save

  public :: set_type_mhd
  public :: generate_flux
  public :: sort_mhd_lines

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

  logical :: DoMultiFluidGMCoupling = .false.
  logical :: DoPassJr  = .false.
  logical :: DoIEPrecip = .false.

  ! Coupling P to SWMF:
  ! The iono footprints from BATS should be stored here:
  real(kind=Real8_), public, allocatable :: IonoMap_DSII(:,:,:,:)

  ! Container for MHD density and pressure at outer boundary.
  ! Indices are (dens:pres,Local Time,species[all, h, he, O])
  real(kind=Real8_), public, allocatable :: MhdDensPres_VII(:,:,:)
  real(kind=Real8_), public, allocatable :: FluxBats_IIS(:,:,:)
  real(kind=Real8_), public, allocatable :: PMhdGhost(:)! MHD pressure at RAM ghost cells.
  real(kind=Real8_), public, allocatable :: FluxBats_anis(:,:,:,:)

  ! Variables for E-field coupling with MHD (E=UxB):
  ! Only need equatorial values.
  real(kind=Real8_), public, allocatable :: ETotal_DII(:,:,:)
  
  ! Variables for B-field coupling with MHD:
  integer, public, parameter   :: nPointsMax = 200
  integer, public :: nRadSWMF, nRadSwmfVar, nRadSwmfVarInner, &
       nLonSWMF, nLinesSWMF, nPoints
  real(kind=Real8_), public, allocatable :: iEnd(:)
  real(kind=Real8_), public, allocatable :: MhdLines_IIV(:,:,:)
  real(kind=Real8_), public, allocatable :: xSWMF(:,:,:), ySWMF(:,:,:), &
       zSWMF(:,:,:), LatSWMF(:,:), LonSWMF(:,:)
  real(kind=Real8_), public, allocatable :: xEqSWMF(:,:), yEqSWMF(:,:), &
       pEqSWMF(:,:), nEqSWMF(:,:)
  real(kind=Real8_), public, allocatable :: uEqSWMF_DII(:,:,:), bEqSWMF_DII(:,:,:)

  ! Variables for coupling to SWMF-IE component:
  integer, public :: nIePhi=0, nIeTheta=0
  real(kind=Real8_), public, allocatable :: SwmfIonoPot_II(:,:)
  real(kind=Real8_), public, allocatable :: SwmfPot_II(:,:)

  real(kind=real8_), allocatable :: JrIono(:,:), energy_fluxIono(:,:), ave_eIono(:,:), &
                                    num_fluxIono(:,:), dis_energy_fluxIono(:,:), dis_ave_eIono(:,:)

contains

!==============================================================================
  subroutine RAMCouple_Allocate

    use ModRamGrids, ONLY: nT, nE, nRExtend, nR

    implicit none

    ALLOCATE(IonoMap_DSII(3,2,nRextend,nT), MhdDensPres_VII(3,nT,4), FluxBats_IIS(nE, nT, 1:4), &
             PMhdGhost(nT), FluxBats_anis(nE,nPa,nT,1:4), iEnd(2*(nRExtend)*nT), xEqSWMF(nRExtend,nT-1), &
             yEqSWMF(nRExtend,nT-1), pEqSWMF(nRExtend,nT-1), nEqSWMF(nRExtend,nT-1), SwmfPot_II(nR+1, nT), &
             uEqSWMF_DII(3,nRExtend,nT-1),bEqSWMF_DII(3,nRExtend,nT-1),ETotal_DII(2,nR,nT))

    SWMFPot_II = 0.
    FluxBats_anis = 0.
    PMhdGhost = 0.
    FluxBats_IIS = 0.
    IonoMap_DSII = 0.
    xEqSWMF = 0.
    yEqSWMF = 0.
    pEqSWMF = 0.
    nEqSWMF = 0.
    uEqSWMF_DII = 0.
    bEqSWMF_DII = 0.

  end subroutine RAMCouple_Allocate

!==============================================================================
  subroutine RAMCouple_Deallocate

    implicit none

    DEALLOCATE(IonoMap_DSII, MhdDensPres_VII, FluxBats_IIS, PMhdGhost, iEnd, &
               FluxBats_anis, xEqSWMF, yEqSWMF, pEqSWMF, nEqSWMF, SwmfPot_II)

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
       case('ux')
          Ux_=i
       case('uy')
          Uy_=i
       case('uz')
          Uz_=i  
       case('bx')
          Bx_=i
       case('by')
          By_=i
       case('bz')
          Bz_=i
       end select
       NameVar = trim( NameVar(nChar + 1:len(NameVar)) )
    end do

    if ( all((/ RhoH_,  RhoHe_, RhoO_ /)  .ne. -1) ) TypeMhd = 'multiS'
    if ( all((/ PresH_, PresHe_,PresO_ /) .ne. -1) ) TypeMhd = 'multiF'
    if ( all((/ PresSw_,PresH_ ,PresO_ /) .ne. -1) ) TypeMhd = 'mult3F'
    if ( Ppar_ .ne. -1) TypeMhd = 'anisoP'

    if(DoTestMe) then
       write(*,*) 'RAM_SCB: TypeMhd = ', TypeMhd
       write(*,*) 'RAM_SCB: MHD Indices = ', &
            'total = ', TotalRho_, TotalPres_, &
            'H+    = ', RhoH_, PresH_, &
            'He+   = ', RhoHe_,PresHe_, &
            'O+    = ', RhoO_, PresO_, &
            'Sw H+ = ', RhoSw_,PresSw_
    end if

  end subroutine set_type_mhd

  !===========================================================================
  subroutine sort_mhd_lines(nVarIn, nPointIn, BufferLine_VI)
    ! Sort MHD field lines & associated state variables into manageable array.
    ! Extend field lines in preparation for calculation of Euler potential
    ! surfaces.  Create equatorial variables.
    
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
    
    NewLine: do iLine=1, nLinesSWMF
       ! Skip missing lines (i.e., inside MHD boundary)
       if (iLine<BufferLine_VI(1, iPointBuff)) then
          iEnd(iLine) = -1
          cycle NewLine
       end if

       ! Grab all points on current line:
       do while(BufferLine_VI(1,iPointBuff)==iLine)
          MhdLines_IIV(iLine, iPointLine, :) = BufferLine_VI(:, iPointBuff)
          ! Do not exceed buffer:
          if (iPointBuff == nPointIn) exit
          ! Continue along line: 
          iPointBuff = iPointBuff + 1
          iPointLine = iPointLine + 1
       end do
       
       ! Check for very long lines or open lines.
       if (iPointLine>nPointsMax-5) then
          write(*,*) NameSub, 'ERROR: nPointsMax exceeded.'
          ! This line is too long.  Why?
          if (abs(MhdLines_IIV(iLine, iPointLine-1, 5)) > 5.0) write(*,*) &
               'Open field line?  Z=',MhdLines_IIV(iLine,iPointLine-1,5)
          call CON_stop(NameSub//'nPointsMax exceeded.')
       endif

       ! Save length of current line, reset.
       iEnd(iLine) = iPointLine - 1
       iPointLine = 1
    end do NewLine

    MhdLines_IIV(:,:,3:5) = MhdLines_IIV(:,:,3:5)/6378100.0 !XYZ in R_E.

    ! Now, extend lines to ionosphere.  Fill in missing lines w/ dipole lines.
    nPoints = maxval(iEnd)+5
    iLon=1; iRad=1
    do iLine=1, nLinesSWMF
       ! If line is missing...
       if(iEnd(iLine)<0)then
          xRam = lz(iRad)*cos(phi(iLon))
          yRam = lz(iRad)*sin(phi(iLon))
          MhdLines_IIV(iLine,1,3) = xRam
          MhdLines_IIV(iLine,1,4) = yRam
          MhdLines_IIV(iLine,1,5) = 0.0
          call get_dipole_trace((/xRam, yRam, 0.0_Real8_/), nPoints-1,&
               MhdLines_IIV(iLine,2:nPoints,3), &
               MhdLines_IIV(iLine,2:nPoints,4), &
               MhdLines_IIV(iLine,2:nPoints,5) )
          ! Southern hemisphere?  Use simple symmetry.
          if(mod(iLine, 2) == 0) MhdLines_IIV(iLine,:,5) = &
               -1.0*MhdLines_IIV(iLine,:,5)
          ! Fill equatorial values with simple conditions:
          ! Magnetic field: use dipole field, SI units are Tesla.
          MhdLines_IIV(iLine,1,Bx_:By_) = 0.0
          MhdLines_IIV(iLine,1,Bz_)     = 7.84E15 / lz(iRad)**3
          ! Flow velocity: corotation
          MhdLines_IIV(iLine,1,Ux_) = corot*lz(iRad)*sin(phi(iLon))
          MhdLines_IIV(iLine,1,Uy_) = corot*lz(iRad)*cos(phi(iLon))
          MhdLines_IIV(iLine,1,Uz_) = 0.0
       else
          i = iEnd(iLine)
          !write(*,*)'Extending!'
          !write(*,*)'Adding ', nPoints-i, 'points.'
          !write(*,'(a,3f9.6)') '...from ', x(iLine,i),y(iLine,i),z(iLine,i)
          call get_dipole_trace((/MhdLines_IIV(iLine,i,3), &
               MhdLines_IIV(iLine,i,4),MhdLines_IIV(iLine,i,5)/), nPoints-i, &
               MhdLines_IIV(iLine, i+1:nPoints, 3), &
               MhdLines_IIV(iLine, i+1:nPoints, 4), &
               MhdLines_IIV(iLine, i+1:nPoints, 5) )
          if(mod(iLine, 2) == 0) MhdLines_IIV(iLine,i+1:nPoints,5) &
               = -1.*MhdLines_IIV(iLine,i+1:nPoints,5)
          !write(*,'(a,3f9.6)') 'Got to:',x(iLine,nPoints), y(iLine, nPoints), &
          !z(iLine,nPoints)
          iEnd(iLine) = nPoints
       end if
       ! Keep track of our lon&rad.
       if(mod(iLine, 2) == 0 )iRad=iRad+1
       if(iRad > nRextend) then
          iLon = iLon+1
          iRad = 1
       end if
    end do

    !write(*,*) MhdLines_IIV(1,:nPoints,3)
    !write(*,*) MhdLines_IIV(nLinesSWMF, :nPoints, 3)

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
          iLine = 2*((nRextend+1)*(j-1) + i)-1 
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
    iRad=nRextend-2 ! Do not include ghost cells from coupling.
    ETotal_DII(1,:,1:nT-1) = uEqSWMF_DII(2,2:iRad,:)*bEqSWMF_DII(3,2:iRad,:) &
         -uEqSWMF_DII(3,2:iRad,:)*bEqSWMF_DII(2,2:iRad,:) ! X-component
    
    ETotal_DII(2,:,1:nT-1) = uEqSWMF_DII(1,2:iRad,:)*bEqSWMF_DII(3,2:iRad,:) &
         -uEqSWMF_DII(3,2:iRad,:)*bEqSWMF_DII(1,2:iRad,:) ! Y-component
    
  end subroutine sort_mhd_lines
  !===========================================================================
  subroutine generate_flux
    ! Convert MHD density and pressure into flux.
    ! Methodology depends on MHD: single, multispecies, or multifluid.

    use ModIoUnit, ONLY: UnitTmp_
    use ModConst,  ONLY: cProtonMass, cElectronCharge, cPi
    implicit none

    real, parameter :: cMass2Num = 1.0E6 * cProtonMass

    real(kind=Real8_) :: ratioHeH, ratioOH, fracH, fracHe, fracO
    real(kind=Real8_) :: factor, eCent, MhdPpar(4), MhdPper(4), dens, pres, ppar
    real(kind=Real8_) :: factor1, kappa, gamma1, gamma2, gamma3
    integer :: iT, iE, iPa
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
       ! Use Young et al. (JGR, 1982, Vol. 87 No. A11) findings to
       ! split mass densit into number density for H, He, O.
       ! Start by using empirical relationships from Young et al.
       ! to obtain number density ratios (O/H and He/H) for given
       ! Kp and F10.7 flux.
       ratioOH = 4.5E-2 * exp(0.17*Kp + 0.01*F107)
       ratioHeH= 0.618182*ratioOH*exp(-0.24*Kp - 0.011*F107) + 0.011*ratioOH
       ! Use number dens. ratios to calculate % of each species to total.
       fracH = 1.0 / (1.0 + ratioHeH + ratioOH)
       fracHe= ratioHeH * fracH
       fracO = ratioOH  * fracH
       
       if (DoTestMe) then
          write(*,*) 'Fractions used: fracH, fracHe, fracO ='
          write(*,*) fracH, fracHe, fracO
       end if

       do iT=1, nT
          ! Convert kg/m^-3 to total #/cm^-3 for given species content.
          MhdDensPres_VII(1,iT,1) = MhdDensPres_VII(1,iT,1) / &
               ( cMass2Num * (fracH + 4*fracHe + 16*fracO) )
          ! Convert P to temp in eV.
          MhdDensPres_VII(2,iT,1) = MhdDensPres_VII(2,iT,1) / &
               ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,iT,1) )
          ! Set other species number density.
          MhdDensPres_VII(1,iT,2) = fracH  * MhdDensPres_VII(1,iT,1)
          MhdDensPres_VII(1,iT,3) = fracHe * MhdDensPres_VII(1,iT,1)
          MhdDensPres_VII(1,iT,4) = fracO  * MhdDensPres_VII(1,iT,1)
          ! All other temps the same in single fluid.
          MhdDensPres_VII(2,iT,2:) = MhdDensPres_VII(2,iT,1)
       end do
       MhdDensPres_VII(2,:,1) = MhdDensPres_VII(2,:,1)/7.0

    case ('anisoP')

       ratioOH = 4.5E-2 * exp(0.17*Kp + 0.01*F107)
       ratioHeH= 0.618182*ratioOH*exp(-0.24*Kp - 0.011*F107) + 0.011*ratioOH
       ! Use number dens. ratios to calculate % of each species to total                                                  
       fracH = 1.0 / (1.0 + ratioHeH + ratioOH)
       fracHe= ratioHeH * fracH
       fracO = ratioOH  * fracH

       if (DoTestMe) then
          write(*,*) 'Fractions used: fracH, fracHe, fracO ='
          write(*,*) fracH, fracHe, fracO
       end if

       do iT=1, nT
          ! Convert kg/m^-3 to total #/cm^-3 for given species content.    
          MhdDensPres_VII(1,iT,1) = MhdDensPres_VII(1,iT,1) / &
               ( cMass2Num * (fracH + 4*fracHe + 16*fracO) )
          ! Convert P to temp in eV.                 
          MhdDensPres_VII(2,iT,1) = MhdDensPres_VII(2,iT,1) / &
               ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,iT,1) )
          ! Set other species number density.            
          MhdDensPres_VII(1,iT,2) = fracH  * MhdDensPres_VII(1,iT,1)
          MhdDensPres_VII(1,iT,3) = fracHe * MhdDensPres_VII(1,iT,1)
          MhdDensPres_VII(1,iT,4) = fracO  * MhdDensPres_VII(1,iT,1)
          ! All other temps the same in single fluid.   
          MhdDensPres_VII(2,iT,2:) = MhdDensPres_VII(2,iT,1)

          ! convert P_par to T_par in eV        
          MhdDensPres_VII(3,iT,1) = MhdDensPres_VII(3,iT,1)/ &
               ( cElectronCharge * 1.0e6 * MhdDensPres_VII(1,iT,1))
          ! all other temps the same           
          MhdDensPres_VII(3,iT,2:)=MhdDensPres_VII(3,iT,1)
       end do
       ! electron temperature is lower     
       MhdDensPres_VII(2,:,1) = MhdDensPres_VII(2,:,1)/7.0
       MhdDensPres_VII(3,:,1) = MhdDensPres_VII(3,:,1)/7.0

    case('multiS')
       do iT=1, nT
          ! Save total mass density to calculate temperature.
          ! Do this with # of protons because that's how single fluid
          ! (yet multi species) MHD is actually tracking the fluid.
          MhdDensPres_VII(1,iT,1) = sum(MhdDensPres_VII(1,iT,2:4))
          ! Convert total and species mass density to number density.
          MhdDensPres_VII(1,iT,2) = MhdDensPres_VII(1,iT,2)!/ 1.0
          MhdDensPres_VII(1,iT,3) = MhdDensPres_VII(1,iT,3) / 4.0
          MhdDensPres_VII(1,iT,4) = MhdDensPres_VII(1,iT,4) /16.0
          ! Convert P to temp in eV.
          MhdDensPres_VII(2,iT,1) = MhdDensPres_VII(2,iT,1) / &
               ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,iT,1) )
          ! All other temps the same in single fluid.
          MhdDensPres_VII(2,iT,2:) = MhdDensPres_VII(2,iT,1)
       end do

    case('multiF')
       ! Convert total and species mass density to number density.
       MhdDensPres_VII(1,:,2) = MhdDensPres_VII(1,:,2)/(cMass2Num)
       MhdDensPres_VII(1,:,3) = MhdDensPres_VII(1,:,3)/(cMass2Num*4.0)
       MhdDensPres_VII(1,:,4) = MhdDensPres_VII(1,:,4)/(cMass2Num*16.0)
       MhdDensPres_VII(1,:,1) = sum(MhdDensPres_VII(1,:,2:4))
       
       ! Convert P to temp in eV.
       MhdDensPres_VII(2,:,1) = MhdDensPres_VII(2,:,1) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,1) )
       MhdDensPres_VII(2,:,2) = MhdDensPres_VII(2,:,2) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,2) )
       MhdDensPres_VII(2,:,3) = MhdDensPres_VII(2,:,3) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,3) )
       MhdDensPres_VII(2,:,4) = MhdDensPres_VII(2,:,4) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,4) )

    case('mult3F')
       ! Convert total and species mass density to number density (#/cm3)
       ! Helium slot is now solar wind protons.
       MhdDensPres_VII(1,:,2) = MhdDensPres_VII(1,:,2)/(cMass2Num)
       MhdDensPres_VII(1,:,3) = MhdDensPres_VII(1,:,3)/(cMass2Num)
       MhdDensPres_VII(1,:,4) = MhdDensPres_VII(1,:,4)/(cMass2Num*16.0)
       MhdDensPres_VII(1,:,1) = sum(MhdDensPres_VII(1,:,2:4))
       
       ! Convert Pa to temp in eV.
       MhdDensPres_VII(2,:,1) = MhdDensPres_VII(2,:,1) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,1) )
       MhdDensPres_VII(2,:,2) = MhdDensPres_VII(2,:,2) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,2) )
       MhdDensPres_VII(2,:,3) = MhdDensPres_VII(2,:,3) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,3) )
       MhdDensPres_VII(2,:,4) = MhdDensPres_VII(2,:,4) / &
            ( cElectronCharge * 1.0E6 * MhdDensPres_VII(1,:,4) )

    case default
       call CON_stop(NameSub//' TypeMhd not recognized: '//TypeMhd)
    end select

    if(mod(TimeRamElapsed, 60.0_Real8_) .eq. 0.0) then
       ! Write boundary dens and pressure to file.
       write(NameFile,'(a,i6.6,a)')&
            PathRamOut//"bound_plasma_t",nint(TimeRamElapsed/60._Real8_),".out"
       open( UnitTmp_, FILE=NameFile, STATUS='replace')
       write(UnitTmp_,'(es13.5,i5,5i3)') &
            TimeRamElapsed, &
            TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
            TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond
       write(UnitTmp_, *) 'Units are Hours, cm-3, and eV'

       if (TypeMhd .eq. 'mult3F') then
          write(UnitTmp_, *) 'lt RhoE RhoH RhoSw RhoO pE pH pSw pO'
       elseif(TypeMhd .eq. 'anisoP')then
          write(UnitTmp_, *) 'lt RhoE RhoH RhoHe RhoO pE pH pHe pO pparE pparH pparHe pparO'
       else
          write(UnitTmp_, *) 'lt RhoE RhoH RhoHe RhoO pE pH pHe pO'
       endif

       do iT=1, nT
          if (TypeMhd .ne. 'anisoP')then
             write(UnitTmp_, '(i2.2, 8(1x, E13.6))') &
                  iT, MhdDensPres_VII(1,iT,1:4), MhdDensPres_VII(2,iT,1:4)
         else
             write(UnitTmp_, '(i2.2, 12(1x, E13.6))') &
                  iT,MhdDensPres_VII(1,iT,1:4),MhdDensPres_VII(2,iT,1:4), &
                  MhdDensPres_VII(3,iT,1:4)
          end if
       end do
       close(UnitTmp_)
    end if

    kappa = 3.
    gamma1 = 6.0    !gamma(kappa+1)  
    gamma2 = 1.3293 !gamma(kappa-0.5)  
    gamma3 = 0.886  ! gamma(1.5)
    factor1 = gamma1/gamma2/gamma3*sqrt(2*cPi)/(2*kappa-3)**1.5

    do iE=1, nE
       eCent = 1000.0 * Ekev(iE) ! Center of energy bin in eV.
       factor=4.0E6*Ekev(iE)     ! Unit conversion factor * eCent in KeV.
       do iT=1, nT
          if (TypeMhd .ne. 'anisoP')then
             if (NameDistrib .eq. 'MAXW') then

                ! Assuming a maxwellian distribution about MHD temperature, divide
                ! densities into energy bins and calculate fluxes.

                FluxBats_IIS(iE, iT, 1) = 1./sqrt(M1(1))*exp(-1.0*eCent/MhdDensPres_VII(2,iT,1))*&
                     factor * MhdDensPres_VII(1,iT,1) * &
                     (MhdDensPres_VII(2,iT,1)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE, iT, 2) = 1./sqrt(M1(2))*exp(-1.0*eCent/MhdDensPres_VII(2,iT,2))*&
                     factor * MhdDensPres_VII(1,iT,2) * &
                     (MhdDensPres_VII(2,iT,2)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE, iT, 3) = 1./sqrt(M1(3))*exp(-1.0*eCent/MhdDensPres_VII(2,iT,3))*&
                     factor * MhdDensPres_VII(1,iT,3) * &
                     (MhdDensPres_VII(2,iT,3)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE, iT, 4) = 1./sqrt(M1(4))*exp(-1.0*eCent/MhdDensPres_VII(2,iT,4))*&
                     factor * MhdDensPres_VII(1,iT,4) * &
                     (MhdDensPres_VII(2,iT,4)/1000.0_Real8_)**(-1.5)
             elseif(NameDistrib .eq.'KAPA')then
                ! assuming a kappa distribution
                
                FluxBats_IIS(iE,iT,1) = 1./sqrt(M1(1))*(1+1.0*eCent/(&
                     (kappa-1.5)*MhdDensPres_VII(2,iT,1)))**(-(kappa+1))*&
                     factor1*factor * MhdDensPres_VII(1,iT,1) * &
                     (MhdDensPres_VII(2,iT,1)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE,iT,2) = 1./sqrt(M1(2))*(1+1.0*eCent/(&
                     (kappa-1.5)*MhdDensPres_VII(2,iT,2)))**(-(kappa+1))*&
                     factor1*factor * MhdDensPres_VII(1,iT,2) * &
                     (MhdDensPres_VII(2,iT,2)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE,iT,3) = 1./sqrt(M1(3))*(1+1.0*eCent/(&
                     (kappa-1.5)*MhdDensPres_VII(2,iT,3)))**(-(kappa+1))*&
                     factor1*factor * MhdDensPres_VII(1,iT,3) * &
                     (MhdDensPres_VII(2,iT,3)/1000.0_Real8_)**(-1.5)
                FluxBats_IIS(iE,iT,4) = 1./sqrt(M1(4))*(1+1.0*eCent/(&
                     (kappa-1.5)*MhdDensPres_VII(2,iT,4)))**(-(kappa+1))*&
                     factor1*factor * MhdDensPres_VII(1,iT,4) * &
                     (MhdDensPres_VII(2,iT,4)/1000.0_Real8_)**(-1.5)
             else
                call CON_stop(NameSub//' NameDistribution not recognized: '//NameDistrib)
             end if
          else  ! anisotropic Mhd
             ! Tpar and Tper in eV; 
             MhdPpar = MhdDensPres_VII(3,iT,1:4)
             MhdPper = (MhdDensPres_VII(2,iT,1:4)*3-MhdPpar)/2.0
             if (NameDistrib .eq. 'MAXW')then
                ! assume a bi-maxwellian distribution                              
                do iPa=1, nPa
                   FluxBats_anis(iE,iPa,iT,1) = 1./sqrt(M1(1)) &
                        * exp(-1.0*eCent*MU(ipa)**2/MhdPpar(1) - 1.0*eCent*(1-mu(ipa)**2)/MhdPper(1)) &
                        * factor * MhdDensPres_VII(1,iT,1) &
                        /(MhdPper(1)/1000.0_Real8_*sqrt(MhdPpar(1)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,2) = 1./sqrt(M1(2)) &
                        * exp(-1.0*eCent*MU(ipa)**2/MhdPpar(2) - 1.0*eCent*(1-MU(ipa)**2)/MhdPper(2)) &
                        * factor * MhdDensPres_VII(1,iT,2) &
                        /(MhdPper(2)/1000.0_Real8_*sqrt(MhdPpar(2)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,3) = 1./sqrt(M1(3)) &
                        * exp(-1.0*eCent*MU(ipa)**2/MhdPpar(3) - 1.0*eCent*(1-MU(ipa)**2)/MhdPper(3)) &
                        * factor * MhdDensPres_VII(1,iT,3) &
                        /(MhdPper(3)/1000.0_Real8_*sqrt(MhdPpar(3)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,4) = 1./sqrt(M1(4)) &
                        * exp(-1.0*eCent*MU(ipa)**2/MhdPpar(4) - 1.0*eCent*(1-MU(ipa)**2)/MhdPper(4)) &
                        * factor * MhdDensPres_VII(1,iT,4) &
                        /(MhdPper(4)/1000.0_Real8_*sqrt(MhdPpar(4)/1000.0_Real8_))
                end do
             elseif(NameDistrib .eq. 'KAPA')then
                ! assume a bi-kappa distribution 
                do iPa=1, nPa
                   FluxBats_anis(iE,iPa,iT,1) = 1./sqrt(M1(1)) * &
                        (1+1.0*eCent*MU(ipa)/((kappa-1.5)*MhdPpar(1) ) + &
                        1.0*eCent*(1-MU(ipa)**2)/((kappa-1.5)*MhdPper(1)))**(-(kappa+1))*&
                        factor1*factor * MhdDensPres_VII(1,iT,1)  &
                        /(MhdPper(1)/1000.0_Real8_*sqrt(MhdPpar(1)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,2) = 1./sqrt(M1(2)) * &
                        (1+1.0*eCent*MU(ipa)/((kappa-1.5)*MhdPpar(2) ) + &
                        1.0*eCent*(1-MU(ipa)**2)/((kappa-1.5)*MhdPper(2)))**(-(kappa+1))*&
                        factor1*factor * MhdDensPres_VII(1,iT,2)  &
                        /(MhdPper(2)/1000.0_Real8_*sqrt(MhdPpar(2)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,3) = 1./sqrt(M1(3)) * &
                        (1+1.0*eCent*MU(ipa)/((kappa-1.5)*MhdPpar(3) ) + &
                        1.0*eCent*(1-MU(ipa)**2)/((kappa-1.5)*MhdPper(3)))**(-(kappa+1))*&
                        factor1*factor * MhdDensPres_VII(1,iT,3)  &
                        /(MhdPper(3)/1000.0_Real8_*sqrt(MhdPpar(3)/1000.0_Real8_))
                   FluxBats_anis(iE,iPa,iT,4) = 1./sqrt(M1(4)) * &
                        (1+1.0*eCent*MU(ipa)/((kappa-1.5)*MhdPpar(4) ) + &
                        1.0*eCent*(1-MU(ipa)**2)/((kappa-1.5)*MhdPper(4)))**(-(kappa+1))*&
                        factor1*factor * MhdDensPres_VII(1,iT,4)  &
                        /(MhdPper(4)/1000.0_Real8_*sqrt(MhdPpar(4)/1000.0_Real8_))
                end do
             else
                call CON_stop(NameSub//' NameDistribution not recognized: '//NameDistrib)
             end if
          end if

       end do
    end do
    
    ! Remove ridiculously small values.
    if (TypeMhd .ne. 'anisoP')then
       where(FluxBats_IIS .lt. 1.0E-30) FluxBats_IIS = 0.0
    else
       where(FluxBats_anis .lt. 1.0E-30) FluxBats_anis = 0.0
    end if

    ! If using two hydrogen species (case of Mult3F), one is saved in 
    ! helium's spot.  Combine H+ fluxes and set He flux to small fraction.
    if (TypeMhd .eq. 'mult3F') then
       FluxBats_IIS(:,:,2) = FluxBats_IIS(:,:,2) + FluxBats_IIS(:,:,3)
       FluxBats_IIS(:,:,3) = 0.01 * FluxBats_IIS(:,:,2) ! 1% of H+ flux.
    end if

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

    character(len=*), parameter :: NameSub = "write_FluxGM"
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
          if(TypeMhd .ne. 'anisoP')then
             write(UnitTmp_,StringFormat) 1.0, iT-1.0, FluxBats_IIS(:,iT,iS)
          else
             do iPa=1,nPa
                write(UnitTmp_,StringFormat) 1.0, iT-1.0, FluxBats_anis(:,ipa,iT,iS)
             end do
          end if
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
    
    use ModConst,   ONLY: cElectronCharge, cPi
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

  !===========================================================================
  subroutine calculate_precip_flux_jr(IS, nTheta, nPhi, rIono, &
       energy_flux, ave_e, num_flux, dis_energy_flux, dis_ave_e,&
       Jr,high_latboundary, DoPrecipFlux, DoJr)

!    use ModRamMain, ONLY: WMU, mu, UPA, PI, EKEV, XNE, IsRestart, &
!         WE, RadiusMax, PathRamOut, TimeRamElapsed, Dt_hI, LZ, FLUX, &
!         TimeRamNow, PathRestartIn, XNE,PPERE,PPARE
    use ModIoUnit,   ONLY: UnitTmp_
!    use Module1, ONLY: nXRaw, nYRaw, radGrid, angleGrid, x, y, z,&
!         r0Start, nthe, npsi, nzeta, paraj, nThetaEquator,&
!         bf,ppar,pper,pnormal, bj, nzetaMidnight
    use ModConst,    ONLY: cElectronMass, cElectronCharge
    use ModNumConst, ONLY: cPi, cHalfPi
    use ModKind

    use ModRamConst,     ONLY: PI
    use ModRamMain,      ONLY: PathRamOut, PathRestartIn
    use ModRamParams,    ONLY: IsRestart
    use ModRamGrids,     ONLY: nR, nT, nE, nPA, RadiusMax
    use ModScbGrids,     ONLY: nthe, npsi, nzeta, nXRaw, nYRaw
    use ModRamVariables, ONLY: WMU, MU, UPA, EKEV, XNE, WE, LZ, FLUX, PPerE, PParE
    use ModScbVariables, ONLY: x, y, z, r0Start, paraj, bf, ppar, pper, pnormal, &
                               bj, nThetaEquator, nZetaMidnight, radGrid, angleGrid
    use ModRamTiming,    ONLY: TimeRamElapsed, Dt_hI, TimeRamNow

    use ModScbSpline,    ONLY: Spline_2D_Point

    implicit none

    character(len=2) :: species
    integer, intent(in) :: IS, nTheta,nPhi
    logical, intent(in) :: DoPrecipFlux, DoJr
    real(real8_), intent(in) :: rIono
    real(real8_),dimension(nTheta,nPhi), intent(out):: energy_flux, &
         ave_e, num_flux, Jr, dis_energy_flux, dis_ave_e
    real(real8_),intent(out) :: high_latboundary(nPhi)
    real(real8_) :: diff_flux_iono(nTheta, nPhi, nE)
    REAL(real8_) :: dydummy,radius, angle, beta=0.1,Nele,jpar
    integer:: nTheta_north
    integer:: i, j, k,kk, l, iDomain,ierr,ntotal, nzetap, idx(nPhi),d
    REAL, allocatable :: colat(:), lon(:)
    real(real8_), allocatable :: xscatter(:),yscatter(:), &
         energy_fluxscatter(:), num_fluxscatter(:), ave_escatter(:), &
         parajscatter(:), diff_fluxscatter(:,:),&
         dis_ave_escatter(:),dis_energy_fluxscatter(:)
    real(real8_):: ave_flux(nR,nT,nE),f(nR,nT,NE,NPA),ave_fluxext(nR+1,nT,nE),num_fluxeq(nR,nT)
    real(real8_):: coordPoint(2),NeExt(nR+1,nT),PperEExt(nR+1,nT),PparEExt(nR+1,nT)
    real(real8_):: radRaw(nXRaw+2), azimRaw(nYRaw), Lpp(nzeta), Lpp_colat(nzeta), Lpp_lon(nzeta)
    real:: rr1,  thangle, thangleOnIono, minl, dTheta, dPhi,Rm,eV,efactor,press
    real, dimension(1:npsi,1:nzeta):: colatGrid, lonGrid,F0
    real(real8_), dimension(1:npsi,1:nzeta):: energy_flux_iono,num_flux_iono,ave_e_iono, &
         XNE_EQ, NeEQ, PperEEQ, PparEEQ, dis_num_flux_iono, dis_ave_e_iono, dis_energy_flux_iono
    real(real8_), dimension(1:npsi,1:nzeta,nE):: ave_fluxEQ
    real(real8_), allocatable :: ave_fluxtmp(:,:), ave_fluxEQtmp(:,:), energy_fluxtmp(:,:), &
         ave_etmp(:,:), num_fluxtmp(:,:), jr_tmp(:,:),dis_energy_fluxtmp(:,:), &
         dis_ave_etmp(:,:), dis_num_fluxtmp(:,:), high_latboundary_tmp(:)
    character(len=100) :: NameFile
    logical            :: DoTest
    !---------------------------------------------------------------------------
    !\
    ! by Yiqun Yu @ 2015
    ! Calculate energy flux and average energy from the precipitation number 
    ! flux at the ionosphere altitude.
    !/
    ! -- 1:  map the precipitation flux at the equator down to the ionosphere 
    !     altitude. Interpolate into the ionospheric grids.
    ! -- 2:  calculate the energy flux and ave. energy.
    !            <E> = \int{f(E)EdE}/\int{f(E)dE}
    !            jE  = \int{f(E)EdE}
    !            Here, f(E) is the energy differential flux, already integrated 
    !                   over the pitch angle.
    ! This is used to pass to IE for the computation of height-integrated Hall 
    ! and Pedeson Conductances, based on the Robinson 1987 formula. 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !\ Y. Yu, 2013, to couple with IE (Region2 Jr)
    ! interpolate the paraj onto the ionosphere grids
    ! like in the swmfiono_II <-- PhiIono(1:npsi,2:zeta) <-- PhiIono(nIeTheta,
    ! nIePhi)
    !/
    ! paraJr is positive when along the field line into the north pole (this is 
    ! opposite to that in GM)
    !----------------------------------------------------------------------------


    energy_flux = 0.0
    ave_e = 0.0
    num_flux = 0.0
    Jr = 0.0
    high_latboundary = 0.5*cPi

    if (.not. DoPrecipFlux .and. .not. DoJr) return

    if (TimeRamElapsed .lt. Dt_hI .or. maxval(x) .eq. 0) then
       ! before the IM component is computed
       if (.not. IsRestart)then
          energy_flux = 0.0
          ave_e = 0.0
          num_flux = 0.0
          Jr = 0.0
          high_latboundary = 0.5*cPi

       else
          if (DoPrecipFlux)then

             NameFile=PathRestartIn//'/restart_precip1.rst'
             open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
             read(UnitTMP_) energy_flux
             close(UnitTMP_)

             NameFile=PathRestartIn//'/restart_precip2.rst'
             open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
             read(UnitTMP_) ave_e
             close(UnitTMP_)

             NameFile=PathRestartIn//'/restart_precip3.rst'
             open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
             read(UnitTMP_) num_flux
             close(UnitTMP_)
          end if
          if (DoPassJr)then
             NameFile=PathRestartIn//'/restart_jr.rst'
             open(unit=UnitTMP_,file=trim(NameFile),status='old',form='unformatted')
             read(UnitTMP_) Jr
             close(UnitTMP_)
          end if

       end if

       return
    end if

    if (.not. DoPrecipFlux .and. .not. DoJr) return
    if (DoPrecipFlux)then
       ! FLUX is saved from ram_all for each species before F2 is evolved.
       F(2:NR,1:nT-1, 2:NE, 2:NPA) = FLUX(IS, 2:NR, 1:NT-1, 2:NE, 2:NPA)
       F(2:NR, nT, 2:NE, 2:NPA) = F(2:NR, 1, 2:NE, 2:NPA)
       F(2:NR, :, 2:NE, 1) = F(2:NR, :, 2:NE, 2)

       ! calculate the averaged flux at the equator for the precipitation flux
       ! (Jordanova 1997).
       ! averaged over the pitch angle (in the loss cone), so to remove the
       ! angle dependence.
       ! then for the low altitude flux, it has no pitch angle dependence (i.e.,
       ! isotropic now).
       do i=2, nR
          do j=1, nT
             num_fluxeq(i,j) = 0.0
             do k=2,nE
                ave_flux(i,j,k) = 0.0
                do l=upa(i), npa
                   ave_flux(i,j,k) = ave_flux(i,j,k)+f(i,j,k,l)*WMU(l)
                end do
                ave_flux(i,j,k) = ave_flux(i,j,k)/(mu(NPA)-mu(UPA(i)))
                num_fluxeq(i,j) = num_fluxeq(i,j) + PI*ave_flux(i,j,k)*WE(k)
             end do
          end do
       end do

       ! interpolate into the scb grid on the equator
       DO j = 1, nXRaw+1
          radRaw(j) = LZ(j) !1.75_dp + (radOut-1.75_dp) * real(J-1)/real(nXRaw+1)
       END DO
       DO k = 1, nYRaw ! starting from midnight (MLT=0)
          azimRaw(k) = 24.* REAL(k-1)/REAL(nYRaw-1)
       END DO

       if (IS .eq. 1) species = 'e_'
       if (IS .eq. 2) species = 'h_'
       if (IS .eq. 3) species = 'he'
       if (IS .eq. 4) species = 'o_'

       !extend a bit to 6.75Re
       do j = nXRaw+1, nXRaw+2
          radRaw(j) = radRaw(nXRaw) + REAL(j-nXRaw)*(radRaw(nXRaw)-radRaw(1))/REAL(nXRaw-1)
       end do

       ave_fluxExt(1:nXRaw+1,:,:) = ave_flux(1:nXRaw+1,:,:)
       do k = 2, nE
          do j = 1, nYRaw
             i = nXRaw+2
             ! polinomial interpolation/extrpolation
             call polint(radRaw(nXRaw-1:nXRaw+1), ave_fluxExt(nXRaw-1:nXRaw+1,j,k), radRaw(i), ave_fluxExt(i,j,k), dydummy)
          end do
       end do
       where(ave_fluxExt .le. 1.0e-31)ave_fluxExt = 1.0e-31

       ! interpolate the Ne [1/cm^3]
       NeExt(1:nXRaw+1,:) = XNE(1:nXRaw+1,:)
       PperEExt(1:nXRaw+1,:) = PPERE(1:nXRaw+1,:)
       PparEExt(1:nXRaw+1,:) = PPARE(1:nXRaw+1,:)
       do j=1,nYRaw
          i = nXRaw+2
          call polint(radRaw(nXRaw-1:nXRaw+1), NeExt(nXRaw-1:nXRaw+1,j), radRaw(i), NeExt(i,j), dydummy)
          call polint(radRaw(nXRaw-1:nXRaw+1), PperEExt(nXRaw-1:nXRaw+1,j), radRaw(i), PperEExt(i,j), dydummy)
          call polint(radRaw(nXRaw-1:nXRaw+1), PparEExt(nXRaw-1:nXRaw+1,j), radRaw(i), PparEExt(i,j), dydummy)
       end do

       !\
       ! debug: write out the flux into the loss cone at the equator
       !/
       DoTest = .true.
       if(DoTest)then
          if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
             write(NameFile,'(a,a,a,i6.6,a)')&
                  PathRamOut//"PrecipFlux_atEquator_",species,"t",nint(TimeRamElapsed/300._Real8_),".dat"
             open( UnitTmp_, FILE=NameFile, STATUS='replace')
             write(UnitTmp_,*)'Time: ',TimeRamElapsed
             write(UnitTmp_,*)'nR, nMLT, nE', NR,NT,nE-1
             write(UnitTmp_,*)'L    MLT  Energy  Flux[/cm^2/s/keV] '
             do i=2, nR+1
                do j=1, nT
                   do k=2, nE
                      write(UnitTmp_,'(f8.4,1x,f8.4,1x,f8.4,1x,E13.6)')radRaw(i), azimRaw(j), EKEV(k), ave_fluxExt(i,j,k)
                   end do
                end do
             end do
             close(UnitTmp_)
          end if

          if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
             write(NameFile,'(a,a,a,i6.6,a)')&
                  PathRamOut//"PrecipNumberFlux_atEquator_",species,"t",nint(TimeRamElapsed/300._Real8_),".dat"
             open( UnitTmp_, FILE=NameFile, STATUS='replace')
             write(UnitTmp_,*)'Time: ',TimeRamElapsed,' nR, nMLT ', NR-1,NT
             write(UnitTmp_,*)'L    MLT  Flux[/cm^2/s] '
             do i=2, nR
                do j=1, nT
                   write(UnitTmp_,'(f8.4,1x,f8.4,1x,E13.6)')radRaw(i), azimRaw(j), num_fluxeq(i,j)
                end do
             end do
             close(UnitTmp_)
          end if
       end if

       DO k = 2, nzeta
          DO j = 1, npsi
             radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
             angle = ASIN(y(nThetaEquator,j,k) / radius) + cPi
             IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
                  angle = 2*cPi - ASIN(y(nThetaEquator,j,k) / radius)
             IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
                  angle = - ASIN(y(nThetaEquator,j,k) / radius)
             radGrid(j,k) = radius
             angleGrid(j,k) = angle ! 0 degree at midnight, index starting from noon.
          END DO
       END DO
       ! interpolate into scb grids in the equator (ram index start from MLT=0)
       ! angleGrid: index start from noon. 2--> 3.14; ...nzeta/2-->6.28,
       ! nzeta/2+1-->0...; 
       ! radGrid(npsi,:) is bit outside 6.75, so not well extroplated
       if (.not. allocated(ave_fluxtmp)) allocate(ave_fluxtmp(1:nXRaw+2,1:nYRaw))
       if (.not. allocated(ave_fluxEQtmp)) allocate(ave_fluxEQtmp(1:npsi,2:nzeta))

       ! interpolate into scb grids in equator
       call Spline_2D_point(radRaw(1:nXRaw+2)**2, cPi/12.*azimRaw, NeExt(1:nXRaw+2, 1:nYRaw), &
            radGrid(1:npsi, 2:nzeta)**2, angleGrid(1:npsi, 2:nzeta), NeEQ(1:npsi, 2:nzeta), iDomain)
       where(NeEQ .le.0.0)NeEQ = 0.0
       where(radGrid(1:npsi,2:nzeta) .gt. radRaw(nXRaw+2))NeEQ(1:npsi,2:nzeta) = 0.0
       call Spline_2D_point(radRaw(1:nXRaw+2)**2, cPi/12.*azimRaw, PperEExt(1:nXRaw+2, 1:nYRaw), &
            radGrid(1:npsi, 2:nzeta)**2, angleGrid(1:npsi, 2:nzeta), PperEEQ(1:npsi, 2:nzeta), iDomain)
       where(PperEEQ .le.0.0)PperEEQ = 0.0
       where(radGrid(1:npsi,2:nzeta) .gt. radRaw(nXRaw+2))PperEEQ(1:npsi,2:nzeta) = 0.0
       call Spline_2D_point(radRaw(1:nXRaw+2)**2, cPi/12.*azimRaw, PparEExt(1:nXRaw+2, 1:nYRaw), &
            radGrid(1:npsi, 2:nzeta)**2, angleGrid(1:npsi, 2:nzeta), PparEEQ(1:npsi, 2:nzeta), iDomain)
       where(PparEEQ .le.0.0)PparEEQ = 0.0
       where(radGrid(1:npsi,2:nzeta) .gt. radRaw(nXRaw+2))PparEEQ(1:npsi,2:nzeta) = 0.0

       do k=1, nE

          ave_fluxtmp(1:nXRaw+2, 1:nYRaw) = ave_fluxExt(1:nXRaw+2,1:nYRaw,k)

          call Spline_2D_point(radRaw(1:nXRaw+2)**2, cPi/12.*azimRaw, ave_fluxtmp(1:nXRaw+2, 1:nYRaw), &
               radGrid(1:npsi, 2:nzeta)**2, angleGrid(1:npsi, 2:nzeta), ave_fluxEQtmp(1:npsi, 2:nzeta), iDomain)

          ! for those negative, set to zeros
          where(ave_fluxEQtmp .le. 0.0)ave_fluxEQtmp = 0.0
          ! for those outside 6.75Re, set to zeros
          where(radGrid(1:npsi,2:nzeta) .gt. radRaw(nXRaw+2))ave_fluxEQtmp(1:npsi,2:nzeta) = 0.0

          ave_fluxEQ(1:npsi,2:nzeta,k) = ave_fluxEQtmp(1:npsi,2:nzeta)

       end do


       dis_num_flux_iono = 0.0
       dis_energy_flux_iono = 0.0
       dis_ave_e_iono = 0.0
       ! calculate the discreate number flux and energy flux at the equator in
       ! the scb grids
       do i=1,npsi
          do j=2,nzeta
             ! (electron) ppar in normalized unit, so *pnormal to nPa
             ! NeEQ in 1/cm^3
             ! bj in muA/m^2
             ! at the low altitude boundary
             press = (2*PperEEQ(i,j)+PparEEQ(i,j))/3.0*0.16*1.0e-9 !keV/cm^3-->nPa--> Pa
             Nele = NeEQ(i,j)*1.0e6 !1/m^3 at equator...

             ! bj(1,i,j) at the pole; bj(ntheequator,i,j) at the equator (!=0?)
             ! !units of paraj is muA/m^2
             jpar = abs(paraj(i,j))*1.0e-6 ! A/m^2 (paraj at the boundary)

!             ! following Zhang 2015
!             F0(i,j) = beta*sqrt(Nele*press)/sqrt(2*PI*cElectronMass) !/m^2/s
!             at equator
!             ! flux at the ionosphere is equal to at equator if E is constant
!             Rm = bf(nthe,i,j)/bf(nThetaEquator,i,j) ! Biono
!             if(jpar/cElectronCharge/F0(i,j) .ge.1 .and.
!             jpar/cElectronCharge/F0(i,j) .le. Rm)then
!                eV =
!                press/Nele*(Rm-1)*log((Rm-1)/(Rm-jpar/cElectronCharge/F0(i,j)))
!                efactor = exp(-1.0*eV/(press/Nele)*(Rm-1))
!                dis_num_flux_iono(i,j) = jpar/cElectronCharge *1.0e-4
!                !/m^2/s-->/cm^2/s
!                dis_energy_flux_iono(i,j) =
!                dis_num_flux_iono(i,j)*(2*press/Nele + &
!                     eV*(1-efactor)/(1+(1-1.0/Rm)*efactor))/1.6e-16
!                     !J/cm^2/s--> keV/cm^2/s
!                dis_ave_e_iono(i,j) =
!                dis_energy_flux_iono(i,j)/dis_num_flux_iono(i,j) !keV
!                write(*,*)'efactor,eV,Rm,F0:',efactor,eV,Rm,F0(i,j),
!                jpar/cElectronCharge/F0(i,j)
!                write(*,*)'press,jpar,Nele:',press,jpar,Nele,dis_num_flux_iono(i,j),dis_energy_flux_iono(i,j),dis_ave_e_iono(i,j)
!                write(*,*)'Te(keV), 1-efactor,
!                1+(1-1/Rm)*efactor:',press/Nele/1.6e-16, 1-efactor,
!                1+(1-1./Rm)*efactor
!             end if

             ! following Raeder 2001 (use parameter at ionosphere altitude
             if (paraj(i,j) .lt. 0 .and. press .gt. 0)then ! upward FACs region in the northern
                eV = cElectronCharge**2*Nele**1.5/sqrt(2*PI*cElectronMass*press)*jpar
                dis_energy_flux_iono(i,j) = 4*eV*jpar/1.6e-12         ! J/m^2/s--> keV/cm^2/s
                dis_ave_e_iono(i,j) = cElectronCharge*eV/1.6e-16      ! J-->keV
!                write(*,*)'paraj(i,j), press, Nele,eV,
!                dis_energy_flux_iono:',paraj(i,j), press, Nele, eV,
!                dis_energy_flux_iono(i,j),dis_ave_e_iono(i,j)
             end if

          end do
       end do

       !\
       ! Calculate the energy flux and averaged energy at 200km (along the B):
       !  -- precip_flux(at 200km) at all directions is equal to ave_flux (at
       !  equator).
       !  -- (the ave_flux is the same for all the pitch angles; 
       !  -- precip_flux(90deg) = ave_flux(at edge of the loss cone)
       !  -- isotropic now across the 200km atm. plane (same at r0Start, along B
       !  line)
       !  -- convert to the flux perpendicular to thep plane: flux*cos(alpha). 
       !  -- integrate over solid anle (--> *pi) and energy (--> energy flux)
       !/

       do i=1, npsi
          do j=2, nzeta
             energy_flux_iono(i,j) = 0.0
             num_flux_iono(i,j) = 0.0
             ave_e_iono(i,j) = 0.0
             do k=1, nE
                energy_flux_iono(i,j) = energy_flux_iono(i,j) + PI*ave_fluxEQ(i,j,k)*EKEV(k)*WE(k)
                num_flux_iono(i,j) = num_flux_iono(i,j) + PI*ave_fluxEQ(i,j,k)*WE(k)
             end do
             if (num_flux_iono(i,j) .eq. 0.0) then
                num_flux_iono(i,j) = 1.0e-31
                ave_e_iono(i,j) = 1.0e-31
             else
                ave_e_iono(i,j) = energy_flux_iono(i,j)/num_flux_iono(i,j)
             end if
          end do
       end do

!       if (DoTest)then
          if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
             write(NameFile,'(a,a,a,i6.6,a)')&
                  PathRamOut//"PrecipFlux_atSCBEquator_",species,"t",nint(TimeRamElapsed/300._Real8_),".dat"
             open( UnitTmp_, FILE=NameFile, STATUS='replace')
             write(UnitTmp_,*)'Time: ',TimeRamElapsed
             write(UnitTmp_,*)'nR, nMLT, nE', npsi,nzeta,nE
             write(UnitTmp_,*)'rad    angle  Energy  Flux[/cm^2/s/keV] '
             do i=1, npsi
                do j=2, nzeta
                   do k=1, nE
                      write(UnitTmp_,'(f8.4,1x,f8.4,1x,f8.4,1x,E14.6)')radGrid(i,j), angleGrid(i,j), EKEV(k), ave_fluxEQ(i,j,k)
                   end do
                end do
             end do
             close(UnitTmp_)
          end if

          if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
             write(NameFile,'(a,a,a,i6.6,a)')&
                  PathRamOut//"PrecipNumberFlux_atSCBEquator_",species,"t",nint(TimeRamElapsed/300._Real8_),".dat"
             open( UnitTmp_, FILE=NameFile, STATUS='replace')
             write(UnitTmp_,*)'Time: ',TimeRamElapsed
             write(UnitTmp_,*)'nR, nMLT', npsi,nzeta
             write(UnitTmp_,*)'rad    angle    Flux[/cm^2/s]   Energy_Flux [ergs/cm^2] disFlux[/cm^2/s] disEnergyFlux[ergs/cm^2]'
             do i=1, npsi
                do j=2, nzeta
                   write(UnitTmp_,'(f8.4,1x,f8.4,4(1x,E13.6))')radGrid(i,j), angleGrid(i,j), &
                        num_flux_iono(i,j), energy_flux_iono(i,j)*1.6e-9, &
                        dis_num_flux_iono(i,j), dis_energy_flux_iono(i,j)*1.6e-9
                end do
             end do

             close(UnitTmp_)
          end if

       end if

       if (allocated(ave_fluxtmp)) deallocate(ave_fluxtmp)
       if (allocated(ave_fluxEQtmp)) deallocate(ave_fluxEQtmp)
!    end if

    ! Now this is mapped along the B field lines to the ionosphere altitude
    ! (~200km as the loss cone is calculated there)
    ! assume around the Earth surface (r0Start): the first grid point in 3D
    ! equli. code?
    i = nthe ! the northest point, r0Start, not at the ionosphere altitude. So need to map to IonoAltitude 
    DO k = 2, nzeta
       DO j = 1, npsi
          rr1 = sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2)
          thangle = asin(z(nthe,j,k)/rr1) ! latitude at the northest point
          thangleOnIono = acos(sqrt( (cos(thangle))**2 * rIono/r0Start )) ! latitude at IonoAltitude
          colatGrid(j,k) = 0.5*cPi - thangleOnIono
       END DO
    END DO
    ! the longitude is the same as in angleGrid (index start from the noon, 0
    ! degree at midnight)
    lonGrid(:,2:nzeta) = angleGrid(:,2:nzeta)


    if(DoTest)then
       !\
       ! find the plasmapause location in scb equator, then the colat in the
       ! ionosphere
       ! interpolate the plasmasphere density into scb grids first
       !/
       call Spline_2D_point(radRaw(1:nXRaw+1)**2, cPi/12.*azimRaw, XNE(1:nXRaw+1, 1:nYRaw), &
            radGrid(1:npsi, 2:nzeta)**2, angleGrid(1:npsi, 2:nzeta), XNE_EQ(1:npsi, 2:nzeta), iDomain)
       do k=2, nzeta
          do i=1, npsi
             if (XNE_EQ(i,k) .ge. 50)then
                Lpp(k) = radGrid(i,k)
                Lpp_colat(k) = colatGrid(i,k)
                Lpp_lon(k) = angleGrid(i,k)
             endif
          end do
       end do
       !\
       ! write out the footprint of plasmapause on the ionosphere
       !/
       if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
          write(NameFile,'(a,i6.6,a)')&
               PathRamOut//"Precip_Lpp_t",nint(TimeRamElapsed/300._Real8_),".dat"
!          NameFile =
!          trim(PathRamout)//RamFileName('Precip_Lpp_','dat',TimeRamNow)
          open( UnitTmp_, FILE=NameFile, STATUS='replace')
          write(UnitTmp_,*)'Time: ',TimeRamElapsed
          write(UnitTmp_,*)'nZeta', nzeta-1
          write(UnitTmp_,*)'Lpp    Lpp_colatitude'
          do k=2, nzeta
             write(UnitTmp_,'(f8.4,1x,f8.4)')Lpp(k), Lpp_colat(k)
          end do
          close(UnitTmp_)
       end if
    end if

    ! now interpolate the scb spatial grid into the ionospheric grids.
    ! (nR, nT) --> (colatgrid, longrid)--> (colat, lon)
    dPhi = 2.0*cPi/real(nPhi-1)
    dTheta = cPi/real(nTheta-1)

    nTheta_north = int(nTheta/2)+1
    ! only find the north hemsiphere grid for interpolation
    if(.not. allocated(colat) .and. .not. allocated(lon)) then
       allocate(colat(nTheta_north), lon(nPhi))
       colat = 0.0
       lon   = 0.0
       do i=2, nPhi ! Longitude goes from 0 to 360, index start from midnight, lon=0 is at midnight
          lon(i) = lon(i-1) + dPhi
       end do
       do i=2, nTheta_north ! Colat goes from 0 to 90.
          colat(i) = colat(i-1) + dTheta ! all the hemisphere
       end do
    end if

   nTotal = npsi*(nzeta-1)+1
   if (.not.allocated(xScatter)) allocate(xscatter(ntotal),stat=ierr)
   if (.not.allocated(yScatter)) allocate(yscatter(ntotal),stat=ierr)
   if (.not.allocated(energy_fluxScatter)) allocate(energy_fluxScatter(ntotal),stat=ierr)
   if (.not.allocated(ave_eScatter)) allocate(ave_eScatter(ntotal),stat=ierr)
   if (.not.allocated(num_fluxScatter)) allocate(num_fluxScatter(ntotal),stat=ierr)
   if (.not.allocated(diff_fluxScatter)) allocate(diff_fluxScatter(nE,ntotal),stat=ierr)
   if (.not.allocated(dis_energy_fluxScatter))allocate(dis_energy_fluxScatter(ntotal),stat=ierr)
   if (.not.allocated(dis_ave_eScatter)) allocate(dis_ave_eScatter(ntotal),stat=ierr)
   if (.not.allocated(parajScatter)) allocate(parajScatter(ntotal),stat=ierr)

   dis_energy_fluxscatter = 0.0
   dis_ave_escatter = 0.0
   k = 1
   do i=1,npsi ! don't include the high boundary grid (often discontinuity) radial direction
      do j=2, nzeta ! azimuthal direction (for nzeta+1: x(or y) is equal to at "1")
         xScatter(k) = rIono*sin(colatGrid(i,j)) * cos(lonGrid(i,j))
         yScatter(k) = rIono*sin(colatGrid(i,j)) * sin(lonGrid(i,j))
         if (DoPrecipFlux) then
            energy_fluxScatter(k) = energy_flux_iono(i,j)
            ave_eScatter(k) = ave_e_iono(i,j)
            num_fluxScatter(k) = num_flux_iono(i,j)

            if(paraj(i,j).lt.0)then ! where Jpar is upward (opposite to Bnorm in northern hemisphere)
               dis_energy_fluxScatter(k) = dis_energy_flux_iono(i,j)
               dis_ave_eScatter(k) = dis_ave_e_iono(i,j)
            end if
            do l=1,nE
               diff_fluxScatter(l,k) = ave_fluxEQ(i,j,l)
            end do
         end if
         if (DoJr) parajScatter(k) = paraj(i,j) ! paraj is the FACs at northern pole (muA/m^2)
         k = k + 1
      end do
   end do

  ! Add a point in the center for proper triangulation
   xScatter(k) = 0.
   yScatter(k) = 0.
   energy_fluxScatter(k) = energy_flux_iono(npsi, nzetaMidnight)
   ave_eScatter(k) = ave_e_iono(npsi, nzetaMidnight)
   num_fluxScatter(k) = num_flux_iono(npsi,nzetaMidnight)
   dis_energy_fluxScatter(k) = dis_energy_flux_iono(npsi, nzetaMidnight)
   dis_ave_eScatter(k) = dis_ave_e_iono(nPsi,nzetaMidnight)
   do l=1, nE
      diff_fluxScatter(l,k) = ave_fluxEQ(npsi,nzetaMidnight,l)
   end do
   if (DoJr)parajScatter(k) = paraj(npsi, nzetaMidnight)

   do j=1,nPhi
      ! find the closest longitude
      idx(j) = 2
      minl = abs(lon(j) - lonGrid(1,2))
      kloop: do k=2,nzeta
         if ( abs(lon(j) - lonGrid(1,k)) .lt. minl)then
            minl = abs(lon(j) - lonGrid(1,k))
            idx(j) = k
         end if
      end do kloop
   end do

   !\
   ! Ionosphere precipitation interpolation onto iono. grids
   !/
   if (DoPrecipFlux) then
      call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), energy_fluxScatter(1:nTotal))
      do j=1, nPhi
         do i=1, nTheta_north
            if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
            ! inside the polar inner or outer boundary of scb
               ! don't do the interpolation
               energy_flux(i,j) = 0.0
            else
               coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
               coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
               CALL nnpntd(coordPoint(1), coordPoint(2), energy_flux(i,j))
            end if
            if (colat(i) .lt. minval(colatGrid(:,idx(j))))k = i
         end do
         ! for high-latitude boundary, gradual fading instead of a
         ! discontinuity...
         !do d=1,3
         !   energy_flux(k+1-d,j) = energy_flux(k+1,j)*exp(-d/1.5)
         !end do
      end do
      call nnpntendd()
      call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), ave_eScatter(1:nTotal))
      do j=1, nPhi
         do i=1, nTheta_north
            if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
               ! inside the polar inner or outer boundary of scb
               ! don't do the interpolation
               ave_e(i,j) = 0.0
            else
               coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
               coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
               CALL nnpntd(coordPoint(1), coordPoint(2), ave_e(i,j))
            end if
            if (colat(i) .lt. minval(colatGrid(:,idx(j))))k = i
         end do
         ! for high-latitude boundary, gradual fading instead of a
         ! discontinuity...
!         do d=1,3
!            ave_e(k+1-d,j) = ave_e(k+1,j)*exp(-d/1.5)
!         end do
      end do
      call nnpntendd()
      call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), num_fluxScatter(1:nTotal))

      do j=1, nPhi
         do i=1, nTheta_north
            if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
               num_flux(i,j) = 0.0
            else
               coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
               coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
               CALL nnpntd(coordPoint(1), coordPoint(2), num_flux(i,j))
            end if
            if (colat(i) .lt. minval(colatGrid(:,idx(j))))k = i
         end do
         ! for high-latitude boundary, gradual fading instead of a
         ! discontinuity...
!         DO d=1,3
!            num_flux(k+1-d,j)   = num_flux(k+1,j)*exp(-d/1.5)
!         END DO
      end do
      call nnpntendd()

      call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), dis_energy_fluxScatter(1:nTotal))
      do j=1, nPhi
         do i=1, nTheta_north
            if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
               dis_energy_flux(i,j) = 0.0
            else
               coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
               coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
               CALL nnpntd(coordPoint(1), coordPoint(2), dis_energy_flux(i,j))
            end if
            if (colat(i) .lt. minval(colatGrid(:,idx(j))))k = i
         end do
         ! for high-latitude boundary, exponentially fading instead of a
         ! discontinuity...
!         DO d=1,3
!            dis_energy_flux(k+1-d,j)   = dis_energy_flux(k+1,j)*exp(-d/1.5)
!         END DO
      end do
      call nnpntendd()

      call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), dis_ave_eScatter(1:nTotal))
      do j=1, nPhi
         do i=1, nTheta_north
            if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
               dis_ave_e(i,j) = 0.0
            else
               coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
               coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
               CALL nnpntd(coordPoint(1), coordPoint(2), dis_ave_e(i,j))
            end if
            if (colat(i) .lt. minval(colatGrid(:,idx(j))))k = i
         end do
         ! for high-latitude boundary, exponentially fading instead of a
         ! discontinuity...
!         DO d=1,3
!            dis_ave_e(k+1-d,j)   = dis_ave_e(k+1,j)*exp(-d/1.5)
!         END DO
      end do
      call nnpntendd()
      do kk=1, nE
         call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), diff_fluxScatter(k,1:nTotal))

         do j=1, nPhi
            do i=1, nTheta_north
               if(colat(i) .lt. minval(colatGrid(:,idx(j))) .or. colat(i) .gt. maxval(colatGrid(:,idx(j))))then
                  diff_flux_iono(i,j,kk) = 0.0
               else
                  coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
                  coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
                  CALL nnpntd(coordPoint(1), coordPoint(2), diff_flux_iono(i,j, kk))
               end if
               if (colat(i) .lt. minval(colatGrid(:,idx(j)))) k = i
            end do
            ! for high-latitude boundary, exponentially fading instead of a
            ! discontinuity...
!            DO d=1,3
!               diff_flux_iono(k+1-d,j,kk)   =
!               diff_flux_iono(k+1,j,kk)*exp(-d/1.5)
!            END DO
         end do
         call nnpntendd()
      end do

      if (allocated(energy_fluxscatter)) deallocate(energy_fluxscatter, stat=ierr)
      if (allocated(ave_escatter)) deallocate(ave_escatter, stat=ierr)
      if (allocated(num_fluxscatter)) deallocate(num_fluxscatter, stat=ierr)
      if (allocated(dis_energy_fluxscatter)) deallocate(dis_energy_fluxscatter, stat=ierr)
      if (allocated(dis_ave_escatter)) deallocate(dis_ave_escatter, stat=ierr)
      if (allocated(diff_fluxscatter)) deallocate(diff_fluxscatter, stat=ierr)

       ! map to the south hemisphere
       do i=1, nTheta_north
          energy_flux(nTheta_north+i-1,:)      = energy_flux(nTheta_north-i+1,:)
          ave_e(nTheta_north+i-1,:)            = ave_e(nTheta_north-i+1,:)
          num_flux(nTheta_north+i-1,:)         = num_flux(nTheta_north-i+1,:)
          dis_energy_flux(nTheta_north+i-1,:)  = dis_energy_flux(nTheta_north-i+1,:)
          dis_ave_e(nTheta_north+i-1,:)        = dis_ave_e(nTheta_north-i+1,:)
          diff_flux_iono(nTheta_north+i-1,:,:) = diff_flux_iono(nTheta_north-i+1,:,:)
       end do

       !\
       ! write out the precipitation on the ionosphere
       !/
       if (Mod(int(TimeRamElapsed) ,300) .eq. 0)then
          write(NameFile,'(a,a,a,i6.6,a)')&
               PathRamOut//"PrecipFlux_",species,"t",nint(TimeRamElapsed/300._Real8_),".dat"
!          NameFile =
!          trim(PathRamout)//RamFileName('PrecipFlux_e','dat',TimeRamNow)
          open( UnitTmp_, FILE=NameFile, STATUS='replace')
          write(UnitTmp_,*)'Time: ',TimeRamElapsed
          write(UnitTmp_,*)'nTheta_north, nPhi', nTheta_north, nPhi
          write(UnitTmp_,*)'Theta    Phi EnergyFlux[ergs/cm^2/s] Ave_eIono [keV] Num_Flux [/cm^2/s](Phi=0 at midnight)',&
               'discreteEnergyFlux[ergs/cm^2/s] discreteAveE[keV]'
          do i=1, nTheta_north
             do j=1, nPhi
                write(UnitTmp_,'(f8.4,1x,f8.4,6(1x,E13.6))')colat(i), lon(j), &
                     Energy_flux(i,j)*1.6e-9, Ave_e(i,j), num_flux(i,j), &
                     dis_energy_flux(i,j)*1.6e-9, dis_ave_e(i,j)
             end do
          end do
          close(UnitTmp_)

          ! write the differential electron flux down to the ionosphere.
!          write(NameFile,'(a,a,a,i6.6,a)')&
!               PathRamOut//"PrecipFlux_diff_",species,"t",nint(TimeRamElapsed/600._Real8_),".dat"
!          NameFile =
!          trim(PathRamout)//RamFileName('PrecipFlux_diff_e','dat',TimeRamNow)
!          open( UnitTmp_, FILE=NameFile, STATUS='replace')
!          write(UnitTmp_,*)'Time: ',TimeRamElapsed
!          write(UnitTmp_,*)'nTheta_north, nPhi', nTheta_north, nPhi
!          write(UnitTmp_,*)'Theta    Phi  Energy[keV]   Num_diff_Flux
!          [/cm^2/s/keV]'
!          do i=1, nTheta_north
!             do j=1, nPhi
!                do k=1, nE
!                   write(UnitTmp_,'(f8.4,1x,f8.4,1x, f8.4,1x,E13.6)')colat(i),
!                   lon(j), EKEV(k), diff_flux_iono(i,j,k)
!                end do
!             end do
!          end do
!          close(UnitTmp_)

       end if

       ! in order to couple to the ionosphere module that starts the index at
       ! noon with phi=0, shift the array about 180 degree
       allocate(Energy_fluxtmp(nTheta,nPhi), Ave_etmp(nTheta,nPhi), &
                Num_fluxtmp(nTheta,nPhi), dis_Energy_fluxtmp(nTheta,nPhi), &
                dis_Ave_etmp(nTheta,nPhi))
       Energy_fluxtmp(:, 1:nPhi/2+1) = Energy_flux(:, 1:nPhi/2+1)
       Energy_flux(:,1:nPhi/2+1)     = Energy_flux(:,nPhi/2+1:nPhi)
       Energy_flux(:,nPhi/2+1:nPhi)  = Energy_fluxtmp(:,1:nPhi/2+1)
       Ave_etmp(:, 1:nPhi/2+1)         = Ave_e(:, 1:nPhi/2+1)
       Ave_e(:,1:nPhi/2+1)           = Ave_e(:,nPhi/2+1:nPhi)
       Ave_e(:,nPhi/2+1:nPhi)        = Ave_etmp(:,1:nPhi/2+1)
       Num_fluxtmp(:, 1:nPhi/2+1)      = Num_flux(:, 1:nPhi/2+1)
       Num_flux(:,1:nPhi/2+1)        = Num_flux(:,nPhi/2+1:nPhi)
       Num_flux(:,nPhi/2+1:nPhi)     = Num_fluxtmp(:,1:nPhi/2+1)

       dis_Energy_fluxtmp(:, 1:nPhi/2+1) = dis_Energy_flux(:, 1:nPhi/2+1)
       dis_Energy_flux(:,1:nPhi/2+1)     = dis_Energy_flux(:,nPhi/2+1:nPhi)
       dis_Energy_flux(:,nPhi/2+1:nPhi)  = dis_Energy_fluxtmp(:,1:nPhi/2+1)
       dis_Ave_etmp(:, 1:nPhi/2+1)         = dis_Ave_e(:, 1:nPhi/2+1)
       dis_Ave_e(:,1:nPhi/2+1)           = dis_Ave_e(:,nPhi/2+1:nPhi)
       dis_Ave_e(:,nPhi/2+1:nPhi)        = dis_Ave_etmp(:,1:nPhi/2+1)
       deallocate(Energy_fluxtmp, Ave_etmp, Num_fluxtmp, dis_energy_fluxtmp, dis_ave_etmp)
    end if

    !\
    ! FAC interpolation onto the iono. grids
    !/
    if (DoJr) then
       call NNPNTINITD(nTotal, xScatter(1:nTotal), yScatter(1:nTotal), parajScatter(1:nTotal))
       do j=1,nPhi
          ! initialize the high latitude boundary
          high_latboundary(j) = 0.5*cPi
          do i=1, nTheta_north
             if (colat(i).lt.minval(colatGrid(:,idx(j))).or. &
                 colat(i).gt.maxval(colatGrid(:,idx(j)))) then
                ! inside the polar inner or outer boundary of scb
                ! don't do the interpolation
                Jr(i,j) = 0.0
                ! save the outer/high-lat boundary for the ionosphere
                if (colat(i) .lt. minval(colatGrid(:,idx(j))))then
                   if (high_latboundary(j) .gt. (0.5*cPi - colat(i))) then
                      high_latboundary(j) = 0.5*cPi-colat(i)
                   end if
                end if
             else
                coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
                coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
                CALL nnpntd(coordPoint(1), coordPoint(2), Jr(i,j))
             end if
          end do
       end do
       if (allocated(parajscatter)) deallocate(parajscatter, stat=ierr)
       ! Jr is positive into the north pole, opposite to that in GM;
       ! to be consistent: -1*Jr
       Jr(1:nTheta_north,:) = -Jr(1:nTheta_north,:)

       ! map to the south hemisphere
       do i=1, nTheta_north
          Jr(nTheta_north+i-1,:) = Jr(nTheta_north-i+1,:)
       end do

       ! writing out the mapped Jr on the ionosphere
!       if (mod(int(TimeRamElapsed),300) .eq. 0)then
!          
!          write(NameFile,'(a,i6.6,a)')&
!               PathRamOut//"IonoJr_t",nint(TimeRamElapsed/300._Real8_),".dat"
!          open( UnitTmp_, FILE=NameFile, STATUS='replace')
!          write(UnitTmp_,*)'Time: ', TimeRamElapsed
!          write(UnitTmp_,*)'nTheta_north, nPhi', nTheta_north, nPhi
!          write(UnitTmp_,*)'Theta    Phi     Jr [A/m^2]'
!          do i=1, nTheta_north
!             do j=1, nPhi
!                write(UnitTmp_,'(f8.4,1x,f8.4,1x,E13.6)')colat(i), lon(j),
!                Jr(i,j)*1.0e-6
!                
!             end do
!          end do
!          close(UnitTmp_)
!       end if

       ! rotate to the ionosphere index
       if(.not.allocated(jr_tmp)) allocate(jr_tmp(nTheta,nPhi))
       Jr_tmp(:,1:nPhi/2+1) = Jr(:,1:nPhi/2+1)
       Jr(:,1:nPhi/2+1)     = Jr(:,nPhi/2+1:nPhi)
       Jr(:,nPhi/2+1:nPhi)  = Jr_tmp(:,1:nPhi/2+1)
       if(allocated(jr_tmp)) deallocate(jr_tmp)

       if (.not.allocated(high_latboundary_tmp)) allocate(high_latboundary_tmp(1:nPhi/2+1))
       high_latboundary_tmp(1:nPhi/2+1) = high_latboundary(1:nPhi/2+1)
       high_latboundary(1:nPhi/2+1) = high_latboundary(nPhi/2+1:nPhi)
       high_latboundary(nPhi/2+1:nPhi) = high_latboundary_tmp(1:nPhi/2+1)
       if(allocated(high_latboundary_tmp))deallocate(high_latboundary_tmp)
    end if

    if (allocated(xscatter)) deallocate(xscatter, stat=ierr)
    if (allocated(yscatter)) deallocate(yscatter, stat=ierr)

  end subroutine calculate_precip_flux_jr

  !===========================================================================



end module ModRamCouple
!==============================================================================
