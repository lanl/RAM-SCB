! Wrapper for Inner Magnetosphere (IM) component
module IM_wrapper
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!******************************************************************************
  
  ! Wrapper for RAM_SCB when used in SWMF/IM component mode.

  implicit none

  private ! except

  public:: IM_set_param
  public:: IM_init_session
  public:: IM_run
  public:: IM_save_restart
  public:: IM_finalize

  ! Coupling with IE
  public:: IM_get_for_ie
  public:: IM_put_from_ie_mpi
  public:: IM_put_from_ie
  public:: IM_put_from_ie_complete

  ! Coupling with GM
  public:: IM_get_for_gm
  public:: IM_put_from_gm
  public:: IM_put_from_gm_line
  public:: IM_put_from_gm_crcm
  public:: IM_put_sat_from_gm

  contains
  !=============================================================================
  subroutine IM_set_param(CompInfo, TypeAction)
    
    use CON_comp_info
    use ModReadParam
    use ModRamMpi
    use ModRamIO
    use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case
    
    implicit none
    
    character (len=*), parameter :: NameSub='IM_set_param'
    
    ! Arguments
    type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
    character (len=*), intent(in)     :: TypeAction ! What to do
    !--------------------------------------------------------------------------
    
    select case(TypeAction)
       
    case('VERSION')
       call put(CompInfo,                         &
            Use=.true.,                           &
            NameVersion='RAM-SCB (Jordanova, Zaharia)', &
            Version=2.0)
       
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       IsFramework = .true.
       
    case('CHECK')
       if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
       RETURN
       !call ram_check ! No sub: ram_check yet.
       
    case('GRID')
       call IM_set_grid
       
    case('READ')
       call IM_set_parameters
       
       ! Items with no current actions:
    case('STDOUT')
    case('FILEOUT')
       
    case default
       call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)
       
    end select
    
  end subroutine IM_set_param

  !============================================================================
  subroutine IM_set_grid
    use ModNumConst, ONLY: cTwoPi
    use ModRamMain,  ONLY: RadiusMin, RadiusMax, NR, NT, PI, LZ, PHI, DL1, &
         DPHI, dR, nRextend, GridExtend
    use CON_coupler, ONLY: set_grid_descriptor, is_proc, IM_
    implicit none
    
    character (len=*), parameter :: NameSub='IM_set_grid'
    logical :: IsInitialized=.false.
    logical :: DoTest, DoTestMe
    real    :: RadMaxFull
    integer :: i, j, iS
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)'IM_set_grid called, IsInitialized=', &
         IsInitialized
    if(IsInitialized) return

    IsInitialized = .true.
    
    RadMaxFull = RadiusMin + nRextend*dR
    
    ! IM grid: the equatorial grid is described by Coord1_I and Coord2_I
    ! Occasional +0.0 is used to convert from single to double precision
    ! Grid information set here is shared with GM for coupling with BATS-R-US.
    call set_grid_descriptor( IM_,           & ! component index
         nDim     = 2,                       & ! dimensionality
         nRootBlock_D = (/1,1/),             & ! number of blocks
         nCell_D =(/nRextend, nT/),          & ! size of equatorial grid
         XyzMin_D=(/RadiusMin+0.0,0.0/),     & ! min coordinates
         XyzMax_D=(/RadMaxFull,cTwoPi/),     & ! max coordinates + radial ghost.
         Coord1_I = GridExtend+0.0,          & ! radial coordinates
         Coord2_I = Phi(1:nT)+0.0,           & ! longitudinal coordinates
         TypeCoord= 'SMG',                   & ! solar magnetic coord
         nVar=7,                             & ! Number of "fluid" vars.
         NameVar = 'rho p ppar Hpp Opp Hprho Oprho')
    ! For nVar/NameVar, RAM-SCB+BATS-R-US coupling does not need this
    ! convention.  Placeholders are used ONLY.
    
    if(DoTestMe)then
       write(*,*)NameSub,' NR =              ', NR
       write(*,*)NameSub,' NR + Ghostcells = ', nRextend
       write(*,*)NameSub,' NT =              ', NT
    end if

  end subroutine IM_set_grid

  !============================================================================
  subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)
    
    ! Provide current for IE
    ! The value should be interpolated from nPoints with
    ! indexes stored in Index and weights stored in Weight
    ! The variables should be put into Buff_V
    
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    implicit none
    character(len=*), parameter :: NameSub='IM_get_for_ie'
    
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    !--------------------------------------------------------------------------
    
  end subroutine IM_get_for_ie

  !============================================================================
  subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)
    
    use ModRamMain,  ONLY: Real8_, TimeRamElapsed, PathRamOut
    use ModRamCouple,ONLY: SwmfIonoPot_II, nIePhi, nIeTheta
    use ModPlotFile, ONLY: save_plot_file
    
    implicit none
    
    integer, intent(in):: nTheta, nPhi
    real,    intent(in):: Potential_II(nTheta, nPhi,1)
    
    real(kind=Real8_), parameter :: ratioGrid = 7.0/18.0
    integer :: iTheta10, iTheta80
    character(len=100):: NameFile
    
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'IM_put_from_ie_mpi'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    ! If testing this function, write the potential to file.
    if(DoTest)then
       write(NameFile,'(a,i5.5,a)') &
            PathRamOut//"potential_t",nint(TimeRamElapsed),".out"
       call save_plot_file(NameFile, &
            StringHeaderIn = 'Ionospheric potential', &
            TimeIn         = TimeRamElapsed+0.0, &
            NameVarIn      = 'Theta Phi Pot', &
            CoordMinIn_D   = (/0.0, 0.0/), &
            CoordMaxIn_D   = (/180.0,360.0/), &
            VarIn_IIV = Potential_II)
    end if
    
    ! We only need potential in northern hemisphere, from 10 to 80 degrees.
    nIePhi   = nPhi
    nIeTheta = ratioGrid     * (nTheta-1) + 1
    iTheta10 = ratioGrid/7.0 * (nTheta-1) + 1
    iTheta80 = iTheta10 + nIeTheta - 1
    if (.not. allocated(SwmfIonoPot_II)) &
         allocate(SwmfIonoPot_II(nIeTheta, nIePhi))
    SwmfIonoPot_II(:,:) = Potential_II(iTheta10:iTheta80,:,1)
    
    ! Ram needs SwmfIonoPot in Volts, not kV, so no conversion necessary.
  end subroutine IM_put_from_ie_mpi
  
  !============================================================================
  subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)
    
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    
    implicit none
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd
    character(len=*), parameter   :: NameSub='IM_put_from_ie'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' Not implemented for RAM_SCB.')
  end subroutine IM_put_from_ie
  
  !============================================================================
  
  subroutine IM_put_from_ie_complete
    
    implicit none
    
    character(len=*), parameter :: NameSub = 'IM_put_from_ie_complete'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' Not implemented for RAM_SCB.')
    
  end subroutine IM_put_from_ie_complete
  
  !============================================================================
  
  subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
    real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
    character (len=*),intent(in)       :: NameVar

    character (len=*),parameter :: NameSub = 'IM_put_from_gm'
    
    call CON_stop(NameSub // ' should not be called in IM/RAM_SCB')
    
  end subroutine IM_put_from_gm
  
  !==========================================================================
  
  subroutine IM_put_from_gm_line(nRadiusIn, nLonIn, Map_DSII, &
       nVarLineIn, nPointLineIn, BufferLine_VI, NameVar)
    
    use ModRamMain,  ONLY: nR, nT, TimeRamElapsed, RadiusMin, RadiusMax, &
         Real8_, PathRamOut, nRextend
    use ModNumConst, ONLY: cRadToDeg
    use ModIoUnit,   ONLY: UnitTmp_
    use ModPlotFile, ONLY: save_plot_file
    use ModRamCouple
    
    implicit none
    
    integer, intent(in) :: nRadiusIn, nLonIn
    real,    intent(in) :: Map_DSII(3,2,nRadiusIn,nLonIn)
    integer, intent(in) :: nVarLineIn, nPointLineIn
    real,    intent(in) :: BufferLine_VI(nVarLineIn,nPointLineIn)
    character(len=*), intent(in) :: NameVar
    
    integer           :: iR, iT, iDir, n
    real(kind=Real8_) :: rotate_VII(3,nT,4) = 0.0
    logical, save     :: IsFirstCall = .true.
    
    ! These variables should either be in a module, OR
    ! there is no need for them, and BufferLine_VI should be put 
    ! into RAM_SCB variables right here. 
    ! Note that this routine is only called on the root processor !!!
    integer :: nVarLine   = 0          ! number of vars per line point
    integer :: nPointLine = 0          ! number of points in all lines
    integer, save :: iLine_III(2,nRextend,nT)           ! line index 
    
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='IM_put_from_gm_line'
    
    ! Variables for testing
    integer :: iPoint
    character(len=100):: NameFile
    !---------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    ! Ensure number of points (including radial ghost cell) match.
    if ((nRextend .ne. nRadiusIn) .or. (nT .ne. nLonIn)) then
       write(*,*)"SWMF nR, nT = ", nRadiusIn, nLonIn
       write(*,*)"RAM  nR, nT = ", nRextend,  nT
       call CON_stop(NameSub//': IM gridsize does not match SWMF description.')
    endif
    
    
    ! Save total number of points along all field lines
    nPointLine = nPointLineIn
    nVarLine   = nVarLineIn
    
    ! On first call, set run-specific information.
    if (IsFirstCall) then
       ! Create index array that converts radial 
       ! and local time index to line index.
       n = 0
       do iT = 1, nT; do iR = 1, nRextend; do iDir = 1, 2
          n = n+1
          iLine_III(iDir,iR,iT) = n
       end do; end do; end do
       
       ! These are used for building Euler potential surfaces:
       nRadSWMF = nRextend
       nLonSWMF = nT
       nLinesSWMF = nRextend*nT*2
       
       ! These bind the max latitude of the SCB domain to the 
       ! footpoints of field lines that pass near L=6.75
       nRadSwmfVar      = nR + 2 ! 2nd ghost cell for >6.75.
       nRadSwmfVarInner = nR + 1 ! 1st ghost cell for =6.75.
       
       ! Allocate array to hold sorted line data:
       if(.not. allocated(MhdLines_IIV)) &
            allocate(MhdLines_IIV(nLinesSWMF, nPointsMax, nVarLine))

       ! Use NameVar to determine what type of MHD 
       ! simulation is being used and set indices.
       call set_type_mhd(NameVar, nVarLine)
       
       IsFirstCall = .false.
       
       if(DoTestMe)then
          write(*,*)'IM: ', NameSub, ' nR, nRextend, nLonSWMF, nLinesSWMF ='
          write(*,*)'IM: ', nR, nRextend, nLonSWMF, nLinesSWMF
          write(*,*)'IM: ', NameSub, ' nRadSwmf, nRadSwmfVar, nRadSwmfVarInner='
          write(*,*)'IM: ', nRadSwmf, nRadSwmfVar, nRadSwmfVarInner
       endif
       
    endif
    
    if(DoTest)then
       write(*,*)NameSub,' nVarLine,nPointLine=',nVarLine,nPointLine
       
       ! Set the file name for the ray trace output.
       write(NameFile,'(a,i5.5,a)') &
            PathRamOut//"ray_data_t",nint(TimeRamElapsed),".out"
       open(UnitTmp_, FILE=NameFile, STATUS="replace")
       ! Same format as in GM/BATSRUS/src/ray_trace_new.f90
       write(UnitTmp_, *) 'nRadius, nLon, nPoint=',nRextend, nT, nPointLine
       write(UnitTmp_, *) 'iLine l x y z rho ux uy uz bx by bz p'
       do iPoint = 1, nPointLine
          write(UnitTmp_, *) BufferLine_VI(:, iPoint)
       end do
       close(UnitTmp_)

       ! Now save the mapping files (+0.0 for real precision)
       write(NameFile,'(a,i5.5,a)') &
            PathRamOut//"map_north_t",nint(TimeRamElapsed),".out"
       
       call save_plot_file( &
            NameFile, &
            StringHeaderIn = 'Mapping to northern ionosphere', &
            TimeIn       = TimeRamElapsed+0.0, &
            NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
            CoordMinIn_D = (/RadiusMin+0.0,   0.0/), &
            CoordMaxIn_D = (/RadiusMax+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,1,:,:))
       
       write(NameFile,'(a,i5.5,a)') &
            PathRamOut//"map_south_t",nint(TimeRamElapsed),".out"
       
       call save_plot_file( &
            NameFile, &
            StringHeaderIn = 'Mapping to southern ionosphere', &
            TimeIn       = TimeRamElapsed+0.0, &
            NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
            CoordMinIn_D = (/RadiusMin+0.0,   0.0/), &
            CoordMaxIn_D = (/RadiusMax+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,2,:,:))
    end if
    
    call sort_mhd_lines(nVarLine, nPointLine, BufferLine_VI)
    
    ! Extract BATS plasma density and pressure at the outer boundary.
    ! Need values at first outer ghost cell.
    iR = nR+2; iT=1  !!! CHECK HERE PLZ.
    do iPoint = 1, nPointLine
       ! Search for first point in line that matches
       ! the line index for our current rad/lon pair:
       if (BufferLine_VI(1,iPoint) .lt. iLine_III(1,iR+1,iT)) cycle
       
       if(DoTest) write(*,'(a,1x,i2.2,1x,i4.4,1x,f5.0,1x, f5.3)') &
            "iT, lineExpected, LineFound, R = ", &
            iT, iLine_III(1,iR+1,iT), BufferLine_VI(1,iPoint), &
            sqrt(BufferLine_VI(5,iPoint)**2 + BufferLine_VI(4,iPoint)**2 &
            + BufferLine_VI(3,iPoint)**2)/ 6378100.0
       
       if(iLine_III(1,iR+1,iT) .ne. BufferLine_VI(1,iPoint)) then
          write(*,*) "ERROR: Expected and actual line numbers don't match!"
          write(*,*) "Expected, Actual = ", &
               iLine_III(1,iR,iT), BufferLine_VI(1,iPoint)
          write(*,*) "Likely cause is open field lines in IM domain."
          call CON_stop('IM_put_from_gm_line')
       endif
       
       ! Grab total dens and pressure
       MhdDensPres_VII(1,iT,1) = BufferLine_VI(TotalRho_ , iPoint)
       MhdDensPres_VII(2,iT,1) = BufferLine_VI(TotalPres_, iPoint)
       PMhdGhost = BufferLine_VI(TotalPres_, iPoint)
       
       if (TypeMhd .eq. 'anisoP')then
          ! anisotropic pressure from GM                           
          MhdDensPres_VII(3,iT,1) = BufferLIne_VI(Ppar_,iPoint)
       end if
       
       if (TypeMhd .ne. 'single') then
          MhdDensPres_VII(1,iT,2) = BufferLine_VI(RhoH_ ,  iPoint)
          MhdDensPres_VII(1,iT,4) = BufferLine_VI(RhoO_ ,  iPoint)
          ! Only MultiS has He; only Mult3F has Sw H+.
          ! Store it in Helium's slot.
          if (TypeMhd .eq. 'mult3F') then
             MhdDensPres_VII(1,iT,3) = BufferLine_VI(RhoSw_ , iPoint)
          else
             MhdDensPres_VII(1,iT,3) = BufferLine_VI(RhoHe_ , iPoint)
          end if
       end if
       
       if (TypeMhd .eq. 'multiF') then
          MhdDensPres_VII(2,iT,2) = BufferLine_VI(PresH_,  iPoint)
          MhdDensPres_VII(2,iT,3) = BufferLine_VI(PresHe_, iPoint)
          MhdDensPres_VII(2,iT,4) = BufferLine_VI(PresO_,  iPoint)
       end if
       if (TypeMhd .eq. 'mult3F') then
          MhdDensPres_VII(2,iT,2) = BufferLine_VI(PresH_ , iPoint)
          MhdDensPres_VII(2,iT,3) = BufferLine_VI(PresSw_, iPoint)
          MhdDensPres_VII(2,iT,4) = BufferLine_VI(PresO_ , iPoint)
       end if
       
       iT = iT + 1 ! Cycle to next logitude.
       if (iT .gt. nT) exit
    end do
    
    ! Rotate the local time coordinates by 12 Hours.
    do iT=1, nT
       n = mod(iT+12-1,24)+1
       rotate_VII(:,iT,:) = MhdDensPres_VII(:,n,:)
    end do
    MhdDensPres_VII = rotate_VII
    
    if(DoTest) then
       write(*,*) 'IM_Wrapper: pressure, in Pa:'
       write(*,*) MhdDensPres_VII(2,:,1)
    end if
    
    ! Call routine that converts density/pressure into flux.
    call generate_flux
    
    ! Map_DSII will be used for mapping the pressure to the ionosphere
    ! so that BATS can use it as if it were RCM pressure.
    ! Start here by converting from colat/radiats to lat/degrees.
    ! Rotate the local time coordinates by 12 Hours.
    !IonoMap_DSII =  Map_DSII
    do iT=1, nT
       n = mod(iT+12-1,24)+1
       IonoMap_DSII(:,:,:,iT) = Map_DSII(:,:,:,n)
    end do
    IonoMap_DSII(2,:,:,:) = 90.0 - cRadtoDeg * IonoMap_DSII(2,:,:,:)
    IonoMaP_DSII(3,:,:,:) =        cRadtoDeg * IonoMap_DSII(3,:,:,:)
    
    
  end subroutine IM_put_from_gm_line

  !=============================================================================

  subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
    ! Puts satellite locations and names from GM into IM variables.

    use ModNumConst,   ONLY: cDegToRad
    
    implicit none
    character (len=*),parameter :: NameSub='IM_put_sat_from_gm'
    
    ! Arguments
    integer, intent(in)            :: nSats
    real, intent(in)               :: Buffer_III(3,2,nSats)
    character(len=100), intent(in) :: Buffer_I(nSats)
    !---------------------------------------------------------------------------
    
  end subroutine IM_put_sat_from_gm
  
  !=============================================================================

  subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)
    
    use CON_time,     ONLY: get_time
    use ModRamMain,   ONLY: Real8_, nR, nT, PAllSum, HPAllSum, &
         OPAllSum, PparSum, BNES, NAllSum, HNAllSum, ONAllSum, &
         TimeRamNow, PathRamOut,nRextend, DoAnisoPressureGMCoupling
    use ModRamCouple, ONLY: IonoMap_DSII, MhdDensPres_VII, PMhdGhost, &
         DoMultiFluidGMCoupling
    use ModIoUnit,    ONLY: UNITTMP_
    use ModRamFunctions, ONLY: ram_sum_pressure
    implicit none
    
    character (len=*),parameter :: NameSub='IM_get_for_gm'
    
    integer, intent(in)                                :: iSizeIn,jSizeIn,nVar
    real, dimension(iSizeIn,jSizeIn,nVar), intent(out) :: Buffer_IIV
    character (len=*),intent(in)                       :: NameVar
    
    logical, save     :: IsFirstCall = .true.
    logical :: DoTest, DoTestMe
    real(kind=Real8_) :: cEnerToPa = 1.6E-10 ! KeV/cm3 to Pascals
    integer :: i, iT, iRot
    integer, parameter :: pres_=1, dens_=2, parpres_=3, bmin_=4,&
         Hpres_=3, Opres_=4, Hdens_=5, Odens_=6
    
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    ! Check to ensure that what we get is what GM wants.
    ! Include RAM's ghost cells for continuity with incoming coupling.
    if ( (iSizeIn .ne. nRextend) .or. (jSizeIn .ne. nT) ) then
       write(*,*) 'nR, nT = ', nRextend, nT
       write(*,*) 'iSizeIn, jSizeIn = ', iSizeIn, jSizeIn
       call CON_stop(NameSub//' ERROR: mismatched array sizes!')
    endif
    
    DoMultiFluidGMCoupling = .false.
    DoAnisoPressureGMCoupling = .false.
    
    if(NameVar == 'p:rho:Hpp:Opp:Hprho:Oprho')then
       DoMultiFluidGMCoupling = .true.
    elseif(NameVar == 'p:rho:ppar:bmin')then                             
       DoAnisoPressureGMCoupling = .true.
    elseif (NameVar /= 'p:rho') then
       call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if
    
    ! Sum pressure:
    call ram_sum_pressure()
    
    ! Fill pressure from computational cells; rotate 180.
    do iT=1, nT
       iRot = mod(iT+11,nT-1) + 1
       Buffer_IIV(1:nR,iRot,pres_) = PAllSum(1:nR,iT)
       Buffer_IIV(1:nR,iRot,dens_) = NAllSum(1:nR,iT)
       if (DoAnisoPressureGMCoupling)then
          Buffer_IIV(1:nR,iRot,parpres_) = PparSum(1:nR,iT)
          Buffer_IIV(1:nR,iRot,bmin_)    = BNES(1:nR, iT)
       end if
       if (DoMultiFluidGMCoupling)then
          Buffer_IIV(1:nR,iRot,Hpres_)   = HPAllSum(1:nR,iT)
          Buffer_IIV(1:nR,iRot,Opres_)   = OPAllSum(1:nR,iT)
          Buffer_IIV(1:nR,iRot,Hdens_)   = HNAllSum(1:nR,iT)
          Buffer_IIV(1:nR,iRot,Odens_)   = ONAllSum(1:nR,iT)
       end if
    end do
    
    ! Azimuthal symmetry.
    Buffer_IIV(1:nR,nT,pres_) = Buffer_IIV(1:nR,1,pres_)
    if (DoAnisoPressureGMCoupling)then
       Buffer_IIV(1:nR,nT,parpres_) = Buffer_IIV(1:nR,1,parpres_)
       Buffer_IIV(1:nR,nT,bmin_) = Buffer_IIV(1:nR,1,bmin_)
    end if
    if (DoMultiFluidGMCoupling)then
       Buffer_IIV(1:nR,nT,Hpres_)   = Buffer_IIV(1:nR,1,Hpres_)
       Buffer_IIV(1:nR,nT,Opres_)   = Buffer_IIV(1:nR,1,Opres_)
       Buffer_IIV(1:nR,nT,Hdens_)   = Buffer_IIV(1:nR,1,Hdens_)
       Buffer_IIV(1:nR,nT,Odens_)   = Buffer_IIV(1:nR,1,Odens_)
    end if
    
    ! Convert Units.
    Buffer_IIV(:,:,pres_) = Buffer_IIV(:,:,pres_)*cEnerToPa
    if (DoAnisoPressureGMCoupling) &
         Buffer_IIV(:,:,parpres_) = Buffer_IIV(:,:,parpres_)*cEnerToPa
    
    ! Fill pressure from ghost cells as -1 (i.e., no coupling).
    Buffer_IIV(nR+1:nRextend,:,pres_) = -1
    if (DoAnisoPressureGMCoupling)then
       Buffer_IIV(nR+1:nRextend,:,parpres_) = -1
       Buffer_IIV(nR+1:nRextend,:,bmin_) = -1
    end if
    if (DoMultiFluidGMCoupling)then
       Buffer_IIV(nR+1:nRextend,:,Hpres_) = -1
       Buffer_IIV(nR+1:nRextend,:,Opres_) = -1
       Buffer_IIV(nR+1:nRextend,:,Hdens_) = -1
       Buffer_IIV(nR+1:nRextend,:,Odens_) = -1
    end if
    
    if(DoTestMe)then
       write(*,*)'Maxvals of Buffer_IIV, PAllSum = '
       do i=1, nR
          write(*,*) maxval(Buffer_IIV(i,:,1)), maxval(PAllSum(i, :)*cEnerToPa)
       end do
    end if
    
  end subroutine IM_get_for_gm

  !============================================================================
  
  subroutine IM_init_session(iSession, TimeSimulation)
    
    use ModRamMain, ONLY: iCal, IsComponent, &
         TimeRamStart, StrRamDescription, TimeMax
    use CON_time,   ONLY: get_time, tSimulationMax
    !use CON_variables, ONLY: StringDescription
    implicit none
    
    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    
    character(len=*), parameter :: NameSub='IM_init_session'
    logical :: DoTest, DoTestMe
    !---------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' called for iSession, TimeSimulation =', &
         iSession, TimeSimulation
    
    ! Set iCal to iteration = 1.
    iCal = 1
    
    ! Set RAM-SCB to component mode
    IsComponent = .true.
    
    ! Get run description from framework.
    StrRamDescription = "RAM_SCB in coupled mode" !StringDescription
    
    ! Set starting time, using start time of simulation.
    call get_time(TimeStartOut=TimeRamStart)
    
    ! Set TimeMax from SWMF
    TimeMax = tSimulationMax
    
    call init_ramscb
    
  end subroutine IM_init_session
  
  !=============================================================================
  
  subroutine IM_finalize(TimeSimulation)
    
    implicit none
    
    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    !---------------------------------------------------------------------------
    call finalize_ramscb
  end subroutine IM_finalize
  
  !=============================================================================
  
  subroutine IM_save_restart(TimeSimulation)
    
    use ModRamIO, ONLY: write_restart
    
    implicit none
    
    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    character(len=*), parameter :: NameSub='IM_save_restart'
    !---------------------------------------------------------------------------
    
    call write_restart
    
  end subroutine IM_save_restart
  
  !BOP =========================================================================
  !ROUTINE: IM_run - run IM
  !INTERFACE:
  subroutine IM_run(TimeSimulation,TimeSimulationLimit)
    
    use ModRamMain, ONLY: DtsMax, DtsFramework, TimeRamElapsed
    
    implicit none
    
    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout) :: TimeSimulation   ! current time of component
    !INPUT ARGUMENTS:
    real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
    
    !LOCAL VARIABLES:
    real :: Dt, DtTotal
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='IM_run'
    !---------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub, &
         ' called for TimeSimulation, TimeSimulationLimit=', &
         TimeSimulation, TimeSimulationLimit
    
    ! Ram does 2 timesteps per call.
    ! Time step should 1/2 of total limit, up to DtsMax.
    DtTotal = TimeSimulationLimit - TimeSimulation
    
    ! Multiply with 1.0 to avoid mixed real precision 
    DtsFramework = min(1.0*DtsMax, 0.5 * DtTotal)
    TimeRamElapsed = TimeSimulation
    
    call run_ramscb
    
    TimeSimulation = TimeRamElapsed
    
  end subroutine IM_run
  
  !=============================================================================
  subroutine IM_put_from_gm_crcm(Integral_IIV, iSize, jSize, nIntegral, &
       Bufferline_VI, nVarLine, nPointLine, NameVar, tSimulation)
    ! Dummy routine to satisfy dependencies in SWMF.
    use ModRamMain, ONLY: Real8_
    
    implicit none
    
    integer,           intent(in) :: iSize, jSize, nIntegral
    integer,           intent(in) :: nVarLine, nPointLine
    real(kind=Real8_), intent(in) :: Integral_IIV(iSize,jSize,nIntegral)
    real(kind=Real8_), intent(in) :: BufferLine_VI(nVarLine, nPointLine)
    character(len=100),intent(in) :: NameVar
    real(kind=Real8_), intent(in) :: tSimulation
    
    character(len=*), parameter :: NameSub = 'IM_putfrom_gm_crcm'
    !---------------------------------------------------------------------------
    call CON_stop(NameSub//': Do not call this routine from RAM-SCB!')
    
  end subroutine IM_put_from_gm_crcm

end module IM_wrapper

!=============================================================================
subroutine IM_write_prefix
  
  implicit none
  
  !---------------------------------------------------------------------------
  
end subroutine IM_write_prefix
  
