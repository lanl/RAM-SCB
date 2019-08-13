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
            NameVersion='RAM-SCB (Jordanova, Zaharia, Engel)', &
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
    use ModNumConst,     ONLY: cTwoPi
    use ModRamGrids,     ONLY: RadiusMin, RadiusMax, NR, NT, dPhi, dR, nRextend
    use ModRamVariables, ONLY: LZ, PHI, DL1, GridExtend
    use CON_coupler,     ONLY: set_grid_descriptor, is_proc, IM_
    implicit none
    
    character (len=*), parameter :: NameSub='IM_set_grid'
    logical :: IsInitialized=.false.
    logical :: DoTest, DoTestMe
    integer :: i, j, iS
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)'IM_set_grid called, IsInitialized=', &
         IsInitialized
    if(IsInitialized) return

    IsInitialized = .true.
    
    ! IM grid: the equatorial grid is described by Coord1_I and Coord2_I
    ! Occasional +0.0 is used to convert from single to double precision
    ! Grid information set here is shared with GM for coupling with BATS-R-US.

    ! If we are on an IM block, we can use the RAM-SCB run-time generated grid.
    if(is_proc(IM_))then
       call set_grid_descriptor( IM_,           & ! component index
            nDim     = 2,                       & ! dimensionality
            nRootBlock_D = (/1,1/),             & ! number of blocks
            nCell_D  =(/nRextend, nT/),         & ! size of equatorial grid
            XyzMin_D =(/RadiusMin+0.0,0.0/),    & ! min coordinates
            XyzMax_D =(/GridExtend(nRextend),cTwoPi/), & ! max coordinates+radial ghost
            Coord1_I = GridExtend+0.0,          & ! radial coordinates
            Coord2_I = Phi(1:nT)+0.0,           & ! longitudinal coordinates
            TypeCoord= 'SMG',                   & ! solar magnetic coord
            nVar=2,                             & ! Number of "fluid" vars.
            NameVar = 'p rho')
            !NameVar = 'rho p ppar Hpp Opp Hprho Oprho')
       ! For nVar/NameVar, RAM-SCB+BATS-R-US coupling does not need this
       ! convention.  Placeholders are used ONLY.
    else
       ! Not on IM? Use dummy vars as we have not allocted grid.
       call set_grid_descriptor( IM_,           & ! component index
            nDim     = 2,                       & ! dimensionality
            nRootBlock_D = (/1,1/),             & ! number of blocks
            nCell_D  =(/nRextend, nT/),         & ! size of equatorial grid
            XyzMin_D =(/RadiusMin+0.0,0.0/),    & ! min coordinates
            XyzMax_D =(/6.5,cTwoPi/),           & ! max coordinates+radial ghost
            Coord1_I =(/RadiusMin, RadiusMax/), & ! radial coordinates
            Coord2_I =(/0.0, cTwoPi/),          & ! longitudinal coordinates
            TypeCoord='SMG',                    & ! solar magnetic coord
            nVar=2,                             & ! Number of "fluid" vars.
            NameVar = 'p rho')
            !NameVar = 'rho p ppar Hpp Opp Hprho Oprho')
    endif

       
    if(DoTestMe)then
       write(*,*)NameSub,' Setting RAM-SCB grid characteristics.'
       write(*,*)NameSub,' Are we on an IM proc?  ', is_proc(IM_)
    end if
    if(DoTestMe .and. is_proc(IM_)) then
       write(*,*)NameSub,' NR =              ', NR
       write(*,*)NameSub,' NR + Ghostcells = ', nRextend
       write(*,*)NameSub,' NT =              ', NT
       write(*,*)NameSub,' Phi grid = ', 360.*Phi/cTwoPi
       write(*,*)NameSub,' Rad grid = ', GridExtend
    end if
    
  end subroutine IM_set_grid

  !============================================================================
  subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)
    
    ! Provide current for IE
    ! The value should be interpolated from nPoints with
    ! indexes stored in Index and weights stored in Weight
    ! The variables should be put into Buff_V
    ! This is not implemented for RAM-SCB.  IM to IE coupling occurs
    ! implicitly through IM-to-GM pressure coupling.
    
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    implicit none
    character(len=*), parameter :: NameSub='IM_get_for_ie'
    
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//' Not implemented for RAM_SCB.')
    
  end subroutine IM_get_for_ie

  !============================================================================
  subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)
    
    use ModRamMain,   ONLY: Real8_,  PathRamOut
    use ModRamTiming, ONLY: TimeRamElapsed
    use ModRamCouple, ONLY: SwmfIonoPot_II, nIePhi, nIeTheta
    use ModPlotFile,  ONLY: save_plot_file
    
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
  
  subroutine IM_put_from_gm(Buffer_IIV,Kp,iSizeIn,jSizeIn,nVarIn,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
    real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
    real, intent(in) :: Kp
    character (len=*),intent(in)       :: NameVar

    character (len=*),parameter :: NameSub = 'IM_put_from_gm'
    
    call CON_stop(NameSub // ' should not be called in IM/RAM_SCB')
    
  end subroutine IM_put_from_gm
  
  !==========================================================================
  
  subroutine IM_put_from_gm_line(nRadiusIn, nLonIn, Map_DSII, &
       nVarLineIn, nPointLineIn, BufferLine_VI, NameVar)
    
    use ModRamMain,   ONLY: Real8_, PathRamOut
    use ModRamGrids,  ONLY: nR, nT, RadiusMin, RadiusMax, RadMaxScb, nRextend
    use ModRamVariables,ONLY: GridExtend, Phi
    use ModRamTiming, ONLY: TimeRamElapsed
    use ModNumConst,  ONLY: cRadToDeg
    use ModIoUnit,    ONLY: UnitTmp_
    use ModPlotFile,  ONLY: save_plot_file
    use ModRamCouple
    use ModNumConst,  ONLY: cTwoPi
    
    implicit none
    
    integer, intent(in) :: nRadiusIn, nLonIn
    real,    intent(in) :: Map_DSII(3,2,nRadiusIn,nLonIn)
    integer, intent(in) :: nVarLineIn, nPointLineIn
    real,    intent(in) :: BufferLine_VI(nVarLineIn,nPointLineIn)
    character(len=*), intent(in) :: NameVar
    
    integer           :: iT, iR, iDir, iRot, iLine
    real(kind=Real8_) :: rotate_VII(3,nT,4), rad(10000)
    logical, save     :: IsFirstCall = .true.

    real, parameter :: sens=1E-5  ! Sensitivity for real differences.
    
    ! These variables should either be in a module, OR
    ! there is no need for them, and BufferLine_VI should be put 
    ! into RAM_SCB variables right here. 
    ! Note that this routine is only called on the root processor !!!
    integer :: nVarLine   = 0          ! number of vars per line point
    integer :: nPointLine = 0          ! number of points in all lines
    
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

    ! Initialize rotate array:
    rotate_VII(3,nT,4) = 0.0
    
    ! Save total number of points along all field lines
    nPointLine = nPointLineIn
    nVarLine   = nVarLineIn
    
    ! On first call, set run-specific information.
    if (IsFirstCall) then
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
            CoordMaxIn_D = (/RadMaxScb+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,1,:,:))
       
       write(NameFile,'(a,i5.5,a)') &
            PathRamOut//"map_south_t",nint(TimeRamElapsed),".out"
       
       call save_plot_file( &
            NameFile, &
            StringHeaderIn = 'Mapping to southern ionosphere', &
            TimeIn       = TimeRamElapsed+0.0, &
            NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
            CoordMinIn_D = (/RadiusMin+0.0,   0.0/), &
            CoordMaxIn_D = (/RadMaxScb+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,2,:,:))
    end if
    
    call sort_mhd_lines(nVarLine, nPointLine, BufferLine_VI)

    ! Extract MHD density and pressure at the outer ghost cell of RAM-SCB
    ! This is at index nR+1 (one outer ghost cell).
    do iT=1, nT
       ! Get the number of the field line trace passing through the
       ! RAM-SCB ghost cell center:
       iLine = 2*(nRextend*(iT-1)+nR+1) - 1 ! Northern hemi trace.

       ! Write some debug info to screen:
       if(DoTest) write(*,'(a,1x,i2.2,1x,i4.4,4(2x,f6.3))') &
            "iT, iLine, X, Y, Z, R = ", iT, iLine, &
            MhdLines_IIV(iLine, 1, 3), &
            MhdLines_IIV(iLine, 1, 4), &
            MhdLines_IIV(iLine, 1, 5), &
            sqrt( MhdLines_IIV(iLine, 1, 3)**2   &
            +     MhdLines_IIV(iLine, 1, 4)**2   &
            +     MhdLines_IIV(iLine, 1, 5)**2)
      
       select case (TypeMhd)
       case('single') 
          ! Grab total dens and pressure, use first point in trace.
          if (iEnd(iLine) == 0) then
             MhdDensPres_VII(1,iT,0) = -1.0
             MhdDensPres_VII(2,iT,0) = -1.0
          else
             MhdDensPres_VII(1,iT,0) = MhdLines_IIV(iLine, 1, TotalRho_)
             MhdDensPres_VII(2,iT,0) = MhdLines_IIV(iLine, 1, TotalPres_)
          endif
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
    end do

    ! Rotate the local time coordinates by 12 Hours.
    do iT=1, nT
       iRot = mod(iT+((nT-1)/2-1),nT-1) + 1 ! Even # of cells, 1 ghost cell.
       rotate_VII(:,iT,:) = MhdDensPres_VII(:,iRot,:)
    end do
    MhdDensPres_VII(:,1:nT,:) = rotate_VII(:,nT:1:-1,:)
    
    if(DoTest) then
       write(*,*) 'IM_Wrapper: pressure, in Pa:'
       write(*,*) MhdDensPres_VII(2,:,0)
    end if
    
    ! Call routine that converts density/pressure into flux.
    call generate_flux
    
    ! DEPRECIATED: This info no longer saved.
    ! Map_DSII will be used for mapping the pressure to the ionosphere
    ! so that BATS can use it as if it were RCM pressure.
    ! Start here by converting from colat/radiats to lat/degrees.
    ! Rotate the local time coordinates by 12 Hours.
    !do iT=1, nT
    !   n = mod(iT+((nT-1)/2-1),nT-1) + 1 ! Even # of cells, 1 ghost cell.
    !   IonoMap_DSII(:,:,:,iT) = Map_DSII(:,:,:,n)
    !end do
    !IonoMap_DSII(2,:,:,:) = 90.0 - cRadtoDeg * IonoMap_DSII(2,:,:,:)
    !IonoMap_DSII(3,:,:,:) =        cRadtoDeg * IonoMap_DSII(3,:,:,:)


    ! Clear Blines_DIII and reallocate to match what is necessary:
    if(allocated(Blines_DIII)) deallocate(Blines_DIII)
    allocate(Blines_DIII(6, nRextend, nT-1, 2*nPoints-1))
    Blines_DIII = 0
    
    ! Sort lines one more time- this time for use by SCB.
    do iT=1, nT-1
       do iR=1, nRextend
          ! Get the index that rotates in longitude 180 degrees.
          iRot = iT !mod(iT+((nT-1)/2-1),nT-1) + 1 ! Even # of cells, 1 ghost cell.
          
          ! Check if the line is open:
          iLine = 2*(nRextend*(iT-1)+iR)
          if ((iEnd(iLine) == 0).or.(iEnd(iLine-1) == 0)) then
             IsClosed_II(iR,iT) = .false.
             cycle
          endif
          if((Map_DSII(1,1,iR,iT)-1>sens).or.(Map_DSII(1,2,iR,iT)-1>sens)) then
             IsClosed_II(iR,iT) = .false.
             if(DoTest) then
                write(*,*)'OPEN FIELD LINE AT R, MLT = ', &
                     GridExtend(iR), Phi(iRot)*24./cTwoPi
                write(*,*)'Line ends at R (north, south) = ', Map_DSII(1,:,iR,iT)
             end if
             cycle
          end if
          isClosed_II(iR,iT) = .true.

          ! Fill points from northern hemisphere:
          iLine = iLine - 1
          
          Blines_DIII(1,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,3) !X
          Blines_DIII(2,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,4) !Y
          Blines_DIII(3,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,5) !Z
          Blines_DIII(4,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,Bx_) !Bx
          Blines_DIII(5,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,By_) !By
          Blines_DIII(6,iR,iRot,nPoints:2*nPoints-1) = MhdLines_IIV(iLine,1:nPoints,Bz_) !Bz


          if(DoTest) write(*,*) 'Northern Hemi trace = ', &
               sqrt(MhdLines_IIV(iLine,1:nPoints,3)**2 + &
                    MhdLines_IIV(iLine,1:nPoints,4)**2 + &
                    MhdLines_IIV(iLine,1:nPoints,5)**2)
          
          ! Fill points from southern hemisphere:
          iLine = iLine + 1
          Blines_DIII(1,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,3) !X
          Blines_DIII(2,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,4) !Y
          Blines_DIII(3,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,5) !Z
          Blines_DIII(4,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,Bx_) !Bx
          Blines_DIII(5,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,By_) !By
          Blines_DIII(6,iR,iRot,1:nPoints) = MhdLines_IIV(iLine,nPoints:1:-1,Bz_) !Bz

          rad(1:2*nPoints-1) = sqrt(Blines_DIII(1,iR,iRot,:)**2 &
                                   +Blines_DIII(2,iR,iRot,:)**2 &
                                   +Blines_DIII(3,iR,iRot,:)**2)
          
          if(DoTest) write(*,*) 'Line runs from R=', rad(1), rad(2*nPoints-1)

       end do
    end do
    
    
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
    
    use CON_time,        ONLY: get_time
    use ModRamMain,      ONLY: Real8_, PathRamOut
    use ModRamParams,    ONLY: DoAnisoPressureGMCoupling
    use ModRamGrids,     ONLY: nR, nT, nRextend
    use ModRamVariables, ONLY: PAllSum, HPAllSum, OPAllSum, PparSum, &
         NAllSum, HNAllSum, ONAllSum, BNES
    use ModRamTiming,    ONLY: TimeRamNow
    use ModRamCouple,    ONLY: MhdDensPres_VII, DoMultiFluidGMCoupling
    use ModIoUnit,       ONLY: UNITTMP_
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

    ! Note that the coupler asks for values out to nRextend, or the
    ! radial extent of SCB.  The radial extent of RAM is nR which is <nRextend.
    ! We fill pressure values up to nR and leave the outer cells as negative
    ! values.  BATS will not couple negative values.
    
    ! Check to ensure that what we get is what GM wants.
    ! Include RAM's ghost cells for continuity with incoming coupling.
    if ( (iSizeIn .ne. nRextend) .or. (jSizeIn .ne. nT) ) then
       write(*,*) 'nR, nT = ', nRextend, nT
       write(*,*) 'iSizeIn, jSizeIn = ', iSizeIn, jSizeIn
       call CON_stop(NameSub//' ERROR: mismatched array sizes!')
    endif

    ! Initialize all values to -1; MHD will not couple negative values.
    Buffer_IIV=-1
    
    DoMultiFluidGMCoupling = .false.
    DoAnisoPressureGMCoupling = .false.

    select case(NameVar)
    case('p:rho')
       ! Sum pressure:
       call ram_sum_pressure()

       ! Fill pressure from computational cells; rotate 180.
       do iT=1, nT
          iRot = mod(iT+((nT-1)/2-1),nT-1) + 1 ! Even # of cells, 1 ghost cell.
          Buffer_IIV(1:nR,iRot,pres_) = PAllSum(1:nR,iT)
          Buffer_IIV(1:nR,iRot,dens_) = NAllSum(1:nR,iT)
       end do
       Buffer_IIV(1:nR,nT,pres_) = Buffer_IIV(1:nR,1,pres_)
       Buffer_IIV(1:nR,nT,dens_) = Buffer_IIV(1:nR,1,dens_)
       Buffer_IIV(:,:,pres_) = Buffer_IIV(:,:,pres_)*cEnerToPa
       Buffer_IIV(nR+1:nRextend,:,pres_) = -1

    case('p:rho:Hpp:Opp:Hprho:Oprho')
       call CON_stop(NameSub//' NameVar not currently supported: '//trim(NameVar))
       DoMultiFluidGMCoupling = .true.
    case('p:rho:ppar:bmin')
       call CON_stop(NameSub//' NameVar not currently supported: '//trim(NameVar))
       DoAnisoPressureGMCoupling = .true.
    case default
       call CON_stop(NameSub//' invalid NameVar='//NameVar)
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
    
    use ModRamMain,   ONLY: iCal
    use ModRamTiming, ONLY: TimeRamStart, TimeMax
    use ModRamParams, ONLY: IsComponent, StrRamDescription, IsComponent
    use ModRamCouple, ONLY: RAMCouple_Allocate
    use ModRamGsl,    ONLY: gsl_initialize
    use ModRamSpecies, ONLY: DefineSpecies
    use ModRamInit,   ONLY: ram_init, ram_allocate, init_input
    use ModScbInit,   ONLY: scb_init, scb_allocate
    use ModSceInit,   ONLY: sce_init, sce_allocate
    use ModRamScb,    ONLY: ramscb_allocate
    use ModRamIO,     ONLY: init_output
    use CON_time,     ONLY: get_time, tSimulationMax
    
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

    if(DoTest) write(*,'(a,i4, 2("-",i2.2),1x, i2.2,2(":",i2.2)))') &
         'IM: TimeRamStart = ', &
         TimeRamStart%iYear, TimeRamStart%iMonth,  TimeRamStart%iDay, &
         TimeRamStart%iHour, TimeRamStart%iMinute, TimeRamStart%iSecond
    
    
    ! Set TimeMax from SWMF
    TimeMax = tSimulationMax

    if(DoTest) write(*,'(a, E12.6, a)')'IM: Simulation duration (s) = ', TimeMax
    
    ! Allocate Arrays
    call DefineSpecies
    call ram_allocate
    call scb_allocate
    call ramscb_allocate
    call sce_allocate

    ! Initialize RAM_SCB
    call ram_init
    call scb_init
    call sce_init
    call GSL_Initialize

    ! Allocate arrays for coupling:
    call RAMCouple_Allocate

    ! Initialize input and output:
    call init_input
    call init_output

  end subroutine IM_init_session
  
  !=============================================================================
  
  subroutine IM_finalize(TimeSimulation)

    use ModRamCouple, ONLY: RAMCouple_Deallocate
    use ModRamInit,   ONLY: ram_deallocate
    use ModScbInit,   ONLY: scb_deallocate
    use ModSceInit,   ONLY: sce_deallocate
    use ModRamScb,    ONLY: ramscb_deallocate

    implicit none
    
    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    !---------------------------------------------------------------------------
    
    ! Deallocate arrays
    call ram_deallocate
    call scb_deallocate
    call ramscb_deallocate
    call sce_deallocate
    call RAMCouple_Deallocate

  end subroutine IM_finalize
  
  !=============================================================================
  
  subroutine IM_save_restart(TimeSimulation)
    
    use ModRamRestart, ONLY: write_restart
    
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
    
    use ModRamTiming, ONLY: DtsMax, DtsFramework, TimeRamElapsed, TimeRestart
    use ModRamScbRun, ONLY: run_ramscb
    
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
    DtsFramework   = min(1.0*DtsMax, 0.5 * DtTotal)
    TimeRamElapsed = TimeSimulation
    
    call run_ramscb
    
    TimeSimulation = TimeRamElapsed
    
  end subroutine IM_run
  
  !=============================================================================
  subroutine IM_put_from_gm_crcm(Integral_IIV, Kp, iSizeIn, jSizeIn, &
       nIntegralIn, BufferLine_VI, nVarLine, nPointLine, NameVar, tSimulation)

    integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
    real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
    real,    intent(in) :: Kp
    integer, intent(in) :: nVarLine, nPointLine
    real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)
    real,    intent(in) :: tSimulation
    character (len=*), intent(in) :: NameVar

    character (len=*), parameter :: NameSub='IM_put_from_gm_crcm'

    !---------------------------------------------------------------------------
    call CON_stop(NameSub//': Do not call this routine from RAM-SCB!')
    
  end subroutine IM_put_from_gm_crcm

end module IM_wrapper

!=============================================================================
subroutine IM_write_prefix
  
  implicit none
  
  !---------------------------------------------------------------------------
  
end subroutine IM_write_prefix
  
