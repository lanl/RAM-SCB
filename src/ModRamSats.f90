!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamSats

  use ModRamMain, ONLY: DP

  implicit none

  private !except...

  public:: read_sat_params
  public:: read_sat_input
  public:: init_sats
  public:: fly_sats
  integer, public :: nRamSat=0               ! Default is no virtual sats.

  ! Info about the sat input files.
  integer, parameter :: MaxRamSat=100, MaxRamSatLines = 100000
  integer            :: nSatPoints_I(MaxRamSat)
  integer            :: iSatTime_I(MaxRamSat) = 1
  character(len=200) :: SatFileName_I(MaxRamSat), SatName_I(MaxRamSat)
  character(len=200) :: SatFileName_O(MaxRamSat)
  character(len=3)   :: TypeSatCoord_I(MaxRamSat)

  ! Info about the sat output netcdf files.
  integer :: iSatRecord(MaxRamSat) = 1  ! Location in netCDF file.
  
  ! Variables that contain the time and location of the sats.
  real(DP), allocatable :: SatTraj_IID(:,:,:), SatTime_II(:,:)

  ! Bad data flag (sat outside domain.)
  real(DP)  :: BadDataFlag=-1.0E10

  character(len=3), parameter :: TypeCoordSystem = 'SMG'

  character(len=*), parameter :: NameMod = 'ModRamSats'

  contains

!============================================================================
  subroutine read_sat_params

    use ModReadParam
    use ModRamParams, ONLY: DoSaveRamSats
    use ModRamTiming, ONLY: DtWriteSat

    implicit none

    logical :: DoTest, DoTestMe
    integer :: iSat, l1, l2
    character(len=*), parameter :: NameSub = NameMod // '::read_sat_params'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call read_var('DtOutput', DtWriteSat)
    if(DtWriteSat < 10.0) DtWriteSat = 10.0

    call read_var('nSatellite', nRamSat)
    if (nRamSat > MaxRamSatLines) &
         call CON_stop(NameSub // ': Number of sat files exceeds maximum.')
    if (nRamSat .eq. 0) then
       return
    else
       DoSaveRamSats = .true.
    end if

    do iSat=1, nRamSat
       call read_var('NameTrajectoryFile', SatFileName_I(iSat))
       
       ! Trim full filename down to just satellite name.
       l1 = index(SatFileName_I(iSat), '/', back=.true.) + 1
       l2 = index(SatFileName_I(iSat), '.') - 1
       if (l1-1 <= 0) l1 = 1
       if (l2+1 <= 0) l2 = len_trim(SatFileName_I(iSat))
       SatName_I(iSat) = SatFileName_I(iSat)(l1:l2)
    end do

    if(DoTest) then
       write(*,'(a,i3.3,a)') NameSub // ': Found ', &
            nRamSat, ' virtual satellites in PARAM.in'
       write(*,*) 'Satellite Name (Trajectory File)'
       write(*,*) '---------------------------------'
       do iSat=1, nRamSat
          write(*,"(a,' (',a,')')") trim(SatName_I(iSat)), &
               trim(SatFileName_I(iSat))
       end do
    end if

  end subroutine read_sat_params

!============================================================================
  subroutine read_sat_input
    ! Based on (stolen from?) the version used in BATS-R-US (Umich, Toth).
    ! Note that, for the time being, parallel execution is disabled.

    use ModRamTiming, ONLY: TimeRamStart

    use ModTimeConvert, ONLY: time_int_to_real
    use ModIoUnit,      ONLY: UnitTmp_
    use CON_axes

    implicit none

    integer :: iError, i, iSat , nPoint

    ! One line of input
    character (len=100) :: line

    integer :: iTime_I(7), MaxPoint
    real(DP) :: Xyz_D(3), DateTime
    real(DP), allocatable :: Time_I(:), Xyz_DI(:,:)
    character(len=200):: NameFile

    character(len=*), parameter :: NameSub = NameMod // '::read_sat_input'

    logical :: DoTest, DoTestMe

    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
!    if(iProc == 0)then ! Run only on head node.
       SATELLITES1: do iSat=1, nRamSat
          NameFile = SatFileName_I(iSat)
          open(UnitTmp_, file=NameFile, status="old", iostat = iError)
          if (iError /= 0) call CON_stop(NameSub // &
               ' ERROR1: unable to open file ' // trim(NameFile))
          nPoint = 0

          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE1: do
             read(UnitTmp_,'(a)', iostat = iError ) line
             if (iError /= 0) EXIT READFILE1
             if(index(line,'#START')>0)then
                READPOINTS1: do
                   read(UnitTmp_,*, iostat=iError) iTime_I, Xyz_D
                   if (iError /= 0) EXIT READFILE1
                   ! Add new point
                   nPoint = nPoint + 1
                end do READPOINTS1
             end if
          end do READFILE1
          close(UnitTmp_)
          MaxPoint = max(MaxPoint, nPoint)
       end do SATELLITES1

!    end if ! End running on head node only.

    ! Tell all processors the maximum number of points
    !call MPI_Bcast(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(3, MaxPoint))
    allocate(SatTraj_IID(nRamSat, MaxPoint, 3))
    allocate(SatTime_II(nRamSat, MaxPoint))
    Time_I = 0.0; XYZ_DI = 0.0; SatTraj_IID = 0.0; SatTime_II = 0.0

    ! Read the trajectories
    SATELLITES: do iSat=1, nRamSat

       ! Read file on the root processor
!       if (iProc == 0) then

          NameFile = SatFileName_I(iSat)

          if(DoTest)then
             write(*,*) 'IM:',NameSub, " reading: ",trim(NameFile)
          end if

          open(UnitTmp_, file=NameFile, status="old", iostat = iError)
          if (iError /= 0) call CON_stop(NameSub // &
               ' ERROR: unable to open file ' // trim(NameFile))
          nPoint = 0

          ! Read the file: read #COOR TypeCoord, #START and points
          ! Default coordinate system for RAM-SCB is SM!!
          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE: do

             read(UnitTmp_,'(a)', iostat = iError ) line

             if (iError /= 0) EXIT READFILE

             if(index(line,'#COOR')>0) &
                  read(UnitTmp_,'(a)') TypeSatCoord_I(iSat)

             if(index(line,'#START')>0)then

                READPOINTS: do

                   read(UnitTmp_,*, iostat=iError) iTime_I, Xyz_D

                   if (iError /= 0) EXIT READFILE

                   ! Add new point
                   nPoint = nPoint + 1

                   ! Store coordinates
                   Xyz_DI(:,nPoint) = Xyz_D

                   ! Convert integer date/time to simulation time
                   call time_int_to_real(iTime_I, DateTime)
                   Time_I(nPoint) = DateTime - TimeRamStart % Time

                enddo READPOINTS

             endif

          enddo READFILE

          close(UnitTmp_)

          if(DoTest)write(*,*) NameSub,': nPoint=',nPoint

          ! Convert the coordinates if necessary
          if(DoTest) write(*,*)'File Coord System is ', TypeSatCoord_I(iSat)
          if(TypeSatCoord_I(iSat) /= TypeCoordSystem)then
             if(DoTest) then 
                write(*,*) 'Switching coords to ', TypeCoordSystem
                write(*,'(a, 3f10.7)') 'Before:', Xyz_DI(:,1)
             end if
             do i = 1, nPoint
                Xyz_DI(:,i) = matmul( &
                     transform_matrix( Time_I(i), &
                     TypeSatCoord_I(iSat), TypeCoordSystem), Xyz_DI(:,i) )
             end do
             if(DoTest) then
                write(*,'(a, 3f10.7)') 'After: ', Xyz_DI(:,1)
                write(*,*) 'Example transformation matrix:'
                write(*,*) transform_matrix( Time_I(i), &
                     TypeSatCoord_I(iSat), TypeCoordSystem)
             end if
          end if

!       end if

       ! Tell the number of points to the other processors
       !call MPI_Bcast(nPoint, 1, MPI_INTEGER, 0, iComm, iError)
       nSatPoints_I(iSat) = nPoint

       ! Tell the other processors the satellite time
       !call MPI_Bcast(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

       ! Tell the other processors the coordinates
       !call MPI_Bcast(Xyz_DI, 3*nPoint, MPI_REAL, 0, iComm, iError)

       ! Store time and positions for satellite iSat on all PE-s
       SatTime_II(iSat, 1:nPoint) = Time_I(1:nPoint)
       do i = 1, nPoint
          SatTraj_IID(iSat, i, :) = Xyz_DI(:, i)
       end do

       if(DoTest)then
          nPoint = min(10,nPoint)
          write(*,*) NameSub,': tSat=', SatTime_II( iSat, 1:nPoint)
          write(*,*) NameSub,': xSat=', SatTraj_IID(iSat, 1:nPoint,1)
          write(*,*) NameSub,': ySat=', SatTraj_IID(iSat, 1:nPoint,2)
          write(*,*) NameSub,': zSat=', SatTraj_IID(iSat, 1:nPoint,3)
       end if

    end do SATELLITES

    deallocate(Time_I, Xyz_DI)

  end subroutine read_sat_input

!============================================================================
  subroutine init_sats    
    ! Either create new virtual satellite output files or resume writing
    ! to existing files upon restart.  Collect file information for future
    ! writing.

    use ModRamMain,      ONLY: PathRamOut
    use ModRamTiming,    ONLY: TimeRamRealStart
    use ModRamFunctions, ONLY: RamFileName

    use ModTimeConvert,  ONLY: time_real_to_int, TimeType

    implicit none

    character(len=200) :: SatFileName
    integer :: i

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = NameMod // '::init_sat_files'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    do i=1, nRamSat
      SatFileName = RamFileName(SatName_I(i),'nc',TimeRamRealStart)
      ! Replace now-useless input file name with output file name.
      SatFileName_O(i) = trim(PathRamOut)//trim(SatFileName)
      call create_sat_file(SatFileName_O(i))
    end do

  end subroutine init_sats

!============================================================================
  subroutine create_sat_file(FileNameIn)
    ! Build a new netCDF satellite output file.  Enter all meta data into
    ! the file.

    use netcdf
    use ModTimeConvert
    use ModRamTiming,    ONLY: DtWriteSat
    use ModRamGrids,     ONLY: nS, NE, NPA
    use ModRamVariables, ONLY: EKEV, WE, WMU, MU, species

    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts

    implicit none

    character(len=200), intent(in) :: FileNameIn

    integer :: iS
    integer :: iFileID, iStatus, iTimeDim, iEnDim, iPaDim, iXyzDim, iBoundDim
    integer :: iTimeVar, iXyzVar, iBVar, iEgridVar, iEwidVar, iPgridVar, iBeVar
    integer :: iFluxVar, iOFluxVar, iDtVar, iEcVar, iPwidVar, iFlagVar

    logical :: DoTest, DoTestMe
    character(len=100) :: NameSub = NameMod // '::create_sat_file'
    !-----------------------------------------------------------------------
    call CON_set_do_test(trim(NameSub), DoTest, DoTestMe)
    
    if(DoTest) write(*,*) trim(NameSub), &
         ' Creating NetCDF File ', trim(FileNameIn)

    ! Open file; check for success.
    iStatus = nf90_create(trim(FileNameIn), nf90_clobber, iFileID)

    ! Create dimensions for all variables.
    iStatus = nf90_def_dim(iFileID, 'xyz',         3,   iXyzDim)
    iStatus = nf90_def_dim(iFileID, 'bound',       27,   iBoundDim)
    iStatus = nf90_def_dim(iFileID, 'energy',      nE,  iEnDim)
    iStatus = nf90_def_dim(iFileID, 'pitch_angle', nPa, iPaDim)
    iStatus = nf90_def_dim(iFileID, 'time',  nf90_unlimited, iTimeDim)


    ! Define variables and attributes.  Attributes follow ncdf 
    ! conventions where possible (see unidata.ucar.edu for details.)
          ! TIME
    iStatus = nf90_def_var(iFileID, 'Time', nf90_float, iTimeDim, iTimeVar)
    iStatus = nf90_put_att(iFileID, iTimeVar, 'title', &
         'Time in seconds from start of simulation')
    iStatus = nf90_put_att(iFileID, iTimeVar, 'units', 'seconds')

          ! BAD DATA FLAG
    iStatus = nf90_def_var(iFileID, 'BadData', nf90_float, iFlagVar)
    iStatus = nf90_put_att(iFileID, iFlagVar, 'title', &
         'Value for missing data, typically used when sat outside domain.')

          ! DTWRITE
    iStatus = nf90_def_var(iFileID, 'DtWrite', nf90_float, iDtVar)
    iStatus = nf90_put_att(iFileID, iDtVar, 'title', &
         'Time in seconds between each attempted write.')
    iStatus = nf90_put_att(iFileID, iDtVar, 'disclaimer', &
         'Gaps may still exist!')

          ! LOCATION
    iStatus = nf90_def_var(iFileID, 'SM_xyz', nf90_float, &
         (/iXyzDim, iTimeDim/), iXyzVar)
    iStatus = nf90_put_att(iFileID, iXyzVar, 'title', &
         'Location of spacecraft in SM coordiantes.')
    iStatus = nf90_put_att(iFileID, iXyzVar, 'units', 'Earth Radii')

          ! TOTAL MAG FIELD
    iStatus = nf90_def_var(iFileID, 'B_xyz', nf90_float, &
         (/iXyzDim, iTimeDim/), iBVar)
    iStatus = nf90_put_att(iFileID, iBVar, 'title', &
         'Total Magnetic field at spacecraft in SM coordinates.')
    iStatus = nf90_put_att(iFileID, iBVar, 'units', 'nT')

          ! EXTERNAL MAG FIELD
    iStatus = nf90_def_var(iFileID, 'Bext_xyz', nf90_float, &
         (/iXyzDim, iTimeDim/), iBeVar)
    iStatus = nf90_put_att(iFileID, iBeVar, 'title', &
         'Total minus dipole field at spacecraft in SM coordinates.')
    iStatus = nf90_put_att(iFileID, iBeVar, 'units', 'nT')

          ! E-CONVECTION
    iStatus = nf90_def_var(iFileID, 'Econv_xyz', nf90_float, &
         (/iXyzDim, iTimeDim/), iEcVar)
    iStatus = nf90_put_att(iFileID, iEcVar, 'title', &
         'Convection Electric field at spacecraft in SM coordinates.')
    iStatus = nf90_put_att(iFileID, iEcVar, 'units', 'mV/m')

          ! FLUX
    do iS = 1, nS
       iStatus = nf90_def_var(iFileID, 'Flux'//species(iS)%s_code, nf90_float, &
                              (/iEnDim, iPaDim, iTimeDim/), iFluxVar)
       iStatus = nf90_put_att(iFileID, iFluxVar, 'title', &
            'Energy and pitch-angle dependent particle flux at spacecraft: '//species(iS)%s_code)
       iStatus = nf90_put_att(iFileID, iFluxVar, 'units', '1/cm2/s/ster/keV')
    enddo

         ! OMNIDIRECTIONAL FLUX
    do iS = 1, nS
       iStatus = nf90_def_var(iFileID, 'omni'//species(iS)%s_code, nf90_float, &
                              (/iEnDim, iTimeDim/), iOFluxVar)
       iStatus = nf90_put_att(iFileID, iOFluxVar, 'title', &
            'Energy and pitch-angle dependent particle flux at spacecraft: '//species(iS)%s_code)
       iStatus = nf90_put_att(iFileID, iOFluxVar, 'units', '1/cm2/s/keV')
    enddo

          ! ENERGY GRID
    iStatus = nf90_def_var(iFileID, 'energy_grid', nf90_float, &
         iEnDim, iEgridVar)
    iStatus = nf90_put_att(iFileID, iEgridVar, 'title', &
         'Energy grid as values at window centers.')
    iStatus = nf90_put_att(iFileID, iEgridVar, 'units', 'KeV')

    iStatus = nf90_def_var(iFileID, 'energy_width', nf90_float, &
         iEnDim, iEwidVar)
    iStatus = nf90_put_att(iFileID, iEwidVar, 'title', &
         'Widths of each energy bin whose centers are listed in energy_grid.')
    iStatus = nf90_put_att(iFileID, iEwidVar, 'units', 'KeV')

          ! PITCH-ANGLE GRID
    iStatus = nf90_def_var(iFileID, 'pa_grid', nf90_float, &
         iPaDim, iPgridVar)
    iStatus = nf90_put_att(iFileID, iPgridVar, 'title', &
         'Pitch angle grid as cosine of equitorial PA at window centers.')
    iStatus = nf90_put_att(iFileID, iPgridVar, 'units', 'unitless')

    iStatus = nf90_def_var(iFileID, 'pa_width', nf90_float, &
         iPaDim, iPwidVar)
    iStatus = nf90_put_att(iFileID, iPwidVar, 'title', &
         'Width of each pitch angle bin whose centers are listed in pa_grid.')
    iStatus = nf90_put_att(iFileID, iPwidVar, 'units', 'unitless')

    ! Write meta-data as Global Attributes.
    call write_ncdf_globatts(iFileID)

    ! Leave the "define mode" of the file.
    iStatus = nf90_enddef(iFileID)

    ! Write static (no time dimension) variables.
    iStatus = nf90_put_var(iFileID, iEgridVar, Ekev)
    iStatus = nf90_put_var(iFileID, iEwidVar,  wE)
    iStatus = nf90_put_var(iFileID, iPgridVar, Mu)
    iStatus = nf90_put_var(iFileID, iPwidVar,  wMu)
    iStatus = nf90_put_var(iFileID, iDtVar,    DtWriteSat)  
    iStatus = nf90_put_var(iFileID, iFlagVar,  BadDataFlag)

    ! Close the file.
    iStatus = nf90_close(iFileID)

  end subroutine create_sat_file

!============================================================================
  subroutine fly_sats
    ! Fly all virtual satellites; take virtual measurments.  Virtually fun.
    ! Locate virtual sats, interpolate result to position, and write data.

    use ModCoordTransform
    use ModConst,        ONLY: cPi
    use ModRamFunctions
    use ModRamTiming,    ONLY: TimeRamElapsed
    use ModRamGrids,     ONLY: nE, nPa, nS
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModRamVariables, ONLY: WMU
    use ModScbVariables, ONLY: x, y, z, bX, bY, bZ, bnormal, EXConv, EYConv, &
                               EZConv, bXIntern, bYIntern, bZIntern
    use ModRamParams, ONLY: DoSaveRamSats

    use ModRamScb, ONLY: flux3DEQ, indexPA

    use ModRamGSL, ONLY: GSL_NN
    use ModTimeConvert,  ONLY: time_real_to_int, TimeType

    implicit none
    integer :: GSLerr, i, iPa, iTime, iSat, iLoc(27), jLoc(27), kLoc(27), iTemp(3), &
               Pai, ix, ii, ij, ik, iS, iE, iT, iA(3)
    real(DP) :: xSat(3), dTime, xNear(27), yNear(27), zNear(27), BtNear(3,27), &
                BeNear(3,27), EcNear(3,27), xyzNear(3,27), &
                xNearT(27), yNearT(27), zNearT(27), &
                rNear, SatB(6), SatEc(3)
    real(DP), ALLOCATABLE :: distance(:,:,:), SatFlux(:,:,:), OmnFlux(:,:), &
                             SatFluxNear(:,:,:,:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = NameMod // '::fly_sats'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(.not. DoSaveRamSats) return
   
    ALLOCATE(distance(nthe,npsi,nzeta-1),SatFlux(nS,nE,nPa),OmnFlux(nS,nE), &
             SatFluxNear(nS,nE,nPa,27))
    distance = 0.0; SatFlux = 0.0; OmnFlux = 0.0; SatFluxNear = 0.0

    SATLOOP: do iSat=1, nRamSat
       distance = 0.0

       ! If current time is outside bounds of orbit file, do not trace.
       if(  (TimeRamElapsed .lt. SatTime_II(iSat, 1)) .or. &
            (TimeRamElapsed .ge. SatTime_II(iSat, nSatPoints_I(iSat))) ) then
          if(DoTest)then
             write(*,*) 'Time out of bounds for satellite', SatName_I(iSat)
             write(*,*) 'TimeRam, Max/Min SatTime = ', TimeRamElapsed, &
                  SatTime_II(iSat,1), SatTime_II(iSat, nSatPoints_I(iSat))
          end if
          cycle SATLOOP
       end if

       ! Find index of time nearest to current time.
       iTime = iSatTime_I(iSat)
       do while( (iTime < nSatPoints_I(iSat)) .and. &
            (SatTime_II(iSat, iTime) < TimeRamElapsed) )
          iTime = iTime + 1
       end do
       iSatTime_I(iSat) = iTime ! Save for next time

       if ((SatTime_II(iSat,iTime)-SatTime_II(iSat,iTime-1)).gt.600) then
          xSat = SatTraj_IID(iSat,iTime-1,:)
          SatB = BadDataFlag; SatEc = BadDataFlag
          SatFlux = BadDataFlag; OmnFlux = BadDataFlag
          xyzNear = BadDataFlag; BtNear = BadDataFlag
          call append_sat_record(SatFileName_O(iSat), iSatRecord(iSat), TimeRamElapsed, &
               xSat, SatB, SatEc, SatFlux, OmnFlux)
          iSatRecord(iSat) = iSatRecord(iSat) + 1
          cycle SATLOOP                 ! don't trace this time.
       else
       ! Set satellite's current position by interpolating to 
       ! current simulation time.
          dTime = (TimeRamElapsed - SatTime_II(iSat, iTime-1)) / &
                  (SatTime_II(iSat, iTime) - SatTime_II(iSat, iTime-1))
          xSat = dTime*(SatTraj_IID(iSat,iTime,:) - SatTraj_IID(iSat,iTime-1,:)) &
               + SatTraj_IID(iSat,iTime-1,:)
       end if
       ! Calculate distance to each point in 3DEQ domain.
       ! Note that we use 1:nzeta-1 as point zeta=1 == zeta=nzeta.
       distance = (x(:,:,1:nzeta-1) - xSat(1))**2 + &
                  (y(:,:,1:nzeta-1) - xSat(2))**2 + &
                  (z(:,:,1:nzeta-1) - xSat(3))**2
       ! Collect indices of nearest neighbor
       iTemp = minloc(distance)
       iLoc(:)=0; jLoc(:)=0; kLoc(:)=0
       xNear(:)=BadDataFlag
       yNear(:)=BadDataFlag
       zNear(:)=BadDataFlag
       xyzNear(:,:)=BadDataFlag
       BtNear(:,:)=BadDataFlag
       BeNear(:,:)=BadDataFlag
       EcNear(:,:)=BadDataFlag
       iT = 1
       iA(1) = 0; iA(2) = -1; iA(3) = 1
       do ii = 1,3
          do ij = 1,3
             do ik = 1,3
                iLoc(iT) = iTemp(1) + iA(ii)
                jLoc(iT) = iTemp(2) + iA(ij)
                kLoc(iT) = iTemp(3) + iA(ik)
                if (((iLoc(iT).gt.nthe).or.(iLoc(iT).lt.1)).or. &
                    ((jLoc(iT).gt.npsi).or.(jLoc(iT).lt.1)).or. &
                    ((kLoc(iT).gt.nzeta).or.(kLoc(iT).lt.1))) then
                   iT = iT
                   if (DoTest) then
                       write(*,*) 'Satellite ', SatName_I(iSat), ' out of bounds!'
                       write(*,*) 'Distance = ', minval(distance)
                       write(*,*) 'Sat location = ', xSat
                   end if
                else
                   xNear(iT) = x(iLoc(iT), jLoc(iT), kLoc(iT))
                   yNear(iT) = y(iLoc(iT), jLoc(iT), kLoc(iT))
                   zNear(iT) = z(iLoc(iT), jLoc(iT), kLoc(iT))
                   rNear = sqrt(xNear(iT)**2 + yNear(iT)**2 + zNear(iT)**2)
                   xyzNear(1,iT) = xNear(iT)
                   xyzNear(2,iT) = yNear(iT)
                   xyzNear(3,iT) = zNear(iT)
                   BtNear(1,iT) = bX(iLoc(iT), jLoc(iT), kLoc(iT))
                   BtNear(2,iT) = bY(iLoc(iT), jLoc(iT), kLoc(iT))
                   BtNear(3,iT) = bZ(iLoc(iT), jLoc(iT), kLoc(iT))
                   BeNear(1,iT) = BtNear(1,iT) - bxintern(iLoc(iT), jLoc(iT), kLoc(iT))
                   BeNear(2,iT) = BtNear(2,iT) - byintern(iLoc(iT), jLoc(iT), kLoc(iT))
                   BeNear(3,iT) = BtNear(3,iT) - bzintern(iLoc(iT), jLoc(iT), kLoc(iT))
                   EcNear(1,iT) = EXConv(iLoc(iT), jLoc(iT), kLoc(iT))
                   EcNear(2,iT) = EYConv(iLoc(iT), jLoc(iT), kLoc(iT))
                   EcNear(3,iT) = EZConv(iLoc(iT), jLoc(iT), kLoc(iT))
                   if (((abs(xNear(iT)).le.1e-9).and.&
                        (abs(yNear(iT)).le.1e-9).and.&
                        (abs(zNear(iT)).le.1e-9)) &
                       .or.(rNear.lt.1.75).or.(rNear.gt.6.75)) then
                      iT = iT
                   else
                      iT = iT + 1
                   endif
                end if
             end do
          end do
       end do
       iT = iT - 1

       if(DoTest) then
          write(*,*) 'Using these indices:'
          write(*,*) 'iLoc = ', iLoc
          write(*,*) 'jLoc = ', jLoc
          write(*,*) 'kLoc = ', kLoc
          
          write(*,*) 'Which gives us these locations (x,y,z):'
          write(*,'(4(f11.7,1x))') xNear
          write(*,'(4(f11.7,1x))') yNear
          write(*,'(4(f11.7,1x))') zNear
          
          write(*,*) 'And this B-field:'
          write(*,'(4(f11.7,1x))') BtNear(1,:)
          write(*,'(4(f11.7,1x))') BtNear(2,:)
          write(*,'(4(f11.7,1x))') BtNear(3,:)
       end if

       ! Interpolate all variables to this point.
       if (iT.le.3) then
        SatB = BadDataFlag
        SatEc = BadDataFlag
        OmnFlux = BadDataFlag
        SatFlux = BadDataFlag
       else
          DO i = 1,3
             ! Magnetic Field:
             CALL GSL_NN(xNear(1:iT),yNear(1:iT),zNear(1:iT),BtNear(i,1:iT), &
                         xSat(1),xSat(2),xSat(3),SatB(i),GSLerr)
             CALL GSL_NN(xNear(1:iT),yNear(1:iT),zNear(1:iT),BeNear(i,1:iT), &
                         xSat(1),xSat(2),xSat(3),SatB(i+3),GSLerr)
             ! Convective E:
             CALL GSL_NN(xNear(1:iT),yNear(1:iT),zNear(1:iT),EcNear(i,1:iT), &
                         xSat(1),xSat(2),xSat(3),SatEc(i),GSLerr)
          END DO
          SatB = SatB * bnormal ! Convert to correct units (nT)
          SatEc = SatEc * 1.0/6.4 ! Convert to correct units (mV/m)

          ! Reset Omnidirectional flux.
          OmnFlux(:,:)   = 0.0
          SatFlux(:,:,:) = BadDataFlag
          ! Flux for all energies, pitch angles, and species.
          do iS=1, nS
             do iE=1, nE
                do iPa=1, nPa
                   ix = 0
                   SatFluxNear(iS,iE,iPa,:) = 0.0
                   xNearT = 0.0; yNearT = 0.0; zNearT = 0.0
                   do i=1, iT
                      Pai = indexPA(iLoc(i), jLoc(i), kLoc(i), iPa)
                      if (Pai.gt.0) then
                         if (flux3DEQ(iS,jLoc(i),kLoc(i),iE,Pai).gt.0.0) then
                            ix = ix + 1
                            SatFluxNear(iS,iE,iPa,ix) = flux3DEQ(iS,jLoc(i),kLoc(i),iE,Pai)
                            xNearT(ix) = xNear(i)
                            yNearT(ix) = yNear(i)
                            zNearT(ix) = zNear(i)
                         end if
                      end if
                   end do
                   if (ix.gt.3) then
                      CALL GSL_NN(xNearT(1:ix),yNearT(1:ix),zNearT(1:ix),SatFluxNear(iS,iE,iPa,1:ix), &
                                  xSat(1),xSat(2),xSat(3),SatFlux(iS,iE,iPa),GSLerr)
                      if (SatFlux(iS,iE,iPa).gt.0.0) then
                         OmnFlux(iS,iE) = OmnFlux(iS,iE) + SatFlux(iS,iE,iPa)*4.0*cPi*wMu(iPa)
                      end if
                   end if
                end do
                if (OmnFlux(iS,iE).le.0.0) OmnFlux(iS,iE) = BadDataFlag
             end do
          end do
       end if
 
       ! Write information to output file.
       call append_sat_record(SatFileName_O(iSat), iSatRecord(iSat), TimeRamElapsed, &
            xSat, SatB, SatEc, SatFlux, OmnFlux) !, SatEi

       ! Increment record location in NCDF file.
       iSatRecord(iSat) = iSatRecord(iSat) + 1

    end do SATLOOP

    DEALLOCATE(distance,SatFlux,OmnFlux,SatFluxNear)

  contains

    !========================================================================
    subroutine append_sat_record(InFileName, iRec, TimeIn, xVec, bVec, &
                                 ecVec, FluxIn, OmnFluxIn)! eiVec
      !Open a netCDF file and write new values.

      use netcdf
      use ModRamGrids, ONLY: nS, nE, nPa
      use ModRamVariables, ONLY: species

      implicit none

      ! Arguments:
      integer, intent(in)           :: iRec
      character(len=200), intent(in):: InFileName
      real(DP), intent(in) :: xVec(:), bVec(:), ecVec(:), TimeIn, FluxIn(:,:,:), &
                              OmnFluxIn(:,:) !, eiVec
      
      ! Local variables:
      real(DP), ALLOCATABLE :: Flux(:,:,:), OmFx(:,:)
      integer :: iS
      integer :: iStatus, iFileID, iStart1D(1), iStart2D(2), iStart3D(3)
      integer :: iTimeVar, iXyzVar, iBVar, iFluxVar, iOmniVar, iEcVar, iBeVar

      logical :: DoTest, DoTestMe
      character(len=100) :: NameSubSub = NameSub // '::append_sat_record'
      !-------------------------------------------------------------------
      call CON_set_do_test(NameSubSub, DoTest, DoTestMe)
      ALLOCATE(Flux(nS,nE,nPa),OmFx(nS,nE))
      Flux = 0.0; OmFx = 0.0

      ! Copy flux to a local variable.
      Flux = FluxIn
      OmFx = OmnFluxIn

      ! Set starting location matrices.
      iStart1D = (/iRec/)
      iStart2D = (/1, iRec/)
      iStart3D = (/1,1,iRec/)

      if(DoTest) then
         write(*,*) 'Dumping data to ', trim(InFileName)
         write(*,*) 'Time = ', TimeIn
         write(*,*) 'XYZ=', xVec
         write(*,*) 'B  =', bVec
      end if

      ! Open the NetCDF file.
      iStatus = nf90_open(trim(InFileName), nf90_write, iFileID)

      ! Collect variable IDs.
      iStatus = nf90_inq_varid(iFileID, 'Time',      iTimeVar)
      iStatus = nf90_inq_varid(iFileID, 'SM_xyz',    iXyzVar)
      iStatus = nf90_inq_varid(iFileID, 'B_xyz',     iBVar)
      iStatus = nf90_inq_varid(iFileID, 'Bext_xyz',  iBeVar)
      iStatus = nf90_inq_varid(iFileID, 'Econv_xyz', iEcVar)

      ! Write new values to file:
           ! Time
      iStatus = nf90_put_var(iFileID, iTimeVar, TimeIn, iStart1D)
           ! Satelite position
      iStatus = nf90_put_var(iFileID, iXyzVar, xVec, iStart2D)
           ! Vector magnetic field.
      iStatus = nf90_put_var(iFileID, iBVar, bVec(1:3), iStart2D)
      iStatus = nf90_put_var(iFileID, iBeVar,bVec(4:6), iStart2D)
           ! Vector convection E-field
      iStatus = nf90_put_var(iFileID, iEcVar, ecVec, iStart2D)


      do iS = 1, nS
         iStatus = nf90_inq_varid(iFileID, 'Flux'//species(iS)%s_code, iFluxVar)
         iStatus = nf90_inq_varid(iFileID, 'omni'//species(iS)%s_code, iOmniVar)
         iStatus = nf90_put_var(iFileID, iFluxVar, Flux(iS,:,:), iStart3D)
         iStatus = nf90_put_var(iFileID, iOmniVar, OmFx(iS,:),   iStart2D)
      enddo
      
      ! Close NCDF file.
      iStatus = nf90_close(iFileID)

      DEALLOCATE(Flux,OmFx)
    end subroutine append_sat_record
    
    !========================================================================
  end subroutine fly_sats

!============================================================================
end module ModRamSats

!============================================================================
