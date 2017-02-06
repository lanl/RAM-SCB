!============================================================================
module ModRamSats
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.

  use ModRamMain, ONLY: Real8_

  implicit none
  save
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
  character(len=3)   :: TypeSatCoord_I(MaxRamSat)

  ! Info about the sat output netcdf files.
  integer :: iSatRecord(MaxRamSat) = 1  ! Location in netCDF file.
  
  ! Variables that contain the time and location of the sats.
  real(kind=Real8_), allocatable :: SatTraj_IID(:,:,:), SatTime_II(:,:)

  ! Bad data flag (sat outside domain.)
  real(kind=Real8_)  :: BadDataFlag=-1.0E10

  character(len=3), parameter :: TypeCoordSystem = 'SMG'

  character(len=*), parameter :: NameMod = 'ModRamSats'

contains

!============================================================================
  subroutine read_sat_params

    use ModReadParam
    use ModRamIO, ONLY: DtWriteSat, DoSaveRamSats

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
    use ModRamMain,     ONLY: TimeRamStart
    use ModTimeConvert, ONLY: time_int_to_real
    use ModIoUnit,      ONLY: UnitTmp_
    use ModRamMpi,      ONLY: iProc
    use CON_axes

    integer :: iError, i, iSat , nPoint

    ! One line of input
    character (len=100) :: line

    integer      :: iTime_I(7)
    integer      :: MaxPoint
    real(kind=Real8_) :: Xyz_D(3)
    real(kind=Real8_) :: DateTime
    real(kind=Real8_), allocatable :: Time_I(:), Xyz_DI(:,:)
    character(len=100):: NameFile

    character(len=*), parameter :: NameSub = NameMod // '::read_sat_input'

    logical :: DoTest, DoTestMe

    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
    if(iProc == 0)then ! Run only on head node.
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

    end if ! End running on head node only.

    ! Tell all processors the maximum number of points
    !call MPI_Bcast(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(3, MaxPoint))
    allocate(SatTraj_IID(nRamSat, MaxPoint, 3))
    allocate(SatTime_II(nRamSat, MaxPoint))

    ! Read the trajectories
    SATELLITES: do iSat=1, nRamSat

       ! Read file on the root processor
       if (iProc == 0) then

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

       end if

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

    use ModRamMain, ONLY: IsRestart, PathRamOut

    character(len=200) :: FileName

    integer :: i

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = NameMod // '::init_sat_files'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)


    do i=1, nRamSat
       ! Replace now-useless input file name with output file name.
       SatFileName_I(i) = trim(PathRamOut)//'/'//trim(SatName_I(i))//'.nc'
       if(IsRestart) then
          ! On restart, attempt to continue existing output file.
          call resume_sat_file(SatFileName_I(i), iSatRecord(i))
       else
          ! else create a new file from scratch.
          call create_sat_file(SatFileName_I(i))
       end if
    end do

  end subroutine init_sats

!============================================================================
  subroutine create_sat_file(FileNameIn)

    ! Build a new netCDF satellite output file.  Enter all meta data into
    ! the file.

    use netcdf
    use ModTimeConvert
    use ModRamMain, ONLY: nE, nPa, Ekev, wE, Mu, wMu, TimeRamStart
    use ModRamIO,   ONLY: ncdf_check, write_ncdf_globatts, DtWriteSat

    character(len=200), intent(in) :: FileNameIn

    integer :: i, iFileID, iStatus, iTimeDim, iEnDim, iPaDim, iXyzDim
    integer :: iTimeVar, iXyzVar, iBVar, iEgridVar, iEwidVar, iPgridVar, iBeVar
    integer :: iHVar, iHeVar, iOVar, ieVar, iStartVar, iDtVar, iEcVar, &
         iEiVar, iPwidVar, ioHVar, ioHeVar, ioOVar, ioeVar, iFlagVar

    logical :: DoTest, DoTestMe
    character(len=100) :: NameSub = NameMod // '::create_sat_file'
    !-----------------------------------------------------------------------
    call CON_set_do_test(trim(NameSub), DoTest, DoTestMe)
    
    if(DoTest) write(*,*) trim(NameSub), &
         ' Creating NetCDF File ', trim(FileNameIn)

    ! Open file; check for success.
    iStatus = nf90_create(FileNameIn, nf90_clobber, iFileID)
    call ncdf_check(iStatus, NameSub)


    ! Create dimensions for all variables.
    iStatus = nf90_def_dim(iFileID, 'xyz',         3,   iXyzDim)
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

          ! E-INDUCED
    !iStatus = nf90_def_var(iFileID, 'Eind_xyz', nf90_float, &
    !     (/iXyzDim, iTimeDim/), iEcVar)
    !iStatus = nf90_put_att(iFileID, iEiVar, 'title', &
    !     'Induced Electric field at spacecraft in SM coordinates.')
    !iStatus = nf90_put_att(iFileID, iEiVar, 'units', 'mV/m')

          ! FLUX
    iStatus = nf90_def_var(iFileID, 'FluxH+', nf90_float, &
         (/iEnDim, iPaDim, iTimeDim/), iHVar)
    iStatus = nf90_put_att(iFileID, iHVar, 'title', &
         'Energy and pitch-angle dependent particle flux at spacecraft: H+.')
    iStatus = nf90_put_att(iFileID, iHVar, 'units', '1/cm2/s/ster/keV')

    iStatus = nf90_def_var(iFileID, 'FluxHe+', nf90_float, &
         (/iEnDim, iPaDim, iTimeDim/), iHeVar)
    iStatus = nf90_put_att(iFileID, iHeVar, 'title', &
         'Energy and pitch-angle dependent particle flux at spacecraft: He+.')
    iStatus = nf90_put_att(iFileID, iHeVar, 'units', '1/cm2/s/ster/keV')

    iStatus = nf90_def_var(iFileID, 'FluxO+', nf90_float, &
         (/iEnDim, iPaDim, iTimeDim/), iOVar)
    iStatus = nf90_put_att(iFileID, iOVar, 'title', &
         'Energy and pitch-angle dependent particle flux at spacecraft: O+.')
    iStatus = nf90_put_att(iFileID, iOVar, 'units', '1/cm2/s/ster/keV')

    iStatus = nf90_def_var(iFileID, 'Fluxe-', nf90_float, &
         (/iEnDim, iPaDim, iTimeDim/), ieVar)
    iStatus = nf90_put_att(iFileID, ieVar, 'title', &
         'Energy and pitch-angle dependent particle flux at spacecraft: e-.')
    iStatus = nf90_put_att(iFileID, ieVar, 'units', '1/cm2/s/ster/keV')

         ! OMNIDIRECTIONAL FLUX
    iStatus = nf90_def_var(iFileID, 'omniH', nf90_float, &
         (/iEnDim, iTimeDim/), ioHVar)
    iStatus = nf90_put_att(iFileID, ioHVar, 'title', &
         'Omnidirectional particle flux at spacecraft: H+.')
    iStatus = nf90_put_att(iFileID, ioHVar, 'units', '1/cm2/s/keV')

    iStatus = nf90_def_var(iFileID, 'omniO', nf90_float, &
         (/iEnDim, iTimeDim/), ioOVar)
    iStatus = nf90_put_att(iFileID, ioOVar, 'title', &
         'Omnidirectional particle flux at spacecraft: O+.')
    iStatus = nf90_put_att(iFileID, ioOVar, 'units', '1/cm2/s/keV')

    iStatus = nf90_def_var(iFileID, 'omniHe', nf90_float, &
         (/iEnDim, iTimeDim/), ioHeVar)
    iStatus = nf90_put_att(iFileID, ioHeVar, 'title', &
         'Omnidirectional particle flux at spacecraft: He+.')
    iStatus = nf90_put_att(iFileID, ioHeVar, 'units', '1/cm2/s/keV')

    iStatus = nf90_def_var(iFileID, 'omnie', nf90_float, &
         (/iEnDim, iTimeDim/), ioeVar)
    iStatus = nf90_put_att(iFileID, ioeVar, 'title', &
         'Omnidirectional particle flux at spacecraft: e-.')
    iStatus = nf90_put_att(iFileID, ioeVar, 'units', '1/cm2/s/keV')

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
    call ncdf_check(iStatus, NameSub)

    ! Write static (no time dimension) variables.
    iStatus = nf90_put_var(iFileID, iEgridVar, Ekev)
    iStatus = nf90_put_var(iFileID, iEwidVar,  wE)
    iStatus = nf90_put_var(iFileID, iPgridVar, Mu)
    iStatus = nf90_put_var(iFileID, iPwidVar,  wMu)
    iStatus = nf90_put_var(iFileID, iDtVar,    DtWriteSat)  
    iStatus = nf90_put_var(iFileID, iFlagVar,  BadDataFlag)
    !call time_int_to_real(TimeRamStart)
    !iStatus = nf90_put_var(iFileID, iStartVar, TimeRamStart % String)
    call ncdf_check(iStatus, NameSub)

    ! Close the file.
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine create_sat_file

!============================================================================
  subroutine resume_sat_file(FileNameIn, iRecordOut)

    ! If restarting, it is possible to resume writing to any netCDF
    ! virtual satellite output file.  This subroutine examines a
    ! virtual satellite output file, determines if it matches the 
    ! current format, locates the record at which to continue writing, 
    ! and saves this information so that writing may continue as this
    ! simulation continues.

    ! The logical path this subroutine takes is as follows:
    ! 1) check to see if FileNameIn exists.
    ! 2) check if FileNameIn matches the format (size, variables) required.
    ! If either of these fail, a new output file is created.
    ! If both succeed, the proper location to begin adding new records
    ! is returned as iRecordOut.

    character(len=200), intent(in) :: FileNameIn
    integer, intent(out)           :: iRecordOut

    character(len=*), parameter :: NameSub = NameMod//'::resume_sat_file'
    !-----------------------------------------------------------------------

    call CON_stop(NameSub // ' This subroutine is not ready for use.')
    return

  end subroutine resume_sat_file

!============================================================================
  subroutine fly_sats
    ! Fly all virtual satellites; take virtual measurments.  Virtually fun.
    ! Locate virtual sats, interpolate result to position, and write data.

    use ModConst,   ONLY: cPi
    use ModRamFunctions
    use ModRamMain, ONLY: TimeRamElapsed, PathRamOut, nE, nPa, wMu
    use ModRamIO,   ONLY: DoSaveRamSats
    use Module_RAM, ONLY: flux3DEQ, indexPA
    use Module1,    ONLY: x,y,z, nthe, npsi, nzeta, bX, bY, bZ, bnormal, &
         EXConv, EYConv, EZConv, bxintern, byintern, bzintern!, &
        !EXIndEq, EYIndEq, EZIndEq, enormal

    integer :: i, iPa, iTime, iSat, iLoc(4), jLoc(4), kLoc(4), iTemp(3)
    real(kind=Real8_) :: xSat(3), dTime, distance(nthe, npsi, nzeta-1), &
         xNear(4), yNear(4), zNear(4), BtNear(3,4), BeNear(3,4), &
         EcNear(3,4)!, EiNear(3,4)
    real(kind=Real8_), parameter :: MaxDist = 0.25

    character(len=200) :: FileName

    ! Buffers to write to file:
    real(kind=Real8_) :: SatB(6), SatEc(3), SatEi(3), &
         SatFlux(4, nE, nPa)=0.0, OmnFlux(4, nE)=0.0

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = NameMod // '::fly_sats'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(.not. DoSaveRamSats) return

    SATLOOP: do iSat=1, nRamSat
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


       ! Set satellite's current position by interpolating to 
       ! current simulation time.
       dTime = (TimeRamElapsed - SatTime_II(iSat, iTime-1)) / &
            (SatTime_II(iSat, iTime) - SatTime_II(iSat, iTime-1))
       xSat = dTime*(SatTraj_IID(iSat,iTime,:) - SatTraj_IID(iSat,iTime-1,:)) &
            + SatTraj_IID(iSat,iTime-1,:)

       ! Calculate distance to each point in 3DEQ domain.
       ! Note that we use 2:nzeta as point zeta=1 == zeta=nzeta.
       distance = (x(:,:,2:nzeta) - xSat(1))**2 + &
            (y(:,:,2:nzeta) - xSat(2))**2 + &
            (z(:,:,2:nzeta) - xSat(3))**2
       
       ! Collect indices of four nearest neighbors.
       ! Also collect data values at those points.
       do i=1, 4
          iTemp = minloc(distance)           ! find closest point.
          if(minval(distance) .gt. MaxDist**2) then !&! if outside domain...
             if(DoTest) then
                write(*,*) 'Satellite ', SatName_I(iSat), ' out of bounds!'
                write(*,*) 'Distance = ', minval(distance)
                write(*,*) 'Sat location = ', xSat
             end if

             ! Fill in with blank values.
             SatB=BadDataFlag; SatEc=BadDataFlag; SatEi=BadDataFlag
             SatFlux=BadDataFlag; OmnFlux=BadDataFlag

             ! Write information to output file.
             FileName = trim(PathRamOut)//'/'//trim(SatName_I(iSat))//'.nc'
             call append_sat_record(FileName, iSatRecord(iSat), &
                  TimeRamElapsed, xSat, SatB, SatEc, SatFlux, OmnFlux) !, SatEi

             ! Increment record location in NCDF file.
             iSatRecord(iSat) = iSatRecord(iSat) + 1

             cycle SATLOOP                 ! don't trace this time.
          end if
          iLoc(i) = iTemp(1)                 !\
          jLoc(i) = iTemp(2)                 ! Save i,j,k indices.
          kLoc(i) = iTemp(3)                 !/

          ! Exclude this point in next iteration of this loop.
          distance(iLoc(i),jLoc(i),kLoc(i)) = 99.0 
          
          !Collect position and B-field values at NNs.
          xNear(i) = x(iTemp(1), iTemp(2), iTemp(3))
          yNear(i) = y(iTemp(1), iTemp(2), iTemp(3))
          zNear(i) = z(iTemp(1), iTemp(2), iTemp(3))
          BtNear(1,i) = bX(iTemp(1), iTemp(2), iTemp(3))
          BtNear(2,i) = bY(iTemp(1), iTemp(2), iTemp(3))
          BtNear(3,i) = bZ(iTemp(1), iTemp(2), iTemp(3))
          BeNear(1,i) = BtNear(1,i) - bxintern(iTemp(1), iTemp(2), iTemp(3))
          BeNear(2,i) = BtNear(2,i) - byintern(iTemp(1), iTemp(2), iTemp(3))
          BeNear(3,i) = BtNear(3,i) - bzintern(iTemp(1), iTemp(2), iTemp(3))
          EcNear(1,i) = EXConv(iTemp(1), iTemp(2), iTemp(3))
          EcNear(2,i) = EYConv(iTemp(1), iTemp(2), iTemp(3))
          EcNear(3,i) = EZConv(iTemp(1), iTemp(2), iTemp(3))
          !EiNear(1,i = EXIndEq(iTemp(1), iTemp(2), iTemp(3))
          !EiNear(2,i = EYIndEq(iTemp(1), iTemp(2), iTemp(3))
          !EiNear(3,i = EZIndEq(iTemp(1), iTemp(2), iTemp(3))
       end do

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
       ! Magnetic field:
       SatB(1) = WeightFit(4, xNear, yNear, zNear, BtNear(1,:), xSat)
       SatB(2) = WeightFit(4, xNear, yNear, zNear, BtNear(2,:), xSat)
       SatB(3) = WeightFit(4, xNear, yNear, zNear, BtNear(3,:), xSat)
       SatB(4) = WeightFit(4, xNear, yNear, zNear, BeNear(1,:), xSat)
       SatB(5) = WeightFit(4, xNear, yNear, zNear, BeNear(2,:), xSat)
       SatB(6) = WeightFit(4, xNear, yNear, zNear, BeNear(3,:), xSat)
       SatB = SatB * bnormal ! Convert to correct units (nT)

       ! Convective E:
       SatEc(1) = WeightFit(4, xNear, yNear, zNear, EcNear(1,:), xSat)
       SatEc(2) = WeightFit(4, xNear, yNear, zNear, EcNear(2,:), xSat)
       SatEc(3) = WeightFit(4, xNear, yNear, zNear, EcNear(3,:), xSat)
       SatEc = SatEc * 1.0/6.4 ! Convert to correct units (mV/m)

       ! Induced E:
       !SatEi(1) = WeightFit(4, xNear, yNear, zNear, EiNear(1,:), xSat)
       !SatEi(2) = WeightFit(4, xNear, yNear, zNear, EiNear(2,:), xSat)
       !SatEi(3) = WeightFit(4, xNear, yNear, zNear, EiNear(3,:), xSat)
       !SatEi = SatEi * enormal ! Convert to correct units (mV/m)


       ! Reset Omnidirectional flux.
       OmnFlux(:,:)=0.0
       ! Flux for all energies, pitch angles, and species.
       ! Currently, we'll just use nearest neighbor.
       do iPa=1, nPa
          if(indexPA(iLoc(1), jLoc(1), kLoc(1), iPa) .eq. -1)then
             SatFlux(:,:,iPa) = -1.0
          else
             SatFlux(:,:,iPa) = flux3DEQ(:, jLoc(1), kLoc(1), :, &
                  indexPA(iLoc(1), jLoc(1), kLoc(1), iPa))
             OmnFlux(:,:) = OmnFlux(:,:) + SatFlux(:,:,iPa)*4.0*cPi*wMu(iPa)
          end if
       end do

       ! Write information to output file.
       FileName = trim(PathRamOut)//'/'//trim(SatName_I(iSat))//'.nc'
       call append_sat_record(FileName, iSatRecord(iSat), TimeRamElapsed, &
            xSat, SatB, SatEc, SatFlux, OmnFlux) !, SatEi

       ! Increment record location in NCDF file.
       iSatRecord(iSat) = iSatRecord(iSat) + 1

    end do SATLOOP

  contains

    !========================================================================
    subroutine append_sat_record(InFileName, iRec, TimeIn, xVec, bVec, &
         ecVec, FluxIn, OmnFluxIn)! eiVec
      !Open a netCDF file and write new values.

      use netcdf
      use ModRamIO,   ONLY: ncdf_check
      use ModRamMain, ONLY: nE, nPa

      ! Arguments:
      integer, intent(in)           :: iRec
      character(len=200), intent(in):: InFileName
      real(kind=Real8_), intent(in) :: xVec(3), bVec(6), ecVec(3)
      real(kind=Real8_), intent(in) :: TimeIn, FluxIn(4,nE,nPa), &
           OmnFluxIn(4,nE) !, eiVec
      
      ! Local variables:
      real(kind=Real8_) :: Flux(4,nE,nPa), OmFx(4,nE)
      integer :: iStatus, iFileID, iStart1D(1), iStart2D(2), iStart3D(3)
      integer :: iTimeVar, iXyzVar, iBVar, iHVar, iHeVar, iOVar, ieVar, &
           iEcVar, iBeVar, ioHVar, ioHeVar, ioOVar, ioeVar !iEiVar

      logical :: DoTest, DoTestMe
      character(len=100) :: NameSubSub = NameSub // '::append_sat_record'
      !-------------------------------------------------------------------
      call CON_set_do_test(NameSubSub, DoTest, DoTestMe)

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
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))

      ! Collect variable IDs.
      iStatus = nf90_inq_varid(iFileID, 'Time',      iTimeVar)
      iStatus = nf90_inq_varid(iFileID, 'SM_xyz',    iXyzVar)
      iStatus = nf90_inq_varid(iFileID, 'B_xyz',     iBVar)
      iStatus = nf90_inq_varid(iFileID, 'Bext_xyz',  iBeVar)
      iStatus = nf90_inq_varid(iFileID, 'FluxH+',    iHVar)
      iStatus = nf90_inq_varid(iFileID, 'FluxHe+',   iHEVar)
      iStatus = nf90_inq_varid(iFileID, 'FluxO+',    iOVar)
      iStatus = nf90_inq_varid(iFileID, 'Fluxe-',    ieVar)
      iStatus = nf90_inq_varid(iFileID, 'Econv_xyz', iEcVar)
      iStatus = nf90_inq_varid(iFileID, 'omniH', ioHVar)
      iStatus = nf90_inq_varid(iFileID, 'omniHe',ioHeVar)
      iStatus = nf90_inq_varid(iFileID, 'omniO', ioOVar)
      iStatus = nf90_inq_varid(iFileID, 'omnie', ioeVar)
      !iStatus = nf90_inq_varid(iFileID, 'Eind_xyz',  iEiVar)

      ! Write new values to file:
           ! Time
      iStatus = nf90_put_var(iFileID, iTimeVar, TimeIn, iStart1D)
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
           ! Satelite position
      iStatus = nf90_put_var(iFileID, iXyzVar, xVec, iStart2D)
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
           ! Vector magnetic field.
      iStatus = nf90_put_var(iFileID, iBVar, bVec(1:3), iStart2D)
      iStatus = nf90_put_var(iFileID, iBeVar,bVec(4:6), iStart2D)
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
           ! Vector convection E-field
      iStatus = nf90_put_var(iFileID, iEcVar, ecVec, iStart2D)
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
           ! Vector induced E-field
      !iStatus = nf90_put_var(iFileID, iEiVar, eiVec, iStart2D)
      !call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
           ! Flux for all species.
      where(Flux+1.0 .eq. Flux)
         Flux = -1.0
      end where
      iStatus = nf90_put_var(iFileID, ieVar,  Flux(1,:,:), iStart3D)
      iStatus = nf90_put_var(iFileID, iHVar,  Flux(2,:,:), iStart3D)
      iStatus = nf90_put_var(iFileID, iHeVar, Flux(3,:,:), iStart3D)
      iStatus = nf90_put_var(iFileID, iOVar,  Flux(4,:,:), iStart3D)
      iStatus = nf90_put_var(iFileID, ioeVar, OmFx(1,:), iStart2D)
      iStatus = nf90_put_var(iFileID, ioHVar, OmFx(2,:), iStart2D)
      iStatus = nf90_put_var(iFileID, ioHeVar,OmFx(3,:), iStart2D)
      iStatus = nf90_put_var(iFileID, ioOVar, OmFx(4,:), iStart2D)
      call ncdf_check(iStatus, NameSubSub, NameFileIn=trim(InFileName))
      
      ! Close NCDF file.
      iStatus = nf90_close(iFileID)

    end subroutine append_sat_record
    
    !========================================================================
  end subroutine fly_sats

!============================================================================
end module ModRamSats

!============================================================================
