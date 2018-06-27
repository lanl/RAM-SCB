!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamIO

  implicit none

  logical :: IsFramework = .false. 

  ! File output names and units
  integer :: iUnitLog ! Logfile IO unit
  character(len=214) :: NameFileLog, NameFileDst

  ! String that contains the date and time for which this simulation
  ! was initialized (used for output metadata.)
  character(len=21) :: StringRunDate
  
  character(len=10), parameter :: NameMod = 'ModRamMain'
  
  contains
!============================!
!===== BASE SUBROUTINES =====!
!============================!
  subroutine write_prefix

    use ModRamParams, ONLY: IsComponent

    implicit none

    character(len=7) :: StringPrefix = 'IM:'

    if(.not. IsComponent) RETURN
    write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix
!==============================================================================
  function RamFileName(PrefixIn, SuffixIn, TimeIn)
    ! Create a file name using the new RAM-SCB I/O filename standards:
    ! FileNameOut = PrefixIn_dYYYYMMDD_tHHMMSS.SuffixIn
    ! PrefixIn and SuffixIn are short, descriptive strings (e.g. 
    ! 'pressure', 'ram_o' for prefixes and 'dat', 'cdf' for suffixes.

    use ModTimeConvert, ONLY: TimeType

    implicit none

    character(len=200)           :: RamFileName
    type(TimeType),   intent(in) :: TimeIn
    character(len=*), intent(in) :: PrefixIn, SuffixIn
    !------------------------------------------------------------------------
    write(RamFileName, "(a,'_d',i4.4,2i2.2,'_t',3i2.2,'.',a)") &
         trim(PrefixIn), &
         TimeIn % iYear, &
         TimeIn % iMonth, &
         TimeIn % iDay, &
         TimeIn % iHour, &
         TimeIn % iMinute, &
         TimeIn % iSecond, &
         trim(SuffixIn)

  end function RamFileName

!============================================================================
  subroutine init_output
    ! Initialize any output that requires such action.

    use ModRamTiming, ONLY: TimeRamRealStart
    use ModRamSats,   ONLY: read_sat_input, init_sats
    use ModRamMain,   ONLY: PathRamOut, PathScbOut
    use ModRamParams, ONLY: DoSaveRamSats

    use ModIoUnit,    ONLY: UNITTMP_

    implicit none
    !------------------------------------------------------------------------
    if (DoSaveRamSats) then
       call read_sat_input
       call init_sats
    end if

    ! Initialize logfile.
    NameFileLog = trim(PathRamOut)//RamFileName('log','log',TimeRamRealStart)
    open(UNITTMP_, FILE=trim(NameFileLog), STATUS='REPLACE')
    write(UNITTMP_,*) 'RAM-SCB Log'
    write(UNITTMP_,*) 'time year mo dy hr mn sc msc dstRam dstBiot ', &
                      'pparh pperh pparo ppero pparhe pperhe ppare ppere'
    close(UNITTMP_)

    ! Initialize Dstfile.
    NameFileDst = trim(PathScbOut)//RamFileName('dst_computed_3deq','txt',TimeRamRealStart)
    open(UNITTMP_, FILE=trim(NameFileDst), STATUS='REPLACE')
    write(UNITTMP_,*) 'SCB Dst File'
    write(UNITTMP_,*) 'time year mo dy hr mn sc msc dstDPS dstDPSGeo ', &
                      'dstBiot dstBiotGeo'
    close(UNITTMP_)

  end subroutine init_output

!==============================================================================
  subroutine handle_output(TimeIn)
    ! For every type of output, check the following:
    ! 1) if that output is activated,
    ! 2) if it's time to write that output
    ! If so, call the proper routines to write the output.

    !!!! Module Variables
    use ModRamMain,      ONLY: DP, S
    use ModRamParams,    ONLY: DoSaveRamSats, NameBoundMag
    use ModRamTiming,    ONLY: Dt_hI, DtRestart, DtWriteSat, TimeRamNow, &
                               TimeRamElapsed, DtW_Pressure, DtW_hI, DtW_Efield, &
                               DtW_MAGxyz, DtLogFile
    use ModRamGrids,     ONLY: NR, NT
    use ModRamVariables, ONLY: PParH, PPerH, PParHe, PPerHe, PParE, PPerE, &
                               PParO, PPerO, VT
    use ModScbVariables, ONLY: DstBiot, DstBiotInsideGeo, DstDPS, DstDPSInsideGeo
    !!!! Module Subroutines/Functions
    use ModRamSats,      ONLY: fly_sats
    use ModRamFunctions, ONLY: ram_sum_pressure, get_ramdst
    use ModRamRestart,   ONLY: write_restart
    use ModScbIO,        ONLY: Write_MAGxyz
    ! Share Modules
    use ModIOUnit, ONLY: UNITTMP_

    implicit none

    real(DP), intent(in) :: TimeIn

    real(DP) :: dst
    logical :: DoTest, DoTestMe
    integer :: iS
    character(len=4) :: StrScbIter
    character(len=*), parameter :: NameSub = 'handle_output'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Write Logfile
    if(abs(mod(TimeIn, DtLogfile))<=1e-9)then
       ! Get current Dst.
       call get_ramdst(dst)
       open(UNITTMP_, FILE=NameFileLog, POSITION='APPEND')
       write(UNITTMP_, '(E13.6,1x,i4,5(1x,i2.2),1x,i3.3,10(1x,E13.6))') &
            TimeIn, TimeRamNow%iYear, TimeRamNow%iMonth, &
            TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
            TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0), &
            dst, DstBiot, sum(PParH)/(nR*nT), sum(PPerH)/(nR*nT), &
            sum(PParO)/(nR*nT), sum(PPerO)/(nR*nT), sum(PParHe)/(nR*nT), &
            sum(PPerHe)/(nR*nT), sum(PParE)/(nR*nT), sum(PPerE)/(nR*nT)
       close(UNITTMP_)
    end if

    ! Write SCB Dst file
    if(abs(mod(TimeIn, DtLogfile))<=1e-9)then
       open(UNITTMP_, FILE=NameFileDst, POSITION='APPEND')
       write(UNITTMP_, *) TimeIn, TimeRamNow%iYear, TimeRamNow%iMonth, &
                          TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
                          TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0), &
                          dstDPS, dstDPSInsideGeo, DstBiot, DstBiotInsideGeo, &
                          maxval(VT), minval(VT)
       close(UNITTMP_)
    end if

    ! Write Pressure File
    if (abs(mod(TimeIn, DtW_Pressure))<=1e-9) then
       write(StrScbIter,'(I4.4)') int(TimeRamElapsed/Dt_hI)
       call ram_sum_pressure
       call ram_write_pressure
    end if

    ! Write MAGxyz File
    if ((abs(mod(TimeIn, DtW_MAGxyz))<=1e-9).and.(NameBoundMag.ne.'DIPL')) call Write_MAGxyz

    ! Write hI File
    if ((abs(mod(TimeIn, DtW_hI))<=1e-9).and.(NameBoundMag.ne.'DIPL')) call ram_write_hI

    ! Write efield file
    if (abs(mod(TimeIn, DtW_Efield)).le.1e-9) call ram_epot_write

    ! Update Satellite Files
    if (DoSaveRamSats .and. (abs(mod(TimeRamElapsed,DtWriteSat)) .le. 1e-9)) call fly_sats

    ! Write restarts.
    if (abs(mod(TimeIn,DtRestart)).le.1e-9) call write_restart()

    ! Write hourly file
    if (abs(mod(TimeIn,3600.0)).le.1e-9) then
       do iS = 1,4
          S = iS
          call ram_hour_write
       enddo
    endif

  end subroutine handle_output

!===========================!
!==== INPUT SUBROUTINES ====!
!===========================!
!==============================================================================
subroutine read_geomlt_file(NameParticle)

  !!!! Module Variables
  use ModRamMain,      ONLY: DP
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamStart
  use ModRamGrids,     ONLY: NEL, NEL_prot, NTL, NBD
  use ModRamParams,    ONLY: boundary, BoundaryPath
  use ModRamVariables, ONLY: flux_SIII, fluxLast_SII, eGrid_SI, &
                             avgSats_SI, lGrid_SI, tGrid_SI
  !!!! Share Modules
  use ModIOUnit,      ONLY: UNITTMP_
  use ModTimeConvert, ONLY: TimeType, time_int_to_real

  implicit none

  character(len=4), intent(in) :: NameParticle

  character(len=200) :: NameFileIn, StringHeader
  integer :: iError, i, j, iSpec, NEL_

  ! Usual debug variables.
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='read_geomlt_file'

  type(TimeType) :: TimeBoundary

  integer :: FileIndexStart, FileIndexEnd, FileIndex
  integer :: NTL_
  integer :: year, month, day, hour, minute, second, msec
  character(len=25)  :: ISOString, tempString
  character(len=3)   :: MLTString, NSCString
  character(len=1)   :: sa

  character(len=25), allocatable :: TimeBuffer(:)
  real(DP), allocatable :: MLTBuffer(:), EnergyBuffer(:)
  real(DP), allocatable :: NSCBuffer(:,:,:), FluxBuffer(:,:,:)

  !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Build file name using current date.  NameParticle decides if we open
  ! electrons or protons.
  if (NameParticle.eq.'prot') then
     NEL_ = NEL_prot
  elseif (NameParticle.eq.'elec') then
     NEL_ = NEL
  else
     call CON_stop(NameSub//': Invalid particle type '//'"'//NameParticle)
  endif

  select case (boundary)
     case('LANL')
        write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a)') trim(BoundaryPath)//'/', &
              TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
              '_mpa-sopa_', NameParticle, '_geomlt_5-min.txt'
     case('PTM')
        write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a)') trim(BoundaryPath)//'/', &
              TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
              '_ptm_', NameParticle, '_geomlt_5-min.txt'
     case('QDKP')
        write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a)') trim(BoundaryPath)//'/', &
              TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
              '_qdm-kp_', NameParticle, '_geomlt_5-min.txt'
     case('QDMVBZ')
        write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a)') trim(BoundaryPath)//'/', &
              TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
              '_qdm-vbz_', NameParticle, '_geomlt_5-min.txt'
  end select

  allocate(TimeBuffer(1))

  open(unit=UNITTMP_, file=NameFileIn, status='OLD', iostat=iError)
  if(iError/=0) call CON_stop(NameSub//': Error opening file '//NameFileIn)
  FileIndexStart = 0
  FileIndexEnd   = 0
  FileIndex      = 0
  Read_BoundaryFile_Dates: DO
     read(UNITTMP_,*,IOSTAT=iError) TimeBuffer(1)
     if ((trim(TimeBuffer(1)).ne.'#').and.(FileIndexStart.eq.0)) FileIndexStart = FileIndex
     if (iError.lt.0) then
        FileIndexEnd = FileIndex
        exit Read_BoundaryFile_Dates
     else
        FileIndex = FileIndex + 1
        cycle Read_BoundaryFile_Dates
     endif
  ENDDO Read_BoundaryFile_Dates
  FileIndex = FileIndexEnd-FileIndexStart-1
  close(UNITTMP_)

  deallocate(TimeBuffer)

  NTL_ = NTL-1
  NBD = FileIndex/NTL_
  allocate(TimeBUffer(NBD), MLTBuffer(NTL_), EnergyBuffer(NEL_), &
           NSCBuffer(NBD,NTL_,NEL_), FluxBuffer(NBD,NTL_,NEL_))

  open(unit=UNITTMP_, file=NameFileIn, status='OLD', iostat=iError)
  if (FileIndexStart.gt.0) then
     do i=1,FileIndexStart
        read(UNITTMP_,*) StringHeader
     enddo
  endif

  read(UNITTMP_,*) ISOString, MLTString, NSCString, EnergyBuffer(:)
  do i=1,NBD
     do j=1,NTL_
        read(UNITTMP_,*) TimeBuffer(i), MLTBuffer(j), NSCBuffer(i,j,:), FluxBuffer(i,j,:)
     enddo
  enddo

  if(DoTest)then
     write(*,*) NameSub//': File read successfully.'
     write(*,*) NameSub//': Energy Grid='
     write(*,*) EnergyBuffer
     write(*,*) NameSub//': LT Grid='
     write(*,*) MLTBuffer
     write(*,*) NameSub//': First data line='
     write(*,*) FluxBuffer(1,1,:)
     write(*,*) NameSub//': Last data line='
     write(*,*) FluxBuffer(NBD,24,:)
  end if

  ! Allocate global arrays
  if(.not.allocated(flux_SIII)) then
     allocate(flux_SIII(2,NBD,NTL_,NEL_), fluxLast_SII(2,NTL_,NEL_), eGrid_SI(2,NEL_), &
              avgSats_SI(2,NBD), tGrid_SI(2,NBD), lGrid_SI(2,0:NTL))
     flux_SIII = 0.
     fluxLast_SII = 0.
     eGrid_SI = 0.
  endif

  ! Set species index.
  iSpec=1 !Protons.
  if(NameParticle.eq.'elec') iSpec=2 !electrons.

  ! Put data into proper array according to NameParticle.
  lGrid_SI(iSpec,1:NTL-1)     = MLTBuffer
  eGrid_SI(iSpec,1:NEL_)      = EnergyBuffer
  flux_SIII(iSpec,:,:,1:NEL_) = FluxBuffer

  ! Wrap grid about LT=0.
  lGrid_SI(iSpec, 0) = 2.*lGrid_SI(iSpec, 1)-lGrid_SI(iSpec, 2)
  lGrid_SI(iSpec,NTL) = 2.*lGrid_SI(iSpec,NTL-1)-lGrid_SI(iSpec,NTL-2)

  ! Convert ISO time format into usable time (seconds from simulation start)
  do i=1,NBD
     tempString = TimeBuffer(i)
     read(tempString,11) year, sa, month, sa, day, sa, hour, sa, minute, sa, &
                         second, sa, msec, sa
     TimeBoundary % iYear   = year
     TimeBoundary % iMonth  = month
     TimeBoundary % iDay    = day
     TimeBoundary % iHour   = hour
     TimeBoundary % iMinute = minute
     TimeBoundary % iSecond = second
     TimeBoundary % FracSecond = 0.0
     call time_int_to_real(TimeBoundary)
     tGrid_SI(iSpec,i) = TimeBoundary % Time - TimeRamStart % Time
  enddo 

  ! Finally, store last entry into fluxLast_SII.
  fluxLast_SII(iSpec,:,:) = flux_SIII(iSpec,NBD,:,:)

  deallocate(TimeBuffer, MLTBuffer, EnergyBuffer, FluxBuffer, NSCBuffer)

11 FORMAT(I4,A,I2,A,I2,A,I2,A,I2,A,I2,A,I3,A)
  return
end subroutine read_geomlt_file

!===========================================================================
  subroutine read_hI_file
     !!!! Module Variables
     use ModRamMain,      ONLY: DP,PathScbIn
     use ModRamTiming,    ONLY: TimeRamNow
     use ModRamParams,    ONLY: NameBoundMag, UseEfInd
     use ModRamGrids,     ONLY: NR, NT, NPA
     use ModRamVariables, ONLY: FNHS, FNIS, BOUNHS, BOUNIS, BNES, HDNS, dIdt, &
                                dBdt, dIbndt, LZ, MU, PAbn, EIR, EIP
     !!!! Module Subroutines/Functions
     use ModRamFunctions, ONLY: COSD, FUNT, FUNI
     !!!! Share Modules
     use ModIoUnit,       ONLY: UNITTMP_

     implicit none

     logical :: THERE
     integer :: I, J, L
     integer :: IPA

     real(DP) :: RRL, PH

     character(len=100) :: HEADER
     character(len=213) :: hIFile

     THERE = .false.
     IF (NameBoundMag .EQ. 'DIPL') THEN
        hIfile=trim(PathScbIn)//'hI_dipole.dat'
     ELSE
        hIFile=trim(PathScbIn)//RamFileName('hI_output','dat',TimeRamNow)
        INQUIRE(FILE=trim(hIFile), EXIST=THERE)
        if (.not.THERE) hIfile=trim(PathScbIn)//'hI_dipole.dat'
     END IF

     OPEN(UNITTMP_,FILE=trim(hIfile),STATUS='OLD')
     READ(UNITTMP_,'(A)') HEADER
     READ(UNITTMP_,'(A)') HEADER
     DO I=2,NR+1
        DO J=1,NT
           DO L=1,NPA
              if (UseEfind) then
                 READ(UNITTMP_,*) RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L), &
                                  FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J), &
                                  HDNS(I,J,L),EIR(I,J),EIP(I,J)
              else
                 READ(UNITTMP_,*) RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L), &
                                  FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J), &
                                  HDNS(I,J,L)
                 EIR(I,J)=0.
                 EIP(I,J)=0.
              endif
              dIdt(I,J,L)=0.
              dIbndt(I,J,L)=0.
           ENDDO
           IF (J.LT.8.or.J.GT.18) THEN
              DO L=15,2,-1
                 if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
                    FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
                 endif
             ENDDO
           ENDIF
           BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
           dBdt(I,J)=0.
        ENDDO
     ENDDO
     CLOSE(UNITTMP_)

     DO J=1,NT ! use dipole B at I=1
        BNES(1,J)=0.32/LZ(1)**3/1.e4
        dBdt(1,J) = 0.
        EIR(1,J) = 0.
        EIP(1,J) = 0.
        DO L=1,NPA
           FNHS(1,J,L) = FUNT(MU(L))
           FNIS(1,J,L) = FUNI(MU(L))
           BOUNHS(1,J,L)=FUNT(COSD(PAbn(L)))
           BOUNIS(1,J,L)=FUNI(COSD(PAbn(L)))
           HDNS(1,J,L)=HDNS(2,J,L)
           dIdt(1,J,L)=0.
           dIbndt(1,J,L)=0.
        ENDDO
     ENDDO

     return

  end subroutine read_hI_file

!==============================================================================
  subroutine read_initial

    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: nR, nT, nE, nPA, RadiusMax, RadiusMin
    use ModRamVariables, ONLY: F2, EkeV, Lz, MLT, Pa, &
                               PParT, PPerT, PPerO, PParO, PPerE, PParE, &
                               PPerHe, PParHe, PPerH, PParH
    use ModRamParams,    ONLY: InitializationPath

    use ModRamGSL, ONLY: GSL_Interpolation_2D
    use ModRamNCDF, ONLY: ncdf_check

    use netcdf
    use nrtype, ONLY: pi_d

    implicit none

    integer :: i, j, k, l, iS, GSLerr, iRDim, iTDim, iEDim, iPaDim, &
               iR, iT, iE, iPa
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, &
               iFileID, iStatus, iPParTVar, iPPerTVar
    real(DP), allocatable :: iF2(:,:,:,:,:), iFNHS(:,:,:), iFNIS(:,:,:), iBOUNHS(:,:,:), &
                             iBOUNIS(:,:,:), iEIR(:,:), iEIP(:,:), iBNES(:,:), &
                             iHDNS(:,:,:), idBdt(:,:), idIdt(:,:,:), idIbndt(:,:,:), &
                             iLz(:), iMLT(:), iEkeV(:), iPaVar(:), radGrid(:,:), &
                             angleGrid(:,:), iPParT(:,:,:), iPPerT(:,:,:)
    real(DP) :: DL1, DPHI

    character(len=100) :: NameFile

    character(len=*), parameter :: NameSub='read_initial'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call write_prefix
    NameFile = trim(InitializationPath)//'/initialization.nc'
    write(*,*) 'Loading initial condition files '//trim(NameFile)
    ! LOAD INITIAL FILE
    iStatus = nf90_open(trim(NameFile),nf90_nowrite,iFileID)
    call ncdf_check(iStatus, NameSub)

    ! GET DIMENSIONS
    iStatus = nf90_inq_dimid(iFileID, 'nR',   iRDim)
    iStatus = nf90_inquire_dimension(iFileID, iRDim,  len = iR)
    iStatus = nf90_inq_dimid(iFileID, 'nT',   iTDim)
    iStatus = nf90_inquire_dimension(iFileID, iTDim,  len = iT)
    iStatus = nf90_inq_dimid(iFileID, 'nE',   iEDim)
    iStatus = nf90_inquire_dimension(iFileID, iEDim,  len = iE)
    iStatus = nf90_inq_dimid(iFileID, 'nPa',  iPaDim)
    iStatus = nf90_inquire_dimension(iFileID, iPaDim, len = iPa)

    ! ALLOCATE ARRAYS FOR DATA
    ALLOCATE(iF2(4,iR,iT,iE,iPa), iLz(iR+1), iMLT(iT), iEkeV(iE), iPaVar(iPa))
    ALLOCATE(iFNHS(iR+1,iT,iPa), iBOUNHS(iR+1,iT,iPa), idBdt(iR+1,iT), iHDNS(iR+1,iT,iPa), &
             iFNIS(iR+1,iT,iPa), iBOUNIS(iR+1,iT,iPa), idIdt(iR+1,iT,iPa), iBNES(iR+1,iT), &
             idIbndt(iR+1,iT,iPa), iEIR(iR+1,iT), iEIP(iR+1,iT), iPParT(4,iR,iT), &
             iPPerT(4,iR,iT))

    ! GET VARIABLE IDS
    !! FLUXES
    iStatus = nf90_inq_varid(iFileID, 'FluxE',  iFluxEVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxH',  iFluxHVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxHe', iFluxHeVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxO',  iFluxOVar)

    !! PRESSURES
    iStatus = nf90_inq_varid(iFileID, 'PParT', iPParTVar)
    iStatus = nf90_inq_varid(iFileID, 'PPerT', iPPerTVar)

    ! READ DATA
    !! FLUXES
    iStatus = nf90_get_var(iFileID, iFluxEVar,  iF2(1,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHVar,  iF2(2,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHeVar, iF2(3,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxOVar,  iF2(4,:,:,:,:))

    !! PRESSURES
    iStatus = nf90_get_var(iFileID, iPParTVar, iPParT(:,:,:))
    iStatus = nf90_get_var(iFileID, iPPerTVar, iPPerT(:,:,:))

    ! CLOSE INITIALIZATION FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

    DL1 = (RadiusMax - RadiusMin)/(iR - 1)
    DO I=1,iR+1
      iLz(I)=2.+(I-2)*DL1
    END DO

    DPHI=2.*PI_d/(iT-1)
    DO J=1,iT
      iMLT(J)=(J-1)*DPHI
    END DO

    iEkeV = EkeV
    iPaVar = Pa

    ! Now we need to check the dimensions of the initialization file and
    ! interpolate if they are different
    if ((nR.eq.iR).and.(nT.eq.iT)) then
       F2     = iF2
       PParT  = iPParT
       PPerT  = iPPerT
    else
       ! Interpolate spatially for each energy and pitch angle (if required)
       ALLOCATE(radGrid(nR+1,nT),angleGrid(nR+1,nT))
       DO i=1,NR+1
          radGrid(i,:) = Lz(i)
       ENDDO
       DO j=1,NT
          angleGrid(:,j) = MLT(j)*2*PI_d/24
       ENDDO
       do iS=1,4
          CALL GSL_Interpolation_2D(iLz(1:iR), iMLT, iPParT(iS,:,:), radGrid(1:nR,:), &
                                    angleGrid(1:nR,:), PParT(iS,:,:), GSLerr)
          CALL GSL_Interpolation_2D(iLz(1:iR), iMLT, iPPerT(iS,:,:), radGrid(1:nR,:), & 
                                    angleGrid(1:nR,:), PPerT(iS,:,:), GSLerr)
       enddo
       do l=1,nPa
          do k=1,nE
             do iS=1,4
                CALL GSL_Interpolation_2D(iLz(1:iR), iMLT, iF2(iS,:,:,k,l), radGrid(1:nR,:), &
                                          angleGrid(1:nR,:), F2(iS,1:nR,:,k,l), GSLerr)
             enddo
          enddo
       enddo
       DEALLOCATE(radGrid,angleGrid)

       ! Interpolate across pitch angle
       if ((nPa.ne.iPa).or.(nE.ne.iE)) then
        call CON_stop('Changing pitch angle and energy resolution not currently supported')
       endif
    endif

    PPerO = PPerT(4,:,:)
    PParO = PParT(4,:,:)
    PPerHe = PPerT(3,:,:)
    PParHe = PParT(3,:,:)
    PPerE = PPerT(1,:,:)
    PParE = PParT(1,:,:)
    PPerH = PPerT(2,:,:)
    PParH = PParT(2,:,:)

    DEALLOCATE(iF2, iEkeV, iMLT, iLz, iPaVar)
    DEALLOCATE(iFNHS, iBOUNHS, iFNIS, iBOUNIS, iBNES, iHDNS, iEIR, iEIP, &
               idBdt, idIbndt, idIdt, iPParT, iPPerT)

  end subroutine read_initial

!============================!
!==== OUTPUT SUBROUTINES ====!
!============================!
!==========================================================================
  subroutine ram_epot_write

    use ModRamMain, ONLY: DP, PathRamOut
    use ModRamTiming, ONLY: TimeRamNow, TimeRamElapsed
    use ModRamGrids, ONLY: nR, nT
    use ModRamVariables, ONLY: Kp, F107, VT, LZ, MLT

    use ModRamFunctions, ONLY: RamFileName

    use ModIOUnit, ONLY: UNITTMP_

    implicit none

    character(len=23)  :: StringDate
    character(len=300) :: NameFileOut
    integer :: i, j
    real(DP) :: T
   
    T = TimeRamElapsed
    write(StringDate,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2))") &
          TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
          TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond

   ! Write the electric potential [kV]
    NameFileOut=trim(PathRamOut)//RamFileName('efield','in',TimeRamNow)
    OPEN(UNIT=UNITTMP_,FILE=NameFileOut,STATUS='UNKNOWN')
    WRITE(UNITTMP_,*)'UT(hr)    Kp   F107   VTmax   VTmin   Date=',StringDate
    WRITE(UNITTMP_,*) T/3600.,KP,F107,maxval(VT)/1e3,minval(VT)/1e3
    WRITE(UNITTMP_,*) ' L     MLT        Epot(kV)'
    DO I=1,NR+1
      DO J=1,NT
        WRITE (UNITTMP_,*) LZ(I),MLT(J),VT(I,J)/1e3
      END DO
    END DO
    CLOSE(UNITTMP_)

    RETURN

  end subroutine ram_epot_write

!==========================================================================
  subroutine ram_hour_write !Previously WRESULT in ram_all

    use ModRamMain,      ONLY: DP, S, PathRamOut
    use ModRamGrids,     ONLY: nR, nE, nPA, nT
    use ModRamVariables, ONLY: F2, FFACTOR, FNHS, WE, WMU, XNN, XND, ENERN,   &
                               ENERD, EkeV, UPA, LNCN, LNCD, Kp, MLT, LZ, &
                               LECD, LECN
    use ModRamTiming,    ONLY: TimeRamNow

    use ModRamFunctions

    use ModIOUnit, ONLY: UNITTMP_, io_unit_new

    implicit none

    integer :: i, j, k, l, jw, iw
    real(DP) :: weight
    real(DP), ALLOCATABLE :: XNNO(:),XNDO(:)
    real(DP), ALLOCATABLE :: F(:,:,:,:),FZERO(:,:,:),ENO(:),EDO(:),AVEFL(:,:,:),BARFL(:)
    character(len=23) :: StringDate
    character(len=2)  :: ST2
    character(len=214) :: ST4, NameFileOut
    character(len=2), dimension(4) :: speciesString = (/'_e','_h','he','_o'/)

    ALLOCATE(F(NR,NT,NE,NPA),FZERO(NR,NT,NE),ENO(NR),EDO(NR),AVEFL(NR,NT,NE),BARFL(NE))
    ALLOCATE(XNNO(NR),XNDO(NR))
    ST2 = speciesString(S)
    ST4 =trim(PathRamOut)
    write(StringDate,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2))") &
          TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
          TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond

    ! Print RAM output at every hour

    DO I=2,NR
      XNNO(I)=XNN(S,I)
      XNDO(I)=XND(S,I)
      XNN(S,I)=0.
      XND(S,I)=0.
      ENO(I)=ENERN(S,I)
      EDO(I)=ENERD(S,I)
      ENERN(S,I)=0.
      ENERD(S,I)=0.
      DO K=2,NE
        DO L=1,NPA
          DO J=1,NT-1
            IF (L.EQ.1) THEN
              F(I,J,K,1)=F2(S,I,J,K,2)/FFACTOR(S,I,K,2)/FNHS(I,J,2)
            ELSE
              F(I,J,K,L)=F2(S,I,J,K,L)/FFACTOR(S,I,K,L)/FNHS(I,J,L)
            ENDIF
            if (f2(S,i,j,k,l).lt.1E-5) f2(S,i,j,k,l)=1E-5
            if (f(i,j,k,l).lt.1E-5) f(i,j,k,l)=1E-5
            IF (L.LT.UPA(I)) THEN
              WEIGHT=F2(S,I,J,K,L)*WE(K)*WMU(L)
              IF (MLT(J).LE.6.OR.MLT(J).GE.18.) THEN
                XNN(S,I)=XNN(S,I)+WEIGHT
                ENERN(S,I)=ENERN(S,I)+EKEV(K)*WEIGHT
              ELSE
                XND(S,I)=XND(S,I)+WEIGHT
                ENERD(S,I)=ENERD(S,I)+EKEV(K)*WEIGHT
              ENDIF
            ENDIF
          END DO
          F(I,NT,K,L)=F(I,1,K,L)
        END DO
      END DO
      LNCN(S,I)=XNNO(I)-XNN(S,I)
      LECN(S,I)=ENO(I)-ENERN(S,I)
      LNCD(S,I)=XNDO(I)-XND(S,I)
      LECD(S,I)=EDO(I)-ENERD(S,I)
    END DO

    ! Write the trapped equatorial flux [1/s/cm2/sr/keV]
    if (NT.EQ.49) JW=6
    if (NT.EQ.25) JW=3
    IW=4
    NameFileOut=trim(PathRamOut)//RamFileName('ram'//st2,'t',TimeRamNow)
    open(unit=UnitTMP_, file=trim(NameFileOut), status='UNKNOWN')
    DO I=4,NR,IW
      DO 25 J=1,NT-1,JW
        WRITE(UNITTMP_,32) StringDate,LZ(I),KP,MLT(J)
        DO 27 K=4,NE-1
27      WRITE(UNITTMP_,30) EKEV(K),(F(I,J,K,L),L=2,NPA-2)
25    CONTINUE
    END DO
    close(UNITTMP_)

! WPI is not currently working, so these write statements are commented out -ME
    ! Write the total precipitating flux [1/cm2/s]
!    IF (DoUseWPI) THEN
!      OPEN(UNIT=UNITTMP_,FILE=trim(ST4)//ST0//ST1//ST2//'.tip',STATUS='UNKNOWN')
!      WRITE(UNITTMP_,71) T/3600,KP
!      DO I=2,NR
!        DO J=1,NT
!          PRECFL=0.
!          DO K=2,NE ! 0.15 - 430 keV
!            AVEFL(I,J,K)=0.
!            DO L=UPA(I),NPA
!              AVEFL(I,J,K)=AVEFL(I,J,K)+F(I,J,K,L)*WMU(L)
!            ENDDO
!            AVEFL(I,J,K)=AVEFL(I,J,K)/(MU(NPA)-MU(UPA(I)))
!            PRECFL=PRECFL+AVEFL(I,J,K)*PI*WE(K)
!          ENDDO
!          WRITE(UNITTMP_,70) LZ(I),PHI(J),PRECFL
!        END DO
!      END DO
!      close(UNITTMP_)
!    ENDIF

    ! Write the RAM flux [1/s/cm2/sr/keV] at given pitch angle
!    IF (DoUseWPI) THEN
!      OPEN(UNIT=UNITTMP_,FILE=trim(PathRamOut)//'outm'//ST0//ST2//ST1//'.dat',STATUS='UNKNOWN')
!      DO I=1,NR
!        DO 822 J=1,NT-1
!          DO 822 K=2,NE
!            WRITE(UNITTMP_,31) T/3600.,LZ(I),KP,MLT(J),EKEV(K),F(I,J,K,27)
!822         CONTINUE
!      ENDDO
!      CLOSE(UNITTMP_)
!    ENDIF

30  FORMAT(F7.2,72(1PE11.3))
32  FORMAT(' EKEV/PA, Date=',a,' L=',F6.2,' Kp=',F6.2,' MLT=',F4.1)

    DEALLOCATE(F,FZERO,ENO,EDO,AVEFL,BARFL,XNNO,XNDO)

    RETURN
  END SUBROUTINE ram_hour_write


!===========================================================================
  subroutine ram_write_pressure
! Creates RAM pressure files (output_ram/pressure_d{Date and Time}.dat)
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRamOut
    use ModRamTiming,    ONLY: TimeRamNow, TimeRamStart, TimeRamElapsed
    use ModRamParams,    ONLY: DoAnisoPressureGMCoupling
    use ModRamGrids,     ONLY: NR, NT
    use ModRamVariables, ONLY: PParH, PPerH, PParHe, PPerHe, PParO, PPerO, &
                               PPerE, PParE, PAllSum, PParSum, KP, LZ, PHI
    !!!! Share Modules
    use ModIOUnit, ONLY: UNITTMP_

    use nrtype, ONLY: pi_d
    implicit none

    character(len=23)            :: StringTime
    character(len=*), parameter  :: NameSub = 'ram_write_pressure'
    character(len=200)           :: FileName
    integer                      :: iError, i, j
    !------------------------------------------------------------------------
    ! Create Ram pressure output file.
    FileName = trim(PathRamOut)//trim(RamFileName('/pressure','dat',TimeRamNow))
    open(unit=UNITTMP_, status='REPLACE', iostat=iError, file=FileName)
    if(iError /= 0) call CON_stop(NameSub//' Error opening file '//FileName)

    ! Prepare date/time for header.
    write(StringTime,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2))") &
             TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
             TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond

    ! Write pressure to file.
    write(UNITTMP_,'(a, a, a3,f8.3,2x,a4,f3.1)') 'Date=', StringTime, ' T=', &
             TimeRamElapsed/3600.0 + TimeRamStart%iHour, ' Kp=', Kp
    if (.not.DoAnisoPressureGMCoupling) then
       write(UNITTMP_,'(2a)') ' Lsh MLT PPER_H PPAR_H PPER_O PPAR_O PPER_He PPAR_He', &
                              ' PPER_E  PPAR_E   PTotal   [keV/cm3]'
    else
       write(UNITTMP_,'(2a)') ' Lsh   MLT    PPER_H   PPAR_H  ', &
                              ' PPER_O   PPAR_O   PPER_He   PPAR_He ', &
                              ' PPER_E   PPAR_E   PTotal   Ppar [keV/cm3]'
    end if
    do i=2, NR
       do j=1, NT
          if (.not.DoAnisoPressureGMCoupling) then
             write(UNITTMP_,*) LZ(I),PHI(J)*12/PI_d,PPERH(I,J),PPARH(I,J), &
                               PPERO(I,J),PPARO(I,J), PPERHE(I,J),PPARHE(I,J), &
                               PPERE(I,J),PPARE(I,J),PAllSum(i,j)
          else
             write(UNITTMP_,*) LZ(I),PHI(J)*12/PI_d,PPERH(I,J),PPARH(I,J), &
                               PPERO(I,J),PPARO(I,J),PPERHE(I,J),PPARHE(I,J), &
                               PPERE(I,J),PPARE(I,J),PAllSum(i,j),PparSum(i,j)
          end if
       enddo
    enddo
    close(UNITTMP_)

  end subroutine ram_write_pressure

!==============================================================================
subroutine ram_write_hI

  use ModRamMain,      ONLY: PathScbOut
  use ModRamGrids,     ONLY: nR, nT, nPA
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed
  use ModRamVariables, ONLY: LZ, MLT, FNHS, BOUNHS, FNIS, BOUNIS, BNES, HDNS

  use ModIOUnit, ONLY: UNITTMP_

  implicit none

  integer :: i, j, L
  character(len=214) :: filenamehI

  if (abs(TimeRamElapsed).le.1e-9) then
     filenamehI = trim(PathScbOut)//'/hI_output_dipole.dat'
  else
     filenamehI=trim(PathScbOut)//RamFileName('/hI_output', 'dat', TimeRamNow)
  endif
  PRINT*, 'Writing to file ', TRIM(filenamehI)
  OPEN(UNITTMP_, file = TRIM(filenamehI), action = 'write', status = 'unknown')
  WRITE(UNITTMP_, *) '   Lsh        MLT    NPA    h(Lsh,MLT,NPA) hBoun(Lsh,MLT,NPA) '//&
             'I(Lsh,MLT,NPA) IBoun(Lsh,MLT,NPA) Bz(Lsh,MLT) HDENS(Lsh,MLT,NPA)'
  DO i = 2, nR+1
     DO j = 1, nT
        DO L = 1, NPA
           WRITE(UNITTMP_, *) LZ(i), MLT(j), L, FNHS(i,j,L), BOUNHS(i,j,L), &
                FNIS(i,j,L), BOUNIS(i,j,L), BNES(i,j), HDNS(i,j,L)
        END DO
     END DO
  END DO
  CLOSE(UNITTMP_)

end subroutine ram_write_hI

!==============================================================================
subroutine write_dsbnd
! Creates boundary flux files (output_ram/Dsbnd/ds_{Species}_d{Date and Time}.dat)
  !!!! Module Variables
  use ModRamMain,      ONLY: PathRamOut, S
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed
  use ModRamGrids,     ONLY: NE, NR, NT
  use ModRamVariables, ONLY: EKEV, FFACTOR, FGEOS, KP, F107
  !!!! Share Modules
  use ModIoUnit,      ONLY: UNITTMP_

  implicit none

  integer :: K, j
  character(len=2), DIMENSION(4) :: ST2 = (/ '_e','_h','he','_o' /)

  character(len=214) :: NameFluxFile

  NameFluxFile=trim(PathRamOut)//RamFileName('Dsbnd/ds'//St2(S),'dat',TimeRamNow)
  OPEN(UNIT=UNITTMP_,FILE=NameFluxFile, STATUS='UNKNOWN')
  WRITE(UNITTMP_,*)'EKEV FGEOSB [1/cm2/s/ster/keV] T=',TimeRamElapsed/3600,Kp,F107
  DO K=2,NE
     WRITE(UNITTMP_,*) EKEV(K),(FGEOS(S,J,K,2)/FFACTOR(S,NR,K,2),J=1,NT)
  END DO
  CLOSE(UNITTMP_)

end Subroutine write_dsbnd

!==============================================================================
    
END MODULE ModRamIO
