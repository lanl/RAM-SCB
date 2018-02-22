!============================================================================
module ModRamIO
!   A set of variables and tools for general I/O handling, file name
!   creation, etc.
!   Copyright (c) 2016, Los Alamos National Security, LLC
!   All rights reserved.
!******************************************************************************

  use ModRamMain,   ONLY: Real8_, S, nIter
  use ModRamConst,  ONLY: PI
  use ModRamParams, ONLY: DoSaveRamSats, UseNewFmt
  use ModRamNCDF,   ONLY: ncdf_check, write_ncdf_globatts
  use ModRamTiming, ONLY: DtLogFile, DtWriteSat

  implicit none
  save

  logical :: IsFramework = .false. 

  ! File output names and units
  integer :: iUnitLog ! Logfile IO unit
  character(len=100) :: NameFileLog, NameFileDst

  ! String that contains the date and time for which this simulation
  ! was initialized (used for output metadata.)
  character(len=21) :: StringRunDate
  
  character(len=10) :: NameMod = 'ModRamMain'
  
contains
!============================!
!===== BASE SUBROUTINES =====!
!============================!
!==============================================================================
  function RamFileName(PrefixIn, SuffixIn, TimeIn)
    ! Create a file name using the new RAM-SCB I/O filename standards:
    ! FileNameOut = PrefixIn_dYYYYMMDD_tHHMMSS.SuffixIn
    ! PrefixIn and SuffixIn are short, descriptive strings (e.g. 
    ! 'pressure', 'ram_o' for prefixes and 'dat', 'cdf' for suffixes.

    use ModTimeConvert

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

    use ModRamSats,   ONLY: read_sat_input, init_sats
    use ModRamMain,   ONLY: PathRamOut, PathScbOut

    use ModIoUnit,    ONLY: UNITTMP_

    character(len=*), parameter :: NameSub = 'init_output'
    !------------------------------------------------------------------------
    ! Initialize logfile.
    write(NameFileLog, '(a,i6.6, a)') PathRamOut//'log_n', nIter, '.log'
    open(UNITTMP_, FILE=NameFileLog, STATUS='REPLACE')
    write(UNITTMP_,*) 'RAM-SCB Log'
    write(UNITTMP_,*) 'time year mo dy hr mn sc msc dstRam dstBiot ', &
                      'pparh pperh pparo ppero pparhe pperhe ppare ppere'
    close(UNITTMP_)

    ! Initialize Dstfile.
    write(NameFileDst, '(a,i6.6, a)') PathScbOut//'dst_computed_3deq_', nIter, '.txt'
    open(UNITTMP_, FILE=NameFileDst, STATUS='REPLACE')
    write(UNITTMP_,*) 'SCB Dst File'
    write(UNITTMP_,*) 'time year mo dy hr mn sc msc dstDPS dstDPSGeo ', &
                      'dstBiot dstBiotGeo'
    close(UNITTMP_)

    if (DoSaveRamSats) then
       call read_sat_input
       call init_sats
    end if

  end subroutine init_output

!==============================================================================
  subroutine handle_output(TimeIn)
    ! For every type of output, check the following:
    ! 1) if that output is activated,
    ! 2) if it's time to write that output
    ! If so, call the proper routines to write the output.

    !!!! Module Variables
    use ModRamParams,    ONLY: DoSaveRamSats
    use ModRamTiming,    ONLY: Dt_hI, DtRestart, DtWriteSat, TimeRamNow, &
                               TimeRamElapsed, DtW_Pressure, DtW_hI, DtW_Efield
    use ModRamGrids,     ONLY: NR, NT
    use ModRamVariables, ONLY: PParH, PPerH, PParHe, PPerHe, PParE, PPerE, &
                               PParO, PPerO, VT
    use ModScbVariables, ONLY: DstBiot, DstBiotInsideGeo, DstDPS, DstDPSInsideGeo
    !!!! Module Subroutines/Functions
    use ModRamSats,      ONLY: fly_sats
    use ModRamFunctions, ONLY: ram_sum_pressure, get_ramdst
    use ModRamRestart,   ONLY: write_restart
    ! Share Modules
    use ModIOUnit, ONLY: UNITTMP_

    implicit none

    real(kind=Real8_), intent(in) :: TimeIn

    real(kind=Real8_) :: dst
    logical :: DoTest, DoTestMe
    integer :: iS
    character(len=4) :: StrScbIter
    character(len=*), parameter :: NameSub = 'handle_output'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Write Logfile
    if(mod(TimeIn, DtLogfile)==0.0)then
       ! Get current Dst.
       call get_ramdst(dst)
       open(UNITTMP_, FILE=NameFileLog, POSITION='APPEND')
       write(UNITTMP_, *) TimeIn, TimeRamNow%iYear, TimeRamNow%iMonth, &
                         TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
                         TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0), &
                         dst, DstBiot, sum(PParH)/(nR*nT), sum(PPerH)/(nR*nT), &
                         sum(PParO)/(nR*nT), sum(PPerO)/(nR*nT), sum(PParHe)/(nR*nT), &
                         sum(PPerHe)/(nR*nT), sum(PParE)/(nR*nT), sum(PPerE)/(nR*nT)
       close(UNITTMP_)
    end if

    ! Write SCB Dst file
    if(mod(TimeIn, Dt_hI)==0.0)then
       open(UNITTMP_, FILE=NameFileDst, POSITION='APPEND')
       write(UNITTMP_, *) TimeIn, TimeRamNow%iYear, TimeRamNow%iMonth, &
                          TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
                          TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0), &
                          dstDPS, dstDPSInsideGeo, DstBiot, DstBiotInsideGeo, &
                          maxval(VT), minval(VT)
       close(UNITTMP_)
    end if

    ! Write Pressure File
    if (mod(TimeIn, DtW_Pressure)==0.0) then
       write(StrScbIter,'(I4.4)') int(TimeRamElapsed/Dt_hI)
       call ram_sum_pressure
       call ram_write_pressure(StrScbIter)
    end if

    ! Write hI File
    if (mod(TimeIn, DtW_hI)==0.0) call ram_write_hI

    ! Write efield file
    if (mod(TimeIn, DtW_Efield).eq.0) call ram_epot_write

    ! Update Satellite Files
    if (DoSaveRamSats .and. (mod(TimeRamElapsed,DtWriteSat) .eq. 0)) call fly_sats

    ! Write restarts.
    if (mod(TimeIn,DtRestart).eq.0) call write_restart()

    ! Write hourly file
    if (mod(TimeIn,3600.0).eq.0) then
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
  use ModRamMain,      ONLY: PathRamIN
  use ModRamTiming,    ONLY: TimeRamNow, Dt_bc, TimeRamStart
  use ModRamGrids,     ONLY: NE, NEL, NEL_prot, NTL, NBD
  use ModRamParams,    ONLY: boundary, BoundaryPath
  use ModRamVariables, ONLY: EKEV, FGEOS, IsInitialized, timeOffset, StringFileDate,  &
                             flux_SIII, fluxLast_SII, eGrid_SI, avgSats_SI, lGrid_SI, &
                             tGrid_SI
  !!!! Share Modules
  use ModIOUnit,      ONLY: UNITTMP_
  use ModTimeConvert, ONLY: TimeType, time_int_to_real

  implicit none

  character(len=4), intent(in) :: NameParticle

  character(len=200) :: NameFileIn, StringHeader
  integer :: iError, i, j, k, iSpec, NEL_

  ! Buffers to hold read data before placing in correct location.
  real(kind=Real8_), allocatable :: Buffer_I(:), Buffer_III(:,:,:), Buffer2_I(:)
  integer, allocatable           :: iBuffer_II(:,:)

  ! Usual debug variables.
  logical :: DoTest, DoTestMe, IsOpened
  character(len=*), parameter :: NameSub='read_geomlt_file'

  character(len=99) :: NEL_char

  type(TimeType) :: TimeBoundary

  integer :: FileIndexStart, FileIndexEnd, FileIndex
  integer :: NTL_
  integer :: year, month, day, hour, minute, second, msec
  character(len=25)  :: ISOString, tempString
  character(len=3)   :: MLTString, NSCString
  character(len=1)   :: sa

  character(len=25), allocatable :: TimeBuffer(:)
  real(kind=Real8_), allocatable :: MLTBuffer(:), EnergyBuffer(:)
  real(kind=Real8_), allocatable :: NSCBuffer(:,:,:), FluxBuffer(:,:,:)

  !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Save file date.
  write(StringFileDate,'(i4.4,i2.2,i2.2)') &
       TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay

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
     case('QDMKP')
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

  deallocate(TimeBuffer, MLTBuffer, EnergyBuffer, FluxBuffer)

11 FORMAT(I4,A,I2,A,I2,A,I2,A,I2,A,I2,A,I3,A)
  return
end subroutine read_geomlt_file

!===========================================================================
  subroutine read_hI_file
     !!!! Module Variables
     use ModRamMain,      ONLY: PathScbIn
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
     save

     logical :: THERE=.false.
     integer :: I, J, K, L
     integer :: IPA, nFive

     real(kind=Real8_) :: RRL, PH

     character(len=100) :: HEADER
     character(len=200) :: hIFile

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

    use ModRamMain,      ONLY: Real8_, PathRamIn
    use ModRamGrids,     ONLY: nR, nT, nE, nPA, RadiusMax, RadiusMin
    use ModRamVariables, ONLY: F2, FNHS, FNIS, BOUNHS, BOUNIS, BNES, HDNS, EIR, &
                               EIP, dBdt, dIdt, dIbndt, EkeV, Lz, MLT, Pa
    use ModRamParams,    ONLY: InitializationPath

    use ModScbSpline,    ONLY: Spline_2D_linear

    use netcdf

    implicit none

    integer :: i, j, k, l, iS, iDomain
    integer :: iRDim, iTDim, iEDim, iPaDim, iR, iT, iE, iPa, iError
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus
   
    real(kind=Real8_), allocatable :: iF2(:,:,:,:,:), iFNHS(:,:,:), iFNIS(:,:,:), &
         iBOUNHS(:,:,:), iBOUNIS(:,:,:), iEIR(:,:), iEIP(:,:), iBNES(:,:), &
         iHDNS(:,:,:), idBdt(:,:), idIdt(:,:,:), idIbndt(:,:,:), iLz(:), iMLT(:), &
         iEkeV(:), iPaVar(:), radGrid(:,:), angleGrid(:,:), tF2(:,:,:,:,:)
    real(kind=Real8_), allocatable :: pF2(:,:,:,:,:), pFNHS(:,:,:), pFNIS(:,:,:), &
         pBOUNHS(:,:,:), pBOUNIS(:,:,:), pEIR(:,:), pEIP(:,:), pBNES(:,:), &
         pHDNS(:,:,:), pdBdt(:,:), pdIdt(:,:,:), pdIbndt(:,:,:), pLz(:), pMLT(:)

    real(kind=Real8_) :: DL1, DPHI

    character(len=100)             :: NameFile, StringLine

    character(len=*), parameter :: NameSub='read_initial'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call write_prefix
    write(*,*) 'Loading initial condition files.'
    ! LOAD INITIAL FILE
    iStatus = nf90_open(trim(InitializationPath)//'/initialization.nc',nf90_nowrite,iFileID)
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
             idIbndt(iR+1,iT,iPa), iEIR(iR+1,iT), iEIP(iR+1,iT))

    ! GET VARIABLE IDS
    !! FLUXES
    iStatus = nf90_inq_varid(iFileID, 'FluxE',  iFluxEVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxH',  iFluxHVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxHe', iFluxHeVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxO',  iFluxOVar)

    !! hI OUTPUTS
    iStatus = nf90_inq_varid(iFileID, 'FNHS',   iHVar)
    iStatus = nf90_inq_varid(iFileID, 'BOUNHS', iBHVar)
    iStatus = nf90_inq_varid(iFileID, 'FNIS',   iIvar)
    iStatus = nf90_inq_varid(iFileID, 'BOUNIS', iBIVar)
    iStatus = nf90_inq_varid(iFileID, 'BNES',   iBNESVar)
    iStatus = nf90_inq_varid(iFileiD, 'HDNS',   iHDNSVar)
    iStatus = nf90_inq_varid(iFileID, 'EIR',    iEIRVar)
    iStatus = nf90_inq_varid(iFileID, 'EIP',    iEIPVar)
    iStatus = nf90_inq_varid(iFileID, 'dBdt',   iDBDTVar)
    iStatus = nf90_inq_varid(iFileID, 'dIdt',   iDIDTVar)
    iStatus = nf90_inq_varid(iFileID, 'dIbndt', iDIBNVar)

    ! READ DATA
    !! FLUXES
    iStatus = nf90_get_var(iFileID, iFluxEVar,  iF2(1,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHVar,  iF2(2,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHeVar, iF2(3,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxOVar,  iF2(4,:,:,:,:))

    !! hI OUTPUTS
    iStatus = nf90_get_var(iFileID, iHVar,    iFNHS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBHVar,   iBOUNHS(:,:,:))
    iStatus = nf90_get_var(iFileID, iIVar,    iFNIS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBIVar,   iBOUNIS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBNESVar, iBNES(:,:))
    iStatus = nf90_get_var(iFileID, iHDNSVar, iHDNS(:,:,:))
    iStatus = nf90_get_var(iFileID, iEIRVar,  iEIR(:,:))
    iStatus = nf90_get_Var(iFileID, iEIPVar,  iEIP(:,:))
    iStatus = nf90_get_var(iFileID, iDBDTVar, idBdt(:,:))
    iStatus = nf90_get_var(iFileID, iDIDTVar, idIdt(:,:,:))
    iStatus = nf90_get_var(iFileID, iDIBNVar, idIbndt(:,:,:))

    ! CLOSE INITIALIZATION FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

    DL1 = (RadiusMax - RadiusMin)/(iR - 1)
    DO I=1,iR+1
      iLz(I)=2.+(I-2)*DL1
    END DO

    DPHI=2.*PI/(iT-1)
    DO J=1,iT
      iMLT(J)=(J-1)*DPHI
    END DO

    iEkeV = EkeV
    iPaVar = Pa

    ! Now we need to check the dimensions of the initialization file and
    ! interpolate if they are different
    if ((nR.eq.iR).and.(nT.eq.iT).and.(nE.eq.iE).and.(nPa.eq.iPa)) then
       F2 = iF2(:,1:iR,:,:,:)
       FNHS = iFNHS
       FNIS = iFNIS
       BOUNHS = iBOUNHS
       BOUNIS = iBOUNIS
       BNES = iBNES
       HDNS = iHDNS
       EIR = iEIR
       EIP = iEIP
       dBdt = idBdt
       dIdt = idIdt
       dIbndt = idIbndt
    else
       ! Interpolate spatially for each energy and pitch angle (if required)
       if ((nR.ne.iR).or.(nT.ne.iT)) then
          ALLOCATE(radGrid(nR+1,nT),angleGrid(nR+1,nT))
          ALLOCATE(pF2(4,iR+1,iT+2,iE,iPa), pLz(iR+2), pMLT(iT+2))
          ALLOCATE(pFNHS(iR+2,iT+2,iPa), pBOUNHS(iR+2,iT+2,iPa), pdBdt(iR+2,iT+2), pHDNS(iR+2,iT+2,iPa), &
                   pFNIS(iR+2,iT+2,iPa), pBOUNIS(iR+2,iT+2,iPa), pdIdt(iR+2,iT+2,iPa), pBNES(iR+2,iT+2), &
                   pdIbndt(iR+2,iT+2,iPa), pEIR(iR+2,iT+2), pEIP(iR+2,iT+2))
          ALLOCATE(tF2(4,nR+1,nT,nE,nPa))
          DO i=1,NR+1
             radGrid(i,:) = Lz(i)
          ENDDO
          DO j=2,NT
             angleGrid(:,j) = MLT(j)*2*PI/24
          ENDDO
          DL1 = (RadiusMax - RadiusMin)/(iR - 1)
          DO I=1,iR+2
             pLz(I)=2.+(I-2)*DL1
          END DO
          DPHI=2.*PI/(iT-1)
          DO J=1,iT+2
             pMLT(J)=(J-2)*DPHI
          END DO
        !!!
          pBNES(1:iR+1,2:iT+1) = iBNES
          pBNES(iR+2,:) = pBNES(iR+1,:)
          pBNES(:,1) = pBNES(:,iT)
          pBNES(:,iT+2) = pBNES(:,3)
        !!!
          pdBdt(1:iR+1,2:iT+1) = idBdt
          pdBdt(iR+2,:) = pdBdt(iR+1,:)
          pdBdt(:,1) = pdBdt(:,iT)
          pdBdt(:,iT+2) = pdBdt(:,3)
        !!!
          pEIP(1:iR+1,2:iT+1) = iEIP
          pEIP(iR+2,:) = pEIP(iR+1,:)
          pEIP(:,1) = pEIP(:,iT)
          pEIP(:,iT+2) = pEIP(:,3)
        !!!
          pEIR(1:iR+1,2:iT+1) = iEIR
          pEIR(iR+2,:) = pEIR(iR+1,:)
          pEIR(:,1) = pEIR(:,iT)
          pEIR(:,iT+2) = pEIR(:,3)
        !!!
          pFNHS(1:iR+1,2:iT+1,:) = iFNHS
          pFNHS(iR+2,:,:) = pFNHS(iR+1,:,:)
          pFNHS(:,1,:) = pFNHS(:,iT,:)
          pFNHS(:,iT+2,:) = pFNHS(:,3,:)
        !!!
          pFNIS(1:iR+1,2:iT+1,:) = iFNIS
          pFNIS(iR+2,:,:) = pFNIS(iR+1,:,:)
          pFNIS(:,1,:) = pFNIS(:,iT,:)
          pFNIS(:,iT+2,:) = pFNIS(:,3,:)
        !!!
          pHDNS(1:iR+1,2:iT+1,:) = iHDNS
          pHDNS(iR+2,:,:) = pHDNS(iR+1,:,:)
          pHDNS(:,1,:) = pHDNS(:,iT,:)
          pHDNS(:,iT+2,:) = pHDNS(:,3,:)
        !!!
          pdIdt(1:iR+1,2:iT+1,:) = idIdt
          pdIdt(iR+2,:,:) = pdIdt(iR+1,:,:)
          pdIdt(:,1,:) = pdIdt(:,iT,:)
          pdIdt(:,iT+2,:) = pdIdt(:,3,:)
        !!!
          pBOUNHS(1:iR+1,2:iT+1,:) = iBOUNHS
          pBOUNHS(iR+2,:,:) = pBOUNHS(iR+1,:,:)
          pBOUNHS(:,1,:) = pBOUNHS(:,iT,:)
          pBOUNHS(:,iT+2,:) = pBOUNHS(:,3,:)
        !!!
          pBOUNIS(1:iR+1,2:iT+1,:) = iBOUNIS
          pBOUNIS(iR+2,:,:) = pBOUNIS(iR+1,:,:)
          pBOUNIS(:,1,:) = pBOUNIS(:,iT,:)
          pBOUNIS(:,iT+2,:) = pBOUNIS(:,3,:)
        !!!
          pdIbndt(1:iR+1,2:iT+1,:) = idIbndt
          pdIbndt(iR+2,:,:) = pdIbndt(iR+1,:,:)
          pdIbndt(:,1,:) = pdIbndt(:,iT,:)
          pdIbndt(:,iT+2,:) = pdIbndt(:,3,:)
        !!!
          CALL Spline_2D_linear(pLz, pMLT, pBNES, radGrid, angleGrid, BNES, iDomain)
          CALL Spline_2D_linear(pLz, pMLT, pdBdt, radGrid, angleGrid, dBdt, iDomain)
          CALL Spline_2D_linear(pLz, pMLT, pEIP, radGrid, angleGrid, EIP, iDomain)
          CALL Spline_2D_linear(pLz, pMLT, pEIR, radGrid, angleGrid, EIR, iDomain)
          do l=1,nPa
             CALL Spline_2D_linear(pLz, pMLT, pFNHS(:,:,l), radGrid, angleGrid, FNHS(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pFNIS(:,:,l), radGrid, angleGrid, FNIS(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pHDNS(:,:,l), radGrid, angleGrid, HDNS(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pdIdt(:,:,l), radGrid, angleGrid, dIdt(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pBOUNHS(:,:,l), radGrid, angleGrid, BOUNHS(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pBOUNIS(:,:,l), radGrid, angleGrid, BOUNIS(:,:,l), iDomain)
             CALL Spline_2D_linear(pLz, pMLT, pdIbndt(:,:,l), radGrid, angleGrid, dIbndt(:,:,l), iDomain)
             do k=1,nE
                do iS=1,4
                   pF2(iS,1:iR,2:iT+1,k,l) = iF2(iS,:,:,k,l)
                   pF2(iS,iR+1,:,k,l) = pF2(iS,iR,:,k,l)
                   pF2(iS,:,1,k,l) = pF2(iS,:,iT,k,l)
                   pF2(iS,:,iT+2,k,l) = pF2(iS,:,3,k,l)
                   CALL Spline_2D_linear(pLz(1:iR+1), pMLT, pF2(iS,:,:,k,l), &
                                        radGrid, angleGrid, tF2(iS,:,:,k,l), iDomain)
                   F2(iS,:,:,k,l) = tF2(iS,1:nR,:,k,l)
                enddo
             enddo
          enddo
          DEALLOCATE(radGrid,angleGrid)
          DEALLOCATE(pF2, tF2, pMLT, pLz)
          DEALLOCATE(pFNHS, pBOUNHS, pFNIS, pBOUNIS, pBNES, pHDNS, pEIR, pEIP, &
                     pdBdt, pdIbndt, pdIdt)
          HDNS = HDNS
       endif

       ! Interpolate across pitch angle
       if ((nPa.ne.iPa).or.(nE.ne.iE)) then
        call CON_stop('Changing pitch angle and energy resolution in RAM is not &
                       currently supported')
       endif
    endif

    DEALLOCATE(iF2, iEkeV, iMLT, iLz, iPaVar)
    DEALLOCATE(iFNHS, iBOUNHS, iFNIS, iBOUNIS, iBNES, iHDNS, iEIR, iEIP, &
               idBdt, idIbndt, idIdt)

  end subroutine read_initial

!============================!
!==== OUTPUT SUBROUTINES ====!
!============================!
!==========================================================================
  subroutine ram_epot_write

    use ModRamMain, ONLY: Real8_, PathRamOut
    use ModRamTiming, ONLY: TimeRamNow, TimeRamElapsed
    use ModRamGrids, ONLY: nR, nT
    use ModRamVariables, ONLY: Kp, F107, VT, LZ, MLT

    use ModRamFunctions, ONLY: RamFileName

    use ModIOUnit, ONLY: UNITTMP_

    implicit none

    character(len=23)  :: StringDate
    character(len=300) :: NameFileOut
    integer :: i, j
    real(kind=Real8_) :: T
   
    T = TimeRamElapsed
    write(StringDate,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2)'.',i3.3)"), TimeRamNow%iYear, &
          TimeRamNow%iMonth, TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
          TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0)

   ! Write the electric potential [kV]
    NameFileOut=trim(PathRamOut)//RamFileName('efield','in',TimeRamNow)
    OPEN(UNIT=UNITTMP_,FILE=NameFileOut,STATUS='UNKNOWN')
    WRITE(UNITTMP_,*)'UT(hr)    Kp   F107   VTmax   VTmin   Date=',StringDate
    WRITE(UNITTMP_,31) T/3600.,KP,F107,maxval(VT)/1e3,minval(VT)/1e3
    WRITE(UNITTMP_,*) ' L     MLT        Epot(kV)'
    DO I=1,NR+1
      DO J=1,NT
        WRITE (UNITTMP_,22) LZ(I),MLT(J),VT(I,J)/1e3
      END DO
    END DO
    CLOSE(UNITTMP_)

22  FORMAT(F5.2,F10.6,E13.4)
31  FORMAT(F6.2,1X,F6.2,2X,F6.2,1X,F7.2,1X,F7.2,1X,1PE11.3)

    RETURN

  end subroutine ram_epot_write

!==========================================================================
  subroutine ram_hour_write !Previously WRESULT in ram_all

    use ModRamMain,      ONLY: Real8_, S, PathRamOut
    use ModRamConst,     ONLY: pi
    use ModRamGrids,     ONLY: nR, nE, nPA, nT
    use ModRamVariables, ONLY: F2, FFACTOR, FNHS, WE, WMU, XNN, XND, ENERN,   &
                               ENERD, EkeV, UPA, LNCN, LNCD, Kp, MLT, MU, LZ, &
                               VT, F107, LECD, LECN
    use ModRamTiming,    ONLY: TimeRamNow

    use ModRamFunctions

    use ModIOUnit, ONLY: UNITTMP_, io_unit_new
    implicit none
    save
    integer :: i, j, k, l, jw, iw
    real(kind=Real8_) :: weight, esum, csum, psum, precfl
    real(kind=Real8_) :: XNNO(NR),XNDO(NR)
    real(kind=Real8_) :: F(NR,NT,NE,NPA),NSUM,FZERO(NR,NT,NE),ENO(NR),EDO(NR),AVEFL(NR,NT,NE),BARFL(NE)
    character(len=23) :: StringDate
    character(len=2)  :: ST2
    character(len=100) :: ST4, NameFileOut
    character(len=2), dimension(4) :: speciesString = (/'_e','_h','he','_o'/)

    ST2 = speciesString(S)
    ST4 =trim(PathRamOut)
    write(StringDate,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2)'.',i3.3)"), TimeRamNow%iYear, &
          TimeRamNow%iMonth, TimeRamNow%iDay, TimeRamNow%iHour, TimeRamNow%iMinute, &
          TimeRamNow%iSecond, floor(TimeRamNow%FracSecond*1000.0)

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

22  FORMAT(F5.2,F10.6,E13.4)
30  FORMAT(F7.2,72(1PE11.3))
31  FORMAT(F6.2,1X,F6.2,2X,F6.2,1X,F7.2,1X,F7.2,1X,1PE11.3)
32  FORMAT(' EKEV/PA, Date=',a,' L=',F6.2,' Kp=',F6.2,' MLT=',F4.1)
40  FORMAT(20(3X,F8.2))
70  FORMAT(F5.2,F10.6,E13.4)
71  FORMAT(2X,3HT =,F8.0,2X,4HKp =,F6.2,2X,'  Total Precip Flux [1/cm2/s]')
96  FORMAT(2HT=,F6.2,4H Kp=,F5.2,4H AP=,F7.2,4H Rs=,F7.2,7H Date= ,A23,'Plasmasphere e- density [cm-3]')

    RETURN
  END


!===========================================================================
  subroutine ram_write_pressure(StringIter)
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

    implicit none

    character(len=23)            :: StringTime
    character(len=*), intent(in) :: StringIter
    character(len=*), parameter  :: NameSub = 'ram_write_pressure'
    character(len=200)           :: FileName
    integer                      :: iError, i, j
    !------------------------------------------------------------------------
    ! Create Ram pressure output file.
    FileName = trim(PathRamOut)//trim(RamFileName('/pressure','dat',TimeRamNow))
    open(unit=UNITTMP_, status='REPLACE', iostat=iError, file=FileName)
    if(iError /= 0) call CON_stop(NameSub//' Error opening file '//FileName)

    ! Prepare date/time for header.
    write(StringTime,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2)'.',i3.3)") &
             TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
             TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond, &
             floor(TimeRamNow%FracSecond*1000.0)

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
             write(UNITTMP_,*) LZ(I),PHI(J)*12/PI,PPERH(I,J),PPARH(I,J), &
                               PPERO(I,J),PPARO(I,J), PPERHE(I,J),PPARHE(I,J), &
                               PPERE(I,J),PPARE(I,J),PAllSum(i,j)
          else
             write(UNITTMP_,*) LZ(I),PHI(J)*12/PI,PPERH(I,J),PPARH(I,J), &
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
  use ModRamTiming,    ONLY: TimeRamNow
  use ModRamVariables, ONLY: LZ, MLT, FNHS, BOUNHS, FNIS, BOUNIS, BNES, HDNS

  use ModIOUnit, ONLY: UNITTMP_

  implicit none

  integer :: i, j, L
  character(len=200) :: filenamehI

  filenamehI=trim(PathScbOut)//RamFileName('/hI_output', 'dat', TimeRamNow)
  PRINT*, 'Writing to file ', TRIM(filenamehI)
  OPEN(UNITTMP_, file = TRIM(filenamehI), action = 'write', status = 'unknown')
  WRITE(UNITTMP_, *) '   Lsh        MLT    NPA    h(Lsh,MLT,NPA) hBoun(Lsh,MLT,NPA) '//&
             'I(Lsh,MLT,NPA) IBoun(Lsh,MLT,NPA) Bz(Lsh,MLT) HDENS(Lsh,MLT,NPA)'
  DO i = 2, nR+1
     DO j = 1, nT
        DO L = 1, NPA
           WRITE(UNITTMP_, *) LZ(i), NINT(MLT(j)), L, FNHS(i,j,L), BOUNHS(i,j,L), &
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

  character(len=200) :: NameFluxFile

  NameFluxFile=trim(PathRamOut)//RamFileName('Dsbnd/ds'//St2(S),'dat',TimeRamNow)
  OPEN(UNIT=UNITTMP_,FILE=NameFluxFile, STATUS='UNKNOWN')
  WRITE(UNITTMP_,*)'EKEV FGEOSB [1/cm2/s/ster/keV] T=',TimeRamElapsed/3600,Kp,F107
  DO K=2,NE
     WRITE(UNITTMP_,*) EKEV(K),(FGEOS(S,J,K,2)/FFACTOR(S,NR,K,2),J=1,NT)
  END DO
  CLOSE(UNITTMP_)

end Subroutine write_dsbnd

!==============================================================================

  subroutine write_fail_file

    !!!! Module Variables
    use ModRamMain,      ONLY: PathRestartOut, PathRestartIn, niter
    use ModRamFunctions, ONLY: RamFileName
    use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamStart, TimeRamNow, DtsNext, &
                               TOld
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP, FGEOS, PParH, PPerH, PParO, &
                               PPerO, PParHe, PPerHe, PParE, PPerE
    use ModRamScb,       ONLY: indexPA, FLUX3DEQ
    use ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap
    use ModScbVariables, ONLY: x, y, z, bX, bY, bZ, bf, alfa, psi, alphaVal, psiVal, &
                               chi, GradRhoSq, GradZetaSq, GradRhoGradZeta, &
                               GradRhoGradTheta, GradThetaGradZeta, dPPerdRho, &
                               dPPerdZeta, dPPerdTheta, dBsqdRho, dBsqdZeta, &
                               dBsqdTheta, sigma, f, fzet
    !!!! Module Subroutines/Functions
    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts
    !!!! Share Modules
    use ModIOUnit, ONLY: UNITTMP_
    !!!! NetCdf Modules
    use netcdf

    implicit none
    
    integer :: stat
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iPParTVar, &
               iPPerTVar, iBxVar, iByVar, iBzVar, iBfVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iGEOVar, iPaIndexVar, iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAlphaVar, iBetaVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, &
               iChiVar, iTOldVar
    integer :: igrgr, igpgp, igrgp, igrgt, igtgp, idpdr, idpdp, idpdt, &
               idbdr, idbdp, idbdt, isigma, ifVar, ifzetVar
    integer :: nRDim, nTDim, nEDim, nPaDim, nSDim, nThetaDim, nPsiDim, &
               nZetaDim, nRPDim, nZetaPDim, iFlux3DVar
    integer, parameter :: iDeflate = 2

    character(len=2), dimension(4):: NameSpecies = (/'e_','h_','he','o_'/)
    character(len=200)            :: NameFile,CWD

    character(len=*), parameter :: NameSub='write_fail_file'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    ! OPEN FILE
    !NameFile = RamFileName(PathRestartOut//'/restart','nc',TimeRamNow)
    NameFile = 'failed.nc'
    iStatus = nf90_create(trim(NameFile), nf90_clobber, iFileID)
    !iStatus = nf90_create(trim(NameFile), nf90_HDF5, iFileID)
    call ncdf_check(iStatus, NameSub)
    call write_ncdf_globatts(iFileID)

    ! CREATE DIMENSIONS
    iStatus = nf90_def_dim(iFileID, 'nR',     nR,     nRDim)
    iStatus = nf90_def_dim(iFileID, 'nRP',    nR+1,   nRPDim)
    iStatus = nf90_def_dim(iFileID, 'nT',     nT,     nTDim)
    iStatus = nf90_def_dim(iFileID, 'nE',     nE,     nEDim)
    iStatus = nf90_def_dim(iFileID, 'nPa',    nPa,    nPaDim)
    iStatus = nf90_def_dim(iFileID, 'nS',     4,      nSDim)
    iStatus = nf90_def_dim(iFileID, 'nTheta', nthe,   nThetaDim)
    iStatus = nf90_def_dim(iFileID, 'nPsi',   npsi,   nPsiDim)
    iStatus = nf90_def_dim(iFileID, 'nZeta',  nzeta,  nZetaDim)
    iStatus = nf90_def_dim(iFileID, 'nZetaP', nzetap, nZetaPDim)

    ! START DEFINE MODE
    !! FLUXES
    iStatus = nf90_def_var(iFileID, 'FluxE', nf90_double,(/nRdim,nTDim,nEDim,nPaDim/), iFluxEVar)
    iStatus = nf90_def_var(iFileID, 'FluxH', nf90_double,(/nRdim,nTDim,nEDim,nPaDim/), iFluxHVar)
    iStatus = nf90_def_var(iFileID, 'FluxO', nf90_double,(/nRdim,nTDim,nEDim,nPaDim/), iFluxOVar)
    iStatus = nf90_def_var(iFileID, 'FluxHe', nf90_double,(/nRdim,nTDim,nEDim,nPaDim/), iFluxHeVar)
    !iStatus = nf90_def_var(iFileID, 'FGEOS', nf90_double,(/nSdim,nTDim,nEDim,nPadim/), iGEOVar)
    !iStatus = nf90_def_var(iFileID, 'Flux3D', nf90_double,(/nSdim,nPsiDim,nZetaDim,nEDim,nPADim/), iFlux3DVar)

    !! PRESSURES
    iStatus = nf90_def_var(iFileID, 'PParT', nf90_double,(/nSdim,nRDim,nTDim/), iPParTVar)
    iStatus = nf90_def_var(iFileID, 'PPerT', nf90_double,(/nSdim,nRDim,nTDim/), iPPerTVar)

    !! MAGNETIC FIELD
    iStatus = nf90_def_var(iFileID, 'Bx', nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), iBxVar)
    iStatus = nf90_def_var(iFileID, 'By', nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), iByVar)
    iStatus = nf90_def_var(iFileID, 'Bz', nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), iBzVar)
    iStatus = nf90_def_var(iFileID, 'B', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iBfVar)

    !! ELECTRIC FIELD/POTENTIAL
    !iStatus = nf90_def_var(iFileID, 'VTOL', nf90_double,(/nRPDim,nTDim/), iVTOLVar)
    !iStatus = nf90_def_var(iFileID, 'VTN', nf90_double,(/nRPDim,nTDim/), iVTNVar)
    !iStatus = nf90_def_var(iFileID, 'VT', nf90_double,(/nRPDim,nTDim/), iVTVar)

    !! hI OUTPUTS
    !iStatus = nf90_def_var(iFileID, 'FNHS', nf90_double,(/nRPDim,nTDim,NPaDim/), iHVar)
    !iStatus = nf90_def_var(iFileID, 'BOUNHS', nf90_double,(/nRPDim,nTDim,NPaDim/), iBHVar)
    !iStatus = nf90_def_var(iFileID, 'FNIS', nf90_double,(/nRPDim,nTDim,NPaDim/), iIVar)
    !iStatus = nf90_def_Var(iFileID, 'BOUNIS', nf90_double,(/nRPDim,nTDim,NPaDim/), iBIVar)
    !iStatus = nf90_def_var(iFileID, 'BNES', nf90_double,(/nrPDim,nTDim/), iBNESVar)
    !iStatus = nf90_def_var(iFileID, 'HDNS', nf90_double,(/nRPDim,nTDim,nPaDim/), iHDNSVar)
    !iStatus = nf90_def_var(iFileID, 'EIR', nf90_double,(/nrPDim,nTDim/), iEIRVar)
    !iStatus = nf90_def_var(iFileID, 'EIP', nf90_double,(/nrPDim,nTDim/), iEIPVar)
    !iStatus = nf90_def_var(iFileID, 'dBdt', nf90_double,(/nrPDim,nTDim/), iDBDTVar)
    !iStatus = nf90_def_var(iFileID, 'dIdt', nf90_double,(/nRPDim,nTDim,NPaDim/), iDIDTVar)
    !iStatus = nf90_def_var(iFileID, 'dIbndt', nf90_double,(/nRPDim,nTDim,NPaDim/), iDIBNVar)

    !! GRID OUTPUTS
    iStatus = nf90_def_var(iFileID, 'x', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iXVar)
    iStatus = nf90_def_var(iFileID, 'y', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iYVar)
    iStatus = nf90_def_var(iFileID, 'z', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iZVar)

    !! ALPHA/BETA
    iStatus = nf90_def_var(iFileID, 'alpha', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iAlphaVar)
    iStatus = nf90_def_var(iFileID, 'beta', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iBetaVar)
    iStatus = nf90_def_var(iFileID, 'chi', nf90_double,(/nThetaDim,nPsiDim,nZetaPDim/), iChiVar)

    !! Gradients/Derivatives
    iStatus = nf90_def_var(iFileID, 'GradRhoSq',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), igrgr)
    iStatus = nf90_def_var(iFileID, 'GradPhiSq',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), igpgp)
    iStatus = nf90_def_var(iFileID, 'GradRhoGradPhi',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), igrgp)
    iStatus = nf90_def_var(iFileID, 'GradRhoGradTheta',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), igrgt)
    iStatus = nf90_def_var(iFileID, 'GradThetaGradPhi',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), igtgp)
    iStatus = nf90_def_var(iFileID, 'dPPerdRho',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idpdr)
    iStatus = nf90_def_var(iFileID, 'dPPerdPhi',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idpdp)
    iStatus = nf90_def_var(iFileID, 'dPPerdTheta',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idpdt)
    iStatus = nf90_def_var(iFileID, 'dBsqdRho',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idbdr)
    iStatus = nf90_def_var(iFileID, 'dBsqdPhi',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idbdp)
    iStatus = nf90_def_var(iFileID, 'dBsqdTheta',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), idbdt)
    iStatus = nf90_def_var(iFileID, 'sigma',nf90_double,(/nThetaDim,nPsiDim,nZetaDim/), isigma)
    iStatus = nf90_def_var(iFileID, 'f',nf90_double,(/nPsiDim/), ifVar)
    iStatus = nf90_def_var(iFileID, 'fzet',nf90_double,(/nZetaPDim/), ifzetVar)

    ! END DEFINE MODE
    iStatus = nf90_enddef(iFileID)
    call ncdf_check(iStatus, NameSub)

    ! START WRITE MODE
    !! FLUXES
    iStatus = nf90_put_var(iFileID, iFluxEVar,  F2(1,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHVar,  F2(2,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHeVar, F2(3,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxOVar,  F2(4,:,:,:,:))
    !iStatus = nf90_put_var(iFileID, iGEOVar,    FGEOS(:,:,:,:))
    !iStatus = nf90_put_var(iFileID, iFlux3DVar, FLUX3DEQ(:,:,:,:,:))

    !! PRESSURES
    iStatus = nf90_put_var(iFileID, iPParTVar, PParT(:,:,:))
    iStatus = nf90_put_var(iFileID, iPPerTVar, PPerT(:,:,:))

    !! MAGNETIC FIELD
    iStatus = nf90_put_var(iFileID, iBxVar, bX(:,:,:))
    iStatus = nf90_put_var(iFileID, iByVar, bY(:,:,:))
    iStatus = nf90_put_var(iFileID, iBzVar, bZ(:,:,:))
    iStatus = nf90_put_var(iFileID, iBfVar, bf(:,:,:))

    !! ELECTRIC FIELD
    !iStatus = nf90_put_var(iFileID, iVTOLVar, VTOL(:,:))
    !iStatus = nf90_put_var(iFileID, iVTNVar,  VTN(:,:))
    !iStatus = nf90_put_var(iFileID, iVTVar,   VT(:,:))

    !! hI OUTPUTS
    !iStatus = nf90_put_var(iFileID, iHVar,    FNHS(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBHVar,   BOUNHS(:,:,:))
    !iStatus = nf90_put_var(iFileID, iIVar,    FNIS(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBIVar,   BOUNIS(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBNESVar, BNES(:,:))
    !iStatus = nf90_put_var(iFileID, iHDNSVar, HDNS(:,:,:))
    !iStatus = nf90_put_var(iFileID, iEIRVar,  EIR(:,:))
    !iStatus = nf90_put_var(iFileID, iEIPVar,  EIP(:,:))
    !iStatus = nf90_put_var(iFileID, iDBDTVar, dBdt(:,:))
    !iStatus = nf90_put_var(iFileID, iDIDTVar, dIdt(:,:,:))
    !iStatus = nf90_put_var(iFileID, iDIBNVar, dIbndt(:,:,:))

    !! GRID OUTPUTS
    iStatus = nf90_put_var(iFileID, iXVar, x(:,:,:))
    iStatus = nf90_put_var(iFileID, iYVar, y(:,:,:))
    iStatus = nf90_put_var(iFileID, iZVar, z(:,:,:))

    !! ALPHA/BETA
    iStatus = nf90_put_var(iFileID, iAlphaVar, psi(:,:,:))
    iStatus = nf90_put_var(iFileID, iBetaVar,  alfa(:,:,:))
    iStatus = nf90_put_var(iFileID, iChiVar,   chi(:,:,:))

    !! Gradients/Derivatives
    iStatus = nf90_put_var(iFileID, igrgr, GradRhoSq(:,:,:))
    iStatus = nf90_put_var(iFileID, igpgp, GradZetaSq(:,:,:))
    iStatus = nf90_put_var(iFileID, igrgp, GradRhoGradZeta(:,:,:))
    iStatus = nf90_put_var(iFileID, igrgt, GradRhoGradTheta(:,:,:))
    iStatus = nf90_put_var(iFileID, igtgp, GradThetaGradZeta(:,:,:))
    iStatus = nf90_put_var(iFileID, idpdr, dPPerdRho(:,:,:))
    iStatus = nf90_put_var(iFileID, idpdp, dPPerdZeta(:,:,:))
    iStatus = nf90_put_var(iFileID, idpdt, dPPerdTheta(:,:,:))
    iStatus = nf90_put_var(iFileID, idbdr, dBsqdRho(:,:,:))
    iStatus = nf90_put_var(iFileID, idbdp, dBsqdZeta(:,:,:))
    iStatus = nf90_put_var(iFileID, idbdt, dBsqdTheta(:,:,:))
    iStatus = nf90_put_var(iFileID, isigma, sigma(:,:,:))
    iStatus = nf90_put_var(iFileID, ifVar, f(:))
    iStatus = nf90_put_var(iFileID, ifzetVar, fzet(:))

    ! END WRITE MODE
    call ncdf_check(iStatus, NameSub)
    
    ! CLOSE FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_fail_file
!==============================================================================

END MODULE ModRamIO
