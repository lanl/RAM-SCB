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
    use ModRamTiming,    ONLY: Dt_hI, DtRestart, DtWriteSat, TimeRamNow, TimeRamElapsed
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
    character(len=4) :: StrScbIter
    character(len=*), parameter :: NameSub = 'handle_output'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Get current Dst.
    call get_ramdst(dst)

    ! Write Logfile
    if(mod(TimeIn, DtLogfile)==0.0)then
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

       write(StrScbIter,'(I4.4)') int(TimeRamElapsed/Dt_hI)
       call ram_sum_pressure
       call ram_write_pressure(StrScbIter)
    end if

    if (DoSaveRamSats .and. (mod(TimeRamElapsed,DtWriteSat) .eq. 0)) call fly_sats

    ! Write restarts.
    if (mod(TimeIn,DtRestart).eq.0) call write_restart()

  end subroutine handle_output

!===========================!
!==== INPUT SUBROUTINES ====!
!===========================!
!==============================================================================
subroutine read_geomlt_file(NameParticle)

  !!!! Module Variables
  use ModRamMain,      ONLY: PathRamIN
  use ModRamTiming,    ONLY: TimeRamNow, Dt_bc
  use ModRamGrids,     ONLY: NE, NEL, NEL_prot, NBD
  use ModRamParams,    ONLY: boundary
  use ModRamVariables, ONLY: EKEV, FGEOS, IsInitialized, timeOffset, StringFileDate, &
                             flux_SIII, fluxLast_SII, eGrid_SI, avgSats_SI, lGrid_SI
  !!!! Share Modules
  use ModIOUnit, ONLY: UNITTMP_

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

  ! First allocate global arrays

  if(.not.allocated(flux_SIII)) then

     ! allocate larger of two sizes
     if(boundary .eq. 'PTM') then
        NEL_ = max(NEL,NEL_prot)
     else
        NEL_ = NEL
     endif

     allocate(flux_SIII(2,NBD,24,NEL_),fluxLast_SII(2,24,NEL_),eGrid_SI(2,NEL_),avgSats_SI(2,NBD))
     flux_SIII = 0.
     fluxLast_SII = 0.
     eGrid_SI = 0.
  endif

  ! Possible different sizes for two boundary input files
  if((NameParticle .eq. 'prot').and.(boundary .eq. 'PTM')) then
     NEL_ = NEL_prot
  else
     NEL_ = NEL
  endif

  ! Now local arrays
  allocate(Buffer_I(NEL_),Buffer_III(NBD,24,NEL_),iBuffer_II(24,NEL_),Buffer2_I(NBD))

    !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Build file name using current date.  NameParticle decides if we open
  ! electrons or protons.
  if( (NameParticle.ne.'prot') .and. (NameParticle.ne.'elec') ) &
       call CON_stop(NameSub//': Invalid particle type '//&
       '"'//NameParticle//'" (can be prot or elec)')
  ! Chris Jeffery mod, 3/11
  if(Dt_bc.gt.60) then
     write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a)') trim(PathRamIn)//'/', &
          TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
          '_mpa-sopa_', NameParticle, '_geomlt_5-min.txt'
  else
     write(NameFileIn, '(a,i4.4,i2.2,i2.2,3a,i2.2,a)') trim(PathRamIn)//'/', &
          TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
          '_mpa-sopa_', NameParticle, '_geomlt_', Int(Dt_bc), &
          '-sec.txt'
  endif

  if(DoTest) write(*,*)NameSub//': Opening '//trim(NameFileIn)

  ! Open file, check status.
  inquire(unit=UNITTMP_,opened=IsOpened)
  if(IsOpened)write(*,*) 'WARNING- File already opened.'
  open(unit=UNITTMP_, file=NameFileIn, status='OLD', iostat=iError)
  if(iError/=0) call CON_stop(NameSub//': Error opening file '//NameFileIn)

  ! Save file date.
  write(StringFileDate,'(i4.4,i2.2,i2.2)') &
       TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay

  ! Set species index.
  iSpec=1 !Protons.
  if(NameParticle.eq.'elec') iSpec=2 !electrons.

  ! Skip header.
  read(UNITTMP_,*) StringHeader
  read(UNITTMP_,*) StringHeader
  read(UNITTMP_,*) StringHeader

  ! Load energy grid.
  write(NEL_char,*) NEL_
  read(UNITTMP_, '(104x,'//trim(adjustl(NEL_char))//'(G18.6))') Buffer_I

  ! Check energy grid against RAM energy grid.
  if ( (Buffer_I(1)>Ekev(1)).or.(Buffer_I(NEL_)<Ekev(nE)) ) then
     write(*,*)NameSub//': LANL Geo energy grid does not cover RAM energy grid'
     write(*,*)'  LANL E       RAM E'
     write(*,'(2(ES12.6,1x))') Buffer_I(1),  Ekev(1)
     write(*,'(2(ES12.6,1x))') Buffer_I(NEL_), Ekev(nE)
     call CON_stop(NameSub//': Energy grid error.')
  end if

  ! Load fluxes.
  do i=1,NBD
     do j=1, 24
        ! Read one time entry worth of data.
        read(UNITTMP_,'(24x,f6.1,2x,'//trim(adjustl(NEL_CHAR))//'(i2),'//trim(adjustl(NEL_CHAR))//'(f18.4))')&
            lGrid_SI(iSpec, j), iBuffer_II(j,:), Buffer_III(i,j,:)
!        read(UNITTMP_,'(10x,f7.2,4x,36e12.4))')&
!             lGrid_SI(iSpec, j), Buffer_III(i,j,:)
     end do
     ! Crunch number of satellites down to one number.
     Buffer2_I(i) = real(sum(iBuffer_II))/(1.0*NEL_)
!     print*,'i, Buffer_III(i,1,1)=',i,Buffer_III(i,1,1)
  end do

  close(UNITTMP_)

  ! Wrap grid about LT=0.
  lGrid_SI(iSpec, 0) = 2.*lGrid_SI(iSpec, 1)-lGrid_SI(iSpec, 2)
  lGrid_SI(iSpec,25) = 2.*lGrid_SI(iSpec,24)-lGrid_SI(iSpec,23)

  if(DoTest)then
     write(*,*)NameSub//': File read successfully.'
     write(*,*)NameSub//': Energy Grid='
     write(*,'(3(ES12.6, 1x))') Buffer_I
     write(*,*)NameSub//': LT Grid='
     write(*,'(2(F6.2, 1x))') lGrid_SI(iSpec,:)
     write(*,*)NameSub//': First data line='
     write(*,'(6(ES12.5, 1x))')Buffer_III(1,1,:)
     write(*,*)NameSub//': Last data line='
     write(*,'(6(ES12.5, 1x))')Buffer_III(NBD,24,:)
  end if

  ! Put data into proper array according to NameParticle.
  eGrid_SI(iSpec,1:NEL_)      = Buffer_I
  flux_SIII(iSpec,:,:,1:NEL_) = Buffer_III
  avgSats_SI(iSpec,:)    = Buffer2_I

  ! Scrub for bad data.  Must do that here because we don't know
  ! before hand were we will start inside of the file.
  do j=1,24; do k=1,NEL_
     ! If the first entry is bad, use fluxLast_SII.
     if(flux_SIII(iSpec,1,j,k).le.0) &
          flux_SIII(iSpec,1,j,k)= fluxLast_SII(iSpec,j,k)
     ! If later entries are bad, revert to last good data point.
     do i=2,NBD
        if(flux_SIII(iSpec,i,j,k).le.0) &
             flux_SIII(iSpec,i,j,k)=flux_SIII(iSpec,i-1,j,k)
     end do
  end do; end do

  ! Finally, store last entry into fluxLast_SII.
  fluxLast_SII(iSpec,:,:) = flux_SIII(iSpec,NBD,:,:)

  deallocate(Buffer_I,Buffer_III,iBuffer_II,Buffer2_I)

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

!============================!
!==== OUTPUT SUBROUTINES ====!
!============================!
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
END MODULE ModRamIO
