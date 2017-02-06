!============================================================================
subroutine ram_gen_efilename(TimeIn,NameOut)
!    Generate the correct filname given TimeIn, start time TimeRamStart,
!    and electric field model Electric.  Returned is a 100-character character
!    with the full file name.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

  use ModTimeConvert
  use ModRamIO
  use ModRamMain, ONLY: event, electric, TimeRamStart, PathRamIn, PathScbOut

  implicit none

  ! Args and return value.
  type(TimeType), intent(in) :: TimeIn
  character(len=200), intent(out):: NameOut

  integer :: nDay, nFile, nFive
  character(len=100) :: NameSubdir
  character(len=*), parameter :: NameSub = 'ram_gen_efilename'
  !--------------------------------------------------------------------------
  ! Get old ram style characters for building file names.
  ! These formats are a candidate for removal.
  call gen_old_timetags(TimeRamStart, TimeIn, nDay=nDay, &
       nFiveDay=nFive, nFiveTotal=nFile)

  write(NameSubdir, "('/Day',i2.2,'Elfld')") nDay

  ! Generate file name based on time of simulation and file type.
  select case(trim(electric))

  case('WESC')
     ! use weimer 2001  along sc field lines.
     RETURN

  case('W5SC')
     ! use weimer 2005 along sc field lines.
     RETURN

  case('IESC') !SWMF IE potential traced to equitorial plane using SCB.
     if(UseNewFmt) then
        NameOut = trim(PathScbOut) // &
             RamFileName('/IE_SCB','in',TimeIn)
     else
        ! Old way uses event and number of five-minute intevals over
        ! entire run.
        write(NameOut,"(3a,'_',i4.4,a)") &
             trim(PathScbOut),'/IE_SCB_', trim(event), nFile, '_RAM.in'
     end if

  case('IE89') !SWMF IE potential traced to equitorial plane using T89.
     if(useNewFmt) then
        NameOut = trim(PathRamIn) // &
             RamFileName('/IE','in',TimeIn)
     else
        ! Old way uses date, time, and nDay.
        write(NameOut,"(a,'/',a,i4.4,2i2.2,'_',3i2.2,a)") &
             trim(PathRamIn), 'IE_', &
             TimeIn % iYear, TimeIn % iMonth, TimeIn % iDay, &
             TimeIn % iHour, TimeIn % iMinute, TimeIn % iSecond, &
             '_RAM.in'
     endif

  case('WE01') !Weimer 2001 electric potential.
     if(useNewFmt) then
        NameOut = trim(PathRamIn) // &
             RamFileName('/weq01','in',TimeIn)
     else
        ! Old way uses nDay and nFiveDay.
        write(NameOut,"(a,'/',a,i4.4,a)") &
             trim(PathRamIn), 'weq01_', &
             nFive+1, '.in'
     endif

  case('VOLS') !Volland-Stern electric potential.
     if(useNewFmt) then
        NameOut = trim(PathRamIn) // &
             RamFileName('/vsc02','in',TimeIn)
     else
        ! Old way uses nDay, and nFiveDay
        write(NameOut,"(a,'/',a,i4.4,a)") &
             trim(PathRamIn), 'vsc02_', &
             nFive+1, '.in'
     end if
  case default
     call CON_stop(NameSub//' ERROR: Unsupported electric field file type')

  end select

end subroutine ram_gen_efilename

!============================================================================
subroutine ram_get_electric(NextEfile, EOut_II)
  ! Open E-field file "NextFile", collect E-field and indices, update
  ! "NextFile" and return all of the above to the caller.

  use ModTimeConvert
  use ModRamIO
  use ModIoUnit,    ONLY: UNITTMP_
  use ModRamMain,   ONLY: Real8_, IsComponent, PathRamIn, PathScbOut, &
       nR, nT, electric, TimeRamElapsed

  implicit none

  ! Arguments
  character(len=200),intent(in)  :: NextEfile
  real(kind=Real8_), intent(out) :: EOut_II(NR+1, NT)
  ! Kp, F10.7 no longer acquired through E-field files.
  real(kind=Real8_):: KpOut, f107Out
  

  type(TimeType) :: TimeNext
  integer :: nDay, nFile, iError, i, j, jw
  character(len=200) :: StringHeader
  real(kind=Real8_) :: day, tth, ap, rsun, RRL, PH, wep(49), phiofs

  character(len=*), parameter :: NameSub = 'ram_get_electric'
  !--------------------------------------------------------------------------
  ! Initialize efield to zero.
  ! NOTE THAT ON RESTART, WE MUST CHANGE THIS.
  EOut_II  = 0.0
  
  IF ((electric .EQ. 'WESC') .or. IsComponent .or. (electric .eq. 'W5SC')) THEN
     ! print*, 'Reading electric potentials from SCB directly'
     RETURN
  END IF

  call write_prefix
  write(*,*) 'Reading EFile ', trim(NextEFile)
  write(*,'(a,f10.1)')'          at Sim Time', TimeRamElapsed

  !\
  ! Read first file and indices.
  !/
  open(file=NextEfile, unit=UNITTMP_, status='OLD', &
       action='READ', IOSTAT=ierror)

  if(ierror .ne. 0) &
       call CON_stop(NameSub//': ERROR opening file '//trim(NextEfile))

  ! Read header, get Kp and F107 (skip the other crap)
  read(UNITTMP_,'(a)') StringHeader
  read(UNITTMP_,*) day, tth, KpOut, f107Out, ap, rsun, phiofs
  read(UNITTMP_,'(a)') StringHeader
  read(UNITTMP_,'(a)') StringHeader

  ! If coupled to SWMF IE, first E-field is zero.
  if(IsComponent .and. electric .eq. 'IESC') then 
     close(UNITTMP_)
     return
  endif

  do i=1,nR+1
     if (electric.EQ.'IESC') then
        do JW=1,NT
           read(UNITTMP_,*) RRL,PH,WEP(JW) ! conv potential in kV
        enddo
     else
        do JW=1,49
           read(UNITTMP_,*) RRL,PH,WEP(JW) ! conv potential in kV
        enddo
     end if
     do J=1,NT
        if(NT.EQ.49 .OR. electric.EQ.'IESC') JW=J
        if(NT.EQ.25 .AND. electric.NE.'IESC') JW=2*J-1
        EOut_II(I,J)=WEP(JW)*1e3			! to convert in [V]
     enddo
  enddo

  close(UNITTMP_)

  do J=1,NT
     do I=NR,1,-1
        if (EOut_II(I,J).EQ.0) EOut_II(I,J) = EOut_II(I+1,J)/2
     enddo
  enddo

end subroutine ram_get_electric

!============================================================================
