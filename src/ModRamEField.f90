MODULE ModRamEField
! Contains subroutines related to getting the electric field for RAM

!  use ModRamMain, ONLY: Real8_
!  use ModRamGrids, ONLY: NR, NT

!  implicit none

!  real(kind=Real8_), dimension(NR+1,NT) :: VT, EIR, EIP, VTOL, VTN
!  real(kind=Real8_) :: TOLV

  use ModRamVariables, ONLY: VT, VTOL, VTN, TOLV, EIR, EIP

  implicit none
  save

  contains

!==============================================================================
subroutine get_electric_field

  use ModRamMain,    ONLY: Real8_
  use ModRamTiming,  ONLY: TimeRamNow, DtEfi, T
  use ModRamConst,   ONLY: RE
  use ModRamParams,  ONLY: electric, IsComponent
  use ModRamGrids,   ONLY: NR, NT
!  use ModRamIndices, ONLY: KP
!  use ModRamInit,    ONLY: PHI, LZ, PHIOFS
  use ModRamVariables, ONLY: KP, PHI, LZ, PHIOFS

  use ModRamCouple,  ONLY: SwmfPot_II

  use ModTimeConvert, ONLY: TimeType, time_real_to_int

  implicit none

  real(kind=Real8_) :: AVS
  integer :: I, J
  type(TimeType) :: TimeNext
  character(len=200) :: NameEfile

  ! Set "new" values of indices and E-field to "old" values. 
  VTOL = VTN
  print*,'RAM: updating E field at time ', T/3600.
  ! Generate name of next electric field file:
  TimeNext = TimeRamNow
  if (electric.ne.'IESC') then
     TimeNext % Time = TimeNext % Time + DtEfi
  else
     TimeNext % Time = TimeNext % Time
  end if
  call time_real_to_int(TimeNext)
  call ram_gen_efilename(TimeNext, NameEfile)
  ! Open this file and save contents.
  if(electric .ne. 'VOLS') call ram_get_electric(NameEfile, vtn)
  ! No interpolation for IESC case.
  if(electric .eq. 'IESC') vtol=vtn
  ! Set time of "old" file for interpolation:
  TOLV=T
  if (IsComponent .OR. electric=='WESC' .or. electric=='W5SC') then
     VTOL = SwmfPot_II
     VTN  = SwmfPot_II
  endif

  DO I=1,NR+1
     DO J=1,NT
        if (electric .ne. 'VOLS') then
           VT(I,J)=VTOL(I,J)+(VTN(I,J)-VTOL(I,J))*(T-TOLV)/DtEfi
        else
           AVS=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE  ! Voll-Stern parameter
           VT(I,J)=AVS*(LZ(I)*RE)**2*SIN(PHI(J)-PHIOFS) ! [V]
        endif
     ENDDO
  ENDDO

end subroutine

!============================================================================
subroutine ram_gen_efilename(TimeIn,NameOut)
!    Generate the correct filname given TimeIn, start time TimeRamStart,
!    and electric field model Electric.  Returned is a 100-character character
!    with the full file name.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

  use ModRamFunctions, ONLY: RamFileName
  use ModRamMain,   ONLY: PathRamIn
  use ModScbMain,   ONLY: PathScbOut
  use ModRamParams, ONLY: electric

  use ModTimeConvert, ONLY: TimeType

  implicit none

  ! Args and return value.
  type(TimeType), intent(in) :: TimeIn
  character(len=200), intent(out):: NameOut

  integer :: nDay, nFile, nFive
  character(len=*), parameter :: NameSub = 'ram_gen_efilename'
  !--------------------------------------------------------------------------

  ! Generate file name based on time of simulation and file type.
  select case(trim(electric))

  case('WESC')
     ! use weimer 2001  along sc field lines.
     RETURN

  case('W5SC')
     ! use weimer 2005 along sc field lines.
     RETURN

  case('IESC') !SWMF IE potential traced to equitorial plane using SCB.
     NameOut = trim(PathScbOut)//RamFileName('/IE_SCB','in',TimeIn)

  case('IE89') !SWMF IE potential traced to equitorial plane using T89.
     NameOut = trim(PathRamIn)//RamFileName('/IE','in',TimeIn)

  case('WE01') !Weimer 2001 electric potential.
     NameOut = trim(PathRamIn)//RamFileName('/weq01','in',TimeIn)

  case('VOLS') !Volland-Stern electric potential.
     NameOut = trim(PathRamIn)//RamFileName('/vsc02','in',TimeIn)
  
  case default
     call CON_stop(NameSub//' ERROR: Unsupported electric field file type')

  end select

end subroutine ram_gen_efilename

!============================================================================
subroutine ram_get_electric(NextEfile, EOut_II)
  ! Open E-field file "NextFile", collect E-field and indices, update
  ! "NextFile" and return all of the above to the caller.

  use ModRamMain,   ONLY: Real8_, PathRamIn, PathScbOut
!  use ModScbMain,   ONLY: PathScbOut
  Use ModRamTiming, ONLY: TimeRamElapsed
  use ModRamParams, ONLY: electric, IsComponent
  use ModRamGrids,  ONLY: NR, NT

  use ModRamVariables, ONLY: PHIOFS

  use ModTimeConvert, ONLY: TimeType
  use ModIOUnit,      ONLY: UNITTMP_
  implicit none

  ! Arguments
  character(len=200),intent(in)  :: NextEfile
  real(kind=Real8_), intent(out) :: EOut_II(NR+1, NT)
  ! Kp, F10.7 no longer acquired through E-field files.
  real(kind=Real8_):: KpOut, f107Out
  

  type(TimeType) :: TimeNext
  integer :: nDay, nFile, iError, i, j, jw
  character(len=200) :: StringHeader
  real(kind=Real8_) :: day, tth, ap, rsun, RRL, PH, wep(49)

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
        EOut_II(I,J)=WEP(JW)*1e3 ! to convert in [V]
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

END MODULE ModRamEField
