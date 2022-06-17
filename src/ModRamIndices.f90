!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamIndices
  ! A set of tools for handling the indices used in RAM-SCB.  Combined with
  ! the PARAM interface, allows the user to select source of indices, gather
  ! and interpolate the values to the current simulation time.

  implicit none

  contains
  !===========================================================================
  subroutine read_index_file(StartTime, EndTime, NameFile)
    ! Read RAMIndices file

    use ModRamMain, ONLY: Real8_
    use ModRamVariables, ONLY: nRawKp, nRawF107, nRawAE, kptime, &
                               timeKp, timeF107, timeAE, rawKp, rawF107, rawAE
    use ModRamParams,   ONLY: DoUseEMIC
    use ModTimeConvert, ONLY: TimeType, time_int_to_real
    use ModIoUnit,      ONLY: UNITTMP_
    use ModUtilities,   ONLY: CON_stop, CON_set_do_test 


    implicit none

    type(timetype),   intent(in) :: StartTime, EndTime
    character(len=*), intent(in) :: NameFile

    integer :: dateIndex
    real(kind=Real8_) :: tmpF107

    integer :: i, j, iError, iYY, iMM, iDD
    character(len=100) :: StringLine, StringFmt
    real(kind=Real8_) :: tmpKp(8)

    integer :: cyy, cmm, cdd, chh, cmin, ij, ierr, iLines
    character(len=100) :: header, fname
    
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='read_index_file'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest) then
       write(*,*)'Loading indices from ', trim(NameFile)
       write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
            StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
            StartTime%iHour, StartTime%iMinute, StartTime%iSecond
       write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'End Time   = ', &
            EndTime%iYear, EndTime%iMonth, EndTime%iDay, &
            EndTime%iHour, EndTime%iMinute, EndTime%iSecond
    end if
    
    !!! Open file and find the starting date of file
    dateIndex = -1
    open(unit=UNITTMP_, FILE=NameFile, STATUS='OLD', IOSTAT=iError)
    if (iError.ne.0) call CON_stop(NameSub//' Error opening file '//NameFile)
    Read_RamIndices_Dates: Do
        read(UNITTMP_, '(i4, i2, i2, 8(1x,f3.1), 1x, f5.1)', IOSTAT=iError) &
                         iYY, iMM, iDD, tmpKp, tmpF107
        if (iError.lt.0) then
           call CON_stop( &
                NameSub//': Start date outside of range of RamIndices file')
        end if
        if ((StartTime%iYear  .eq. iYY)  .and. &
            (StartTime%iMonth .eq. iMM)  .and. &
            (StartTime%iDay   .eq. iDD)) then
           EXIT Read_RamIndices_Dates
        else
           dateIndex = dateIndex + 1
           CYCLE Read_RamIndices_Dates
        end if
    end do Read_RamIndices_Dates
    close(UNITTMP_)

    !!! Open file and fast forward to starting line
    open(unit=UNITTMP_, FILE=NameFile, STATUS='OLD', IOSTAT=iError)
    if (dateIndex.eq.0) then
      read(UNITTMP_, *) StringLine
      if (DoTest) then
        write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
                  StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
                  StartTime%iHour, StartTime%iMinute, StartTime%iSecond
        write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'End Time   = ', &
                  EndTime%iYear, EndTime%iMonth, EndTime%iDay, &
                  EndTime%iHour, EndTime%iMinute, EndTime%iSecond
        write(*,*) NameSub//': Line before requested date:'
        write(*,*) StringLine
      end if
    elseif (dateIndex.gt.0) then
      ! Fast forward to current position.
      write(StringFmt, "('(',i10,'(/)a)')") dateIndex
      read(UNITTMP_, trim(StringFmt)) StringLine
      if (DoTest) then
        write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
                  StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
                  StartTime%iHour, StartTime%iMinute, StartTime%iSecond
        write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'End Time   = ', &
                  EndTime%iYear, EndTime%iMonth, EndTime%iDay, &
                  EndTime%iHour, EndTime%iMinute, EndTime%iSecond
        write(*,*) NameSub//': Line before requested date:'
        write(*,*) StringLine
      end if
    end if

    !!! Allocate arrays based on number of days required.
    nRawF107 = ceiling((EndTime%Time-StartTime%Time)/86400.0)+1
    nRawKp   = 8*nRawF107
    allocate(timeKp(nRawKp))
    allocate(rawKp(nRawKp))
    allocate(rawF107(nRawF107))
    allocate(timeF107(nRawF107))
    timeKp = 0.0; rawKp = 0.0; rawF107 = 0.0; timeF107 = 0.0

    if(DoTest)then
       write(*,*) NameSub//': Number of lines to read = ', nRawF107
       write(*,*) NameSub//': Number of KP vals to read=', nRawKp
    end if

    !!! Read and store all data for entire interval.
    do i=1, nRawF107
       read(UNITTMP_, '(i4, i2, i2, 8(1x,f3.1), 1x, f5.1)',IOSTAT=iError) &
            iYY, iMM, iDD, tmpKp, rawF107(i)
       if (iError.lt.0) then
          call CON_stop( &
            NameSub//': End date outside of range of RamIndices file')
       elseif (iError.gt.0) then
          call CON_stop( &
            NameSub//': Error in formating of RamIndices file')
       end if
       call time_int_to_real((/iYY, iMM, iDD, 0,0,0,0/), timeF107(i))
       do j=1, 8
          rawKp(8*(i-1)+j) = tmpKp(j)
          call time_int_to_real((/iYY, iMM, iDD, kptime(j),30,0,0/),&
               timeKp(8*(i-1)+j))
       end do
    end do
    
    close(UNITTMP_)

    If (DoUseEMIC) then
       ! Read AE index from the file
       fname = 'AEindex.txt'
       open(UNITTMP_, FILE=fname, status = 'UNKNOWN', action='READ')
       ij = 0
       ierr = 0
       do while (ierr==0)
          read(UNITTMP_, *, iostat=ierr)header
          ij = ij + 1
       end do
       ilines = ij - 1
       rewind(UNITTMP_)
       
       nRawAE = iLines
       allocate(timeAE(nRawAE), rawAE(nRawAE))
       do ij=1, iLines
          read(UNITTMP_, '(i4,i2,i2,1x,i2,1x,i2,1x,i4)') &
               cyy, cmm, cdd, chh, cmin, rawAE(ij)
          call time_int_to_real((/cyy,cmm,cdd,chh,cmin,0,0/), timeAE(ij))
       end do
       
       close(UNITTMP_)
    end If

  end subroutine read_index_file

  !===========================================================================
  subroutine init_indices(StartTime, EndTime)
    ! Based on source of indices, prepare indices for this simulation.
    use ModRamParams,    ONLY: NameIndexFile
    use ModRamVariables, ONLY: NameIndexSource

    use ModTimeConvert, ONLY: TimeType
    use ModUtilities,   ONLY: CON_stop, CON_set_do_test 

    implicit none

    type(timetype), intent(in) :: StartTime, EndTime

    character(len=*), parameter :: NameSub='init_indices'
    !------------------------------------------------------------------------
    select case(NameIndexSource)
    case('file')
       call read_index_file(StartTime, EndTime, NameIndexFile)
    case default
       call CON_stop( &
            NameSub//'Unsupported geomagnetic index source: '//NameIndexSource)
    end select

  end subroutine init_indices

  !===========================================================================
  subroutine update_indices(timeNow)
    ! Interpolate Kp to current time.
    ! Use f10.7 according to current day.
    ! Input time format should be floating point used in ModTimeConvert.
    use ModRamGrids,     ONLY: nS
    use ModRamParams,    ONLY: FixedComposition, OfracN
    use ModRamVariables, ONLY: nRawKp, nRawF107, nRawAE, Kp, Kpmax12, Kpmax24, &
                               F107, AE, timeKp, timeF107, timeAE, rawKp, &
                               rawF107, rawAE, species
    use ModRamParams,    ONLY: DoUseEMIC
    use ModRamMain, ONLY: Real8_
    
    implicit none

    real(kind=Real8_), intent(in) :: timeNow
    
    integer :: iTime, i
    real(kind=Real8_) :: dTime, dateNow, BEXP, AHE0, AHE1, GEXP, Operc

    !------------------------------------------------------------------------
    ! NOTE: AS MORE SOURCES ARE ADDED, USE CASE STATEMENTS TO 
    ! CREATE DIFFERENT METHODS FOR OBTAINING THE INDICES AT timeNow.
    ! Find points in time about current time.  
    ! After this loop, iTime = position in timeKp after timeNow.
    iTime=1
    do while( (iTime < nRawKp) .and. (timeKp(iTime) < timeNow))
       iTime=iTime+1
    end do

    ! Interpolate Kp to current time.
    dTime = (timeNow - timeKp(iTime-1))/(timeKp(iTime) - timeKp(iTime-1))
    Kp = dTime*(rawKp(iTime) - rawKp(iTime-1)) + rawKp(iTime-1)
    Kpmax12 = max(Kp, maxval(rawKp(iTime-4:iTime-1)))
    Kpmax24 = max(Kpmax12, maxval(rawKp(iTime-8:iTime-5)))

    ! F10.7 index is not interpolated; merely use the value at the
    ! current day.
    dateNow=timeNow - mod(timeNow, 86400.0)
    do iTime=1, nRawF107
       if (abs(timeF107(iTime)-dateNow) .le. 1e-9) then
          f107 = rawF107(iTime)
          exit
       end if
    end do

    if(DoUseEMIC) then
       ! Find points in timeAE about current time.
       iTime = 1
       do while ( (iTime .lt. nRawAE) .and. (timeAE(iTime) .le. timeNow))
          iTime = iTime + 1
       end do
       ! Interpolate AE index to current time 
       dTime = (timeNow - timeAE(iTime-1))/(timeAE(iTime)-timeAE(iTime-1))
       AE = int(dTime*(rawAE(iTIme)-rawAE(iTime-1))+rawAE(iTime-1))
    else
       AE = 0
    end if


    ! If using the Young et al. composition model, recalculate the composition
    ! fractions based on new Kp and F10.7
    if (.not.FixedComposition) then
       BEXP=(4.5E-2)*EXP(0.17*KP+0.01*F107)
       AHE0=6.8E-3
       AHE1=0.084
       GEXP=0.011*EXP(0.24*KP+0.011*F107)
       GEXP=BEXP*(AHE0/GEXP+AHE1)
       do i = 1, nS
          select case(species(i)%s_name)
            case("Hydrogen")
              species(i)%s_comp = 4./(4.+BEXP+2.*GEXP)
            case("HeliumP1")
              species(i)%s_comp = 2.*GEXP/(4.+BEXP+2.*GEXP)
            case("OxygenP1")
              Operc = BEXP/(4.+BEXP+2.*GEXP)
              ! Subtract nitrogen percent from oxygen percent
              species(i)%s_comp = (1-OfracN)*Operc
            case("Nitrogen")
              Operc = BEXP/(4.+BEXP+2.*GEXP)
              ! Take a percent of the oxygen composition for nitrogen
              species(i)%s_comp = OfracN*Operc
            case default
              species(i)%s_comp = 1.0
          end select
       enddo
    endif
  end subroutine update_indices
  !===========================================================================

  subroutine test_update_indices(verbose)
    ! Tests:
    ! 1. Load of RamIndices.txt
    ! 2. Interpolation of Kp to specific times by update_indices
    use ModRamMain,      ONLY: make_time, Real8_, test_neq_abs, nTestPassed, nTestRun
    use ModTimeConvert,  ONLY: TimeType
    use ModRamVariables, ONLY: F107, Kp, Kpmax24, NameIndexSource
    use ModRamParams,    ONLY: DoUseEMIC, FixedComposition

    implicit none
    logical, intent(in) :: verbose
    logical :: failure
    type(TimeType) :: StartTime, EndTime, TestTime

    call make_time(2013, 3, 16, 0, 0, 0, 0.0, StartTime)
    call make_time(2013, 3, 17, 0, 0, 0, 0.0, TestTime)
    call make_time(2013, 3, 17, 5, 0, 0, 0.0, EndTime)

    NameIndexSource = 'file'
    DoUseEMIC = .FALSE.
    FixedComposition = .TRUE.
    ! Implicit test of index file intialization
    call init_indices(StartTime, EndTime)
    ! Explicit test of index interpolation
    call update_indices(TestTime%Time)
    call test_neq_abs(Kp, 1.65, 0.001, failure)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_update_indices(1): Interp. Kp. Expected: 1.65; Got: ", Kp
    if (.not.failure) nTestPassed = nTestPassed + 1
    ! Explicit test of index interpolation (#2)
    call make_time(2013, 3, 17, 4, 30, 0, 0.0, TestTime)
    call update_indices(TestTime%Time)
    call test_neq_abs(Kp, 2.3, 0.001, failure)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_update_indices(2): Interp. Kp. Expected: 2.3; Got: ", Kp
    if (.not.failure) nTestPassed = nTestPassed + 1
    ! Test max Kp in previous 24 (#3)
    call make_time(2013, 3, 17, 12, 00, 0, 0.0, TestTime)
    call update_indices(TestTime%Time)
    call test_neq_abs(Kpmax24, 6.7, 0.001, failure)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_update_indices(3): Test Kpmax24. Expected: 6.7; Got: ", Kpmax24
    if (.not.failure) nTestPassed = nTestPassed + 1
    ! Test of F10.7 retrieval (#4)
    call test_neq_abs(F107, 123.9, 0.001, failure)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_update_indices(4): Get F10.7 at noon. Expected: 123.9; Got: ", F107
    if (.not.failure) nTestPassed = nTestPassed + 1

  end subroutine test_update_indices
  !===========================================================================
end module ModRamIndices
