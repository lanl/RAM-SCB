!============================================================================
module ModRamIndices
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================
  
  ! A set of tools for handling the indices used in RAM-SCB.  Combined with
  ! the PARAM interface, allows the user to select source of indices, gather
  ! and interpolate the values to the current simulation time.

  use ModRamMain, ONLY: Real8_

  implicit none
  save

  ! Source of indices
  character(len=4) :: NameIndexSource = 'file'
  character(len=200) :: NameIndexFile = 'RamIndices.txt'
  character(len=200) :: NameOmniFile  = 'omni.txt'

  integer :: nRawKp, nRawF107 ! Number of entries read from file
  real(kind=Real8_), allocatable :: timeKp(:),timeF107(:),rawKp(:),rawF107(:)

contains
  !===========================================================================
  subroutine read_index_file(StartTime, EndTime, NameFile)

    use ModTimeConvert
    use ModIoUnit, ONLY: UNITTMP_

    implicit none

    type(timetype),   intent(in) :: StartTime, EndTime
    character(len=*), intent(in) :: NameFile

    integer :: is, fmday
    character(len=100) :: header
    real(kind=Real8_) :: tmpF107

    integer :: i, j, nline, iError, cmday, iYY, iMM, iDD
    character(len=100) :: StringLine, StringFmt
    real(kind=Real8_) :: tmpKp(8)

    integer, parameter :: monthday(12) = &
         (/0,30,59,90,120,151,181,212,243,273,304,334/)
    integer, parameter :: leapday(12) = &
         (/0,30,60,91,121,152,182,213,244,274,305,335/)
    integer, parameter :: kptime(8) = &
         (/1, 4, 7, 10, 13, 16, 19, 22/)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='read_index_file'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest) write(*,*)'Loading indices from ', trim(NameFile)

    ! Open file and find the starting date of file
    open(unit=UNITTMP_, FILE=NameFile, STATUS='OLD', IOSTAT=iError)
    if(iError /=0 ) call CON_stop &
         (NameSub//' Error opening file '//NameFile)
    read(UNITTMP_,*) header
    i = 1
    read(UNITTMP_, '(i4, i2, i2, 8(1x,f3.1), 1x, f5.1)') &
                     iYY, iMM, iDD, tmpKp, tmpF107

    if ((StartTime%iYear.eq.iYY).and. &
        (StartTime%iMonth.eq.iMM).and. &
        (StartTime%iDay.eq.iDD)) then
       is = 2
       if (DoTest) then
          write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
                   StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
                   StartTime%iHour, StartTime%iMinute, StartTime%iSecond
          write(*,*) NameSub//': Requested start date matches RamIndices start date'
       end if
    else
       is = 1
       ! Find starting point in kyoto file by determining the number
       ! of lines up to starting position.
       if( mod(StartTime%iYear,4)==0 )then
          cmday = leapday(StartTime%iMonth)
       else
          cmday = monthday(StartTime%iMonth)
       end if
       if( mod(iYY,4)==0 )then
          fmday = leapday(iMM)
       else
          fmday = monthday(iMM)
       end if
       nline = 365*(StartTime%iYear - iYY) + (StartTime%iyear - iYY + 1)/4 &
               + (cmday + StartTime%iDay) - (fmday + iDD)

       if (nline.lt.1) then
          call CON_stop( &
            NameSub//': Requested start date is before RamIndices start date')
       elseif (nline.eq.1) then
          if(DoTest)then
             write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
                   StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
                   StartTime%iHour, StartTime%iMinute, StartTime%iSecond
             write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'End Time   = ', &
                   EndTime%iYear, EndTime%iMonth, EndTime%iDay, &
                   EndTime%iHour, EndTime%iMinute, EndTime%iSecond
             write(*,*) NameSub//': Line before requested date is start of RamIndices'
          end if
       elseif (nline.eq.2) then
          read(UNITTMP_, *, IOSTAT=iError) StringLine
          if (iError.lt.0) then
             call CON_stop( &
               NameSub//': End date outside of range of RamIndices file')
          elseif (iError.gt.0) then
             call CON_stop( &
               NameSub//': Error in formating of RamIndices file')
          end if
          if(DoTest)then
             write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'Start Time = ', &
                      StartTime%iYear, StartTime%iMonth, StartTime%iDay, &
                      StartTime%iHour, StartTime%iMinute, StartTime%iSecond
             write(*,'(a, i4, 2("-",i2.2),1x, i2.2,2(":",i2.2))') 'End Time   = ', &
                      EndTime%iYear, EndTime%iMonth, EndTime%iDay, &
                      EndTime%iHour, EndTime%iMinute, EndTime%iSecond
             write(*,*) NameSub//': Line before requested date:'
             write(*,*) StringLine
          end if
       else
          ! Fast forward to current position.
          write(StringFmt, "('(',i10,'(/)a)')") nLine-2
          read(UNITTMP_, trim(StringFmt),IOSTAT=iError) StringLine
          if (iError.lt.0) then
             call CON_stop( &
               NameSub//': End date outside of range of RamIndices file')
          elseif (iError.gt.0) then
             call CON_stop( &
               NameSub//': Error in formating of RamIndices file')
          end if
          if(DoTest)then
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
    end if

    ! Allocate arrays based on number of days required.
    nRawF107 = ceiling((EndTime%Time-StartTime%Time)/86400.0)+1
    nRawKp   = 8*nRawF107
    allocate(timeKp(nRawKp))
    allocate(rawKp(nRawKp))
    allocate(rawF107(nline))
    allocate(timeF107(nline))

    if(DoTest)then
       write(*,*) NameSub//': Number of lines to read = ', nRawF107
       write(*,*) NameSub//': Number of KP vals to read=', nRawKp
    end if

    if (is.eq.2) then
       rawF107(1) = tmpF107
       call time_int_to_real((/iYY, iMM, iDD, 0,0,0,0/), timeF107(i))
       do j=1, 8
          rawKp(8*(i-1)+j) = tmpKp(j)
          call time_int_to_real((/iYY, iMM, iDD, kptime(j),30,0,0/),&
               timeKp(8*(i-1)+j))
       end do
    end if

    ! Read and store all data for entire interval.
    do i=is, nRawF107
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
    
  end subroutine read_index_file

  !===========================================================================
  subroutine init_indices(StartTime, EndTime)
    ! Based on source of indices, prepare indices for this simulation.
    use ModTimeConvert

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
  subroutine get_indices(timeNow, kpNow, f10Now)
    ! Interpolate Kp to current time.
    ! Use f10.7 according to current day.
    ! Input time format should be floating point used in ModTimeConvert.
    use ModRamMain, ONLY: Real8_
    
    implicit none

    real(kind=Real8_), intent(in) :: timeNow
    real(kind=Real8_), intent(out):: kpNow, f10Now
    
    integer :: iTime
    real(kind=Real8_) :: dTime, dateNow, nSecInDay=86400.0

    character(len=*), parameter :: NameSub='get_indices'
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
    dTime = (timeNow - timeKp(iTime-1)) / &
         (timeKp(iTime) - timeKp(iTime-1))
    kpNow = dTime*(rawKp(iTime) - rawKp(iTime-1)) + rawKp(iTime-1)

    ! F10.7 index is not interpolated; merely use the value at the
    ! current day.
    dateNow=timeNow - mod(timeNow, nSecInDay)
    do iTime=1, nRawF107
       if (timeF107(iTime) .eq. dateNow) then
          f10Now = rawF107(iTime)
          exit
       end if
    end do
    
  end subroutine get_indices
  !===========================================================================
  subroutine get_ramdst(dstOut)
    ! Use simple DPS relationship to calculate Dst from RAM domain.
    ! Sums over all species, corrects for internal currents (factor of 1.3).
    use ModRamMain, ONLY: nR, nE, nPa, nT, f2, factor, &
         upa, we, wmu, ekev, Real8_

    implicit none

    real(kind=Real8_), intent(out) :: dstOut

    ! Factor2 includes conversions and factor of 1.3.
    real(kind=Real8_)           :: sumEnergy=0.0, factor2=-5.174E-30
    integer                     :: i, j, k, l, s
    character(len=*), parameter :: NameSub = 'get_ramdst'
    !------------------------------------------------------------------------
    sumEnergy = 0.0
    dstOut = 0.0
    ! Sum energy over whole domain and all species.
    do s=1,4; do i=2,nR; do k=2,nE; do l=1,nPa
       if(l.ge.uPa(i))cycle
       do j=1, nT-1
          sumEnergy=sumEnergy+f2(s,i,j,k,l)*wE(k)*wMu(L)*eKeV(k)
       end do
    end do; end do; end do; end do

    dstout = factor2 * factor * sumEnergy

  end subroutine get_ramdst
  !===========================================================================
end module ModRamIndices
