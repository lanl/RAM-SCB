!============================================================================
module ModRamGeomlt
!    A set of tools for reading, interpolating, and setting plasma boundary
!    conditions from "geomlt" flux files created from LANL geosynchronous 
!    spacecraft.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModRamMain, ONLY: Real8_

  implicit none
  save

  logical :: IsInitialized = .false.

  ! Start time offset in seconds.
  real(kind=Real8_) :: timeOffset

  ! Date of file:
  character(len=8) :: StringFileDate

  ! Values that are read and stored from file:
  real(kind=Real8_), allocatable :: flux_SIII(:,:,:,:), fluxLast_SII(:,:,:), eGrid_SI(:,:), avgSats_SI(:,:)
  real(kind=Real8_) :: lGrid_SI(2,0:25)

contains
  !===========================================================================
subroutine read_geomlt_file(NameParticle)
  
  use ModIoUnit,  ONLY: UnitTMP_
  use ModRamMain, ONLY: TimeRamNow, PathRamIn, nE, Ekev, Dt_bc, NEL, NEL_prot, NBD, boundary

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
subroutine get_geomlt_flux(NameParticleIn, fluxOut_II)
  
  use ModRamMain, ONLY: TimeRamNow, TimeRamElapsed, electrons, &
       nT, nE, Ekev, MLT, Dt_bc, NEL, NEL_prot, NBD, boundary

  implicit none

  character(len=4), intent(in) :: NameParticleIn
  real(kind=Real8_),intent(out):: fluxOut_II(nT, nE)

  character(len=8) :: StringDate
  integer :: iError, iTime, iSpec=1, j, k, NEL_
  real(kind=Real8_) :: y, my_var
  real(kind=Real8_), allocatable :: flux_II(:,:), logFlux_II(:,:), logELan(:), logERam(:)

  ! Usual debug variables.
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='get_geomlt_flux'
  !------------------------------------------------------------------------

  ! Possible different sizes for two boundary input files
  if((NameParticleIn .eq. 'prot').and.(boundary .eq. 'PTM')) then
     NEL_ = NEL_prot
  else
     NEL_ = NEL
  endif

  if(allocated(flux_II)) then
     deallocate(flux_II,logFlux_II,logELan,logERam)
  endif

  allocate(flux_II(0:25,NEL_),logFlux_II(nT,NEL_),logELan(NEL_),logERam(nE))

  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,'(a, i2.2,":",i2.2)')NameSub//': Getting flux at t= ', &
       TimeRamNow%iHour, TimeRamNow%iMinute

  ! Initialize both species on first call.
  if(.not.IsInitialized)then
     if(DoTest)write(*,*)NameSub//': Initializing fluxes...'
     call read_geomlt_file('prot')
     if(electrons) call read_geomlt_file('elec')
     timeOffset=3600.0*TimeRamNow%iHour + 60.0*TimeRamNow%iMinute + &
          TimeRamNow%iSecond + TimeRamNow%FracSecond
     IsInitialized=.true.
  end if
  
  ! Check date of file against current date.
  write(StringDate, '(i4.4,i2.2,i2.2)') &
       TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay
  if(DoTest)write(*,*)'Nowdate, Filedate = ', Stringdate, '  ', StringFileDate
  if(StringDate.ne.StringFileDate)then
     call read_geomlt_file('prot')
     if(electrons) call read_geomlt_file('elec')
  end if

  ! Get the position in the array by determining 
  ! how many 5-min intervals have passed.
  ! Chris Jeffery** Change this from 300 to couple with Dt_bc
  !  iTime=nint(mod(TimeRamElapsed+timeOffset,86400.0)/300.0)+1
  
  ! Also, don't use timeOffset if I'm running PTM
  if(boundary .eq. 'PTM') then
     iTime=nint(mod(TimeRamElapsed,86400.0)/Dt_bc)+1
  else
     iTime=nint(mod(TimeRamElapsed+timeOffset,86400.0)/Dt_bc)+1
  endif

  if(DoTest) write(*,'(a,i3.3)')NameSub//': Index iTime = ', iTime

  ! Set species index.
  iSpec=1  ! Keep this or iSpec will be saved for next call.
  if(NameParticleIn.eq.'elec') iSpec=2

  ! If at end of full day of flux, set fluxLast.
  if (iTime==NBD) fluxLast_SII(iSpec,:,:)=flux_SIII(iSpec,iTime,:,:)

  ! Create array that is periodic in L.
  flux_II(1:24,:) = flux_SIII(iSpec,iTime, :,1:NEL_)
  flux_II(0,:)    = flux_SIII(iSpec,iTime,24,1:NEL_)
  flux_II(25,:)   = flux_SIII(iSpec,iTime,1,1:NEL_)

  if(DoTest)then
     write(*,*)NameSub//'Raw flux at E-min='
     write(*,'(3(ES12.5,1x))')flux_II(:,1)
  end if

  ! Interpolate in LT-space, get log for E-interpolation...
  do j=1,nT; do k=1,NEL_
     call lintp(lGrid_SI(iSpec,:), flux_II(:,k), 26, MLT(j), y, iError)
     logFlux_II(j,k)=log10(y)
  end do; end do

  if(DoTest)then
     write(*,*)NameSub//'Log10-LT-interpolated flux at E-min='
     write(*,'(3(ES12.5,1x))')LogFlux_II(:,1)
  end if

  ! Place energy grids into log space.
  logELan=log10(eGrid_SI(iSpec,1:NEL_))
  logERam=log10(Ekev)

  ! Interpolate in energy space; return to normal units.
  do j=1,nT; do k=1,nE
     call lintp(logELan, logFlux_II(j,:), NEL_, logERam(k), y, iError)
     if (y.gt.8)  then 
       y=8
       write(*,*) ' in ModRamGeomlt: limit flux to 1e8'
     end if
     fluxOut_II(j,k)=10.**y
  end do; end do
 
  do j=2,nT; do k=1,nE
     my_var = fluxOut_II(j,k)
    if (my_var /= my_var) then     ! doesn't work
      write(*,*) ' FluxLanl is NaN'
      fluxOut_II(j,k)=fluxOut_II(j,k-1)
      write(*,*) ' new flux =', fluxOut_II(j,k)
    end if
    if (my_var == 10.**8) then
 !     fluxOut_II(j,k)=fluxOut_II(j,k-1)
      fluxOut_II(j,k)=fluxOut_II(j-1,k)
 !     write(*,*) ' new flux =', fluxOut_II(j,k)
    endif
 end do; end do

 do k=1,nE
    fluxOut_II(1,k)=fluxOut_II(nT,k)
 end do


  ! With fluxOut filled, we are finished.
  if(DoTest) then
     write(*,*)NameSub//': Finished.  Flux at min-E about LT='
     write(*,'(3(ES12.5,1x))')fluxOut_II(:,1)
     write(*,*)NameSub//': Check output_ram/Dsbnd/ files for more info.'
  end if

end subroutine get_geomlt_flux
  !===========================================================================

end module ModRamGeomlt
!============================================================================
