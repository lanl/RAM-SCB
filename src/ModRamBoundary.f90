MODULE ModRamBoundary
! Contains subroutines responsible for calculating the boundary flux for RAM

  use ModRamVariables, ONLY: FGEOS

  implicit none
  save

  contains

!==============================================================================
subroutine get_boundary_flux

  use ModRamMain,   ONLY: S
  use ModRamTiming, ONLY: TimeRamElapsed
  use ModRamParams, ONLY: boundary

  implicit none
 
  integer :: iS
  character(len=4) :: ST7no
  ! Update boundary conditions and wave diffusion matrix
  do iS=1,4
     S = iS
     IF (S.EQ.1) THEN
        if (boundary .EQ. 'LANL') print*, &
           'RAM: calling GEOSB at time (hr) = ', TimeRamElapsed/3600.
     ENDIF
     CALL GEOSB
  end do

end subroutine

!==============================================================================
subroutine get_geomlt_flux(NameParticleIn, fluxOut_II)

  use ModRamMain,      ONLY: Real8_, S
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed, TimeRamStart, Dt_bc
  use ModRamParams,    ONLY: electrons, boundary
  use ModRamGrids,     ONLY: NT, NE, NEL, NEL_prot, NBD
  use ModRamVariables, ONLY: EKEV, MLT, FGEOS, timeOffset, StringFileDate, &
                             flux_SIII, fluxLast_SII, eGrid_SI, avgSats_SI, lGrid_SI

  use ModRamIO, ONLY: read_geomlt_file

  implicit none

  character(len=4), intent(in) :: NameParticleIn
  real(kind=Real8_),intent(out):: fluxOut_II(nT, nE)

  character(len=8) :: StringDate
  integer :: iError, iTime, iSpec=1, j, k, NEL_
  real(kind=Real8_) :: y, my_var
  real(kind=Real8_), allocatable :: flux_II(:,:),logFlux_II(:,:),logELan(:),logERam(:)

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
!  if(.not.IsInitialized)then
  if (S.eq.1) then
     if(DoTest)write(*,*)NameSub//': Initializing fluxes...'
     call read_geomlt_file('prot')
     if(electrons) call read_geomlt_file('elec')
     timeOffset=3600.0*TimeRamStart%iHour + 60.0*TimeRamStart%iMinute + &
          TimeRamStart%iSecond + TimeRamStart%FracSecond
!     IsInitialized=.true.
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
       write(*,*) ' in ModRamBoundary: limit flux to 1e8'
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

!==============================================================================
! *************************************************************************
!                                GEOSB
!                       Boundary conditions set up
!**************************************************************************
  SUBROUTINE GEOSB

    use ModRamMain,      ONLY: Real8_, S, PathSwmfOut, PathRamOut
    use ModRamParams,    ONLY: DoAnisoPressureGMCoupling, IsComponent, boundary, &
                               DoMultiBcsFile
    use ModRamGrids,     ONLY: NTL, NEL, NT, NE, NR
    use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamNow, TimeRamStart
    use ModRamVariables, ONLY: FFACTOR, UPA, EKEV, Kp, F107
    use ModRamCouple,    ONLY: FluxBats_IIS, generate_flux, TypeMHD, MhdDensPres_VII, &
                               FluxBats_anis

    use ModRamIO, ONLY: write_dsbnd

    use ModIoUnit, ONLY: UNITTMP_
    USE ModConst,  ONLY: cProtonMass

    implicit none
    save
    integer             :: ij, ik, j, jw, k, l, nLines
    real(kind=Real8_)   :: bexp, ahe0, ahe1, gexp, doy, azir, &
                           fracComposition, T
    character(len=200)  :: NameFluxFile
    character(len=90)   :: HEADER
    real(kind=Real8_)   :: RELAN(NTL,NEL),FLAN(NT,NEL),FluxLanl(nT,nE)

    T = TimeRamElapsed

    ! Create Young et al. composition ratios.
    BEXP=(4.5E-2)*EXP(0.17*KP+0.01*F107)
    AHE0=6.8E-3
    AHE1=0.084
    GEXP=0.011*EXP(0.24*KP+0.011*F107)
    GEXP=BEXP*(AHE0/GEXP+AHE1)
    select case(s)
      case(2)! Fraction H+
        fracComposition=4./(4.+BEXP+2.*GEXP)
      case(3)! Fraction He+
        fracComposition=2.*GEXP/(4.+BEXP+2.*GEXP)
      case(4)! Fraction O+
        fracComposition=BEXP/(4.+BEXP+2.*GEXP)
      case default ! Electrons.
        fracComposition=1.0
    end select

    ! Zero out FGEOS at this species.
    FGEOS(S,:,:,:)=0.

    ! Read LANL flux (1/cm2/s/sr/keV) assumed isotropic
    IF ((boundary .eq. 'LANL') .OR. (boundary .eq. 'PTM')) THEN
      ! LANL interpolated flux files.
      if (s==1) then
        call get_geomlt_flux('elec', FluxLanl)
      else
        call get_geomlt_flux('prot', FluxLanl)
      end if
      ! fix corrupted bc, NAN check doesn't work? limit ion flux:
      if (S.gt.1) then
        do ij=2,nT
          do ik=1,nE
            if (FluxLanl(iJ,iK).gt.1e7) FluxLanl(iJ,iK)=FluxLanl(iJ-1,iK)
            if (FluxLanl(iJ,iK).gt.1e7) then
              FluxLanl(iJ,iK)=1e7
              write(*,*) ' in GEOSB: limit ion flux to 1e7'
            endif
          end do
        end do
        do ik=1,nE
          FluxLanl(1,iK)=FluxLanl(nT,iK)
        end do
      end if
      ! corrupted bc, end of fix
      ! Adjust flux for composition
      FluxLanl=FluxLanl*fracComposition
      do ik=1,NE
        do ij=1,nT
          do L=2,upa(NR)-1
            FGEOS(s,iJ,iK,L)=FluxLanl(iJ,iK) * FFACTOR(S,NR,IK,L)
          end do
        end do
      end do
    ELSEIF (boundary .EQ. 'SWMF') THEN
      ! Read SWMF flux (1/cm2/s/sr/keV) assumed isotropic
      if (IsComponent) then ! If SWMF component, get flux from Bats.
        write(*,*) 'RAM_SCB: Getting flux from BATS-R-US'
        do iK=1, nE
          do iJ=1,nT
            do L=2,UPA(NR)-1
              if (.not.DoAnisoPressureGMCoupling) then
                FGEOS(s,iJ,iK,L)=FluxBats_IIS(iK,iJ,s)*FFACTOR(S,NR,IK,L)
              else
                ! anisotropy coupling is pitch angle dependent.
                FGEOS(s,iJ,iK,L)=FluxBats_anis(iK,L,iJ,s)*FFACTOR(S,NR,iK,L)
              endif
            end do
          end do
        end do
      else ! Open file if running stand alone.
!        if (DoMultiBcsFile) then
!          write(NameFluxFile, "(a,'/',a,'_',i1,'.swf')") trim(PathSwmfOut), ST7, S
!          nlines = 25
!        else
!          write(NameFluxFile, "(a,'/',a,'.swf')") trim(PathSwmfOut), ST7
!          nlines=49
!        endif
!        write(*,'(2a)')'RAM_SCB: Loading flux from ',trim(NameFluxFile)
!        OPEN(UNIT=UNITTMP_,FILE=trim(NameFluxFile),STATUS='OLD')
!        READ(UNITTMP_,*) HEADER
!        DO IJ=1,nLines
!          READ(UNITTMP_,*) DOY,AZIR,(RELAN(IJ,IK),IK=1,NE)
!        ENDDO
!        CLOSE(UNITTMP_)
!        DO IK=1,NE
!          DO IJ=1,NT
!            IF (NT.EQ.49) JW=IJ
!            if ((NT.EQ.25) .and. DoMultiBcsFile) then
!              JW=IJ
!            else
!              JW=2*IJ-1
!            end if
!            FLAN(IJ,IK) = RELAN(JW,IK)
            ! same mass density at geo as without composition
!            IF ((S.EQ.2) .and. .not.DoMultiBcsFile) FLAN(IJ,IK) = FLAN(IJ,IK)/(1+16*BEXP+4*GEXP)
!            IF ((S.EQ.3) .and. .not.DoMultiBcsFile) FLAN(IJ,IK) = FLAN(IJ,IK)*GEXP/(2+32*BEXP+8*GEXP)
!            IF ((S.EQ.4) .and. .not.DoMultiBcsFile) FLAN(IJ,IK) = FLAN(IJ,IK)*BEXP/(4+64*BEXP+16*GEXP)
!            DO L=2,UPA(NR)-1
!              FGEOS(S,IJ,IK,L)=FLAN(IJ,IK)*FFACTOR(S,NR,IK,L)
!            ENDDO
!          ENDDO
!        ENDDO
      end if
    END IF

    ! Write interpolated dfluxes to file.
    call write_dsbnd

    RETURN
  END

END MODULE ModRamBoundary
