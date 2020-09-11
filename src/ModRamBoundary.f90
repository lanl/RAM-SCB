!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamBoundary
! Contains subroutines responsible for calculating the boundary flux for RAM

  use ModRamVariables, ONLY: FGEOS

  implicit none

  contains

!==============================================================================
subroutine get_boundary_flux
! Called by ModRamScbRun.f90

!!!! Module Variables
  use ModRamGrids,  ONLY: nS
  use ModRamTiming, ONLY: TimeRamElapsed
  use ModRamParams, ONLY: BoundaryFiles
  use ModRamIO,     ONLY: write_prefix

  implicit none
 
  integer :: iS

  ! Update boundary conditions and wave diffusion matrix
  do iS=1,nS
     IF (iS.EQ.1) THEN
        write(*,*) ''
        call write_prefix
        write(*,'(a,F10.2)') 'Calling GEOSB at time:    ', TimeRamElapsed
     ENDIF
     CALL GEOSB(iS)
  end do

end subroutine get_boundary_flux

!==============================================================================
subroutine get_geomlt_flux(NameParticleIn, fluxOut_II)
! Converts a LANL geomlt file into RAM boundary flux

!!!! Module Variables
  use ModRamMain,      ONLY: DP
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed, TimeRamStart
  use ModRamParams,    ONLY: electrons, ProtonFluxCap, ElectronFluxCap
  use ModRamGrids,     ONLY: NT, NE, NEL, NEL_prot, NBD, NTL
  use ModRamVariables, ONLY: EKEV, MLT, timeOffset, StringFileDate,  &
                             flux_SIII, fluxLast_SII, eGrid_SI, &
                             lGrid_SI, tGrid_SI

!!!! Module Subroutines and Functions
  use ModRamGSL, ONLY: GSL_Interpolation_1D
  use ModRamIO,  ONLY: read_geomlt_file


  implicit none

  integer :: GSLerr
  character(len=4), intent(in) :: NameParticleIn
  real(DP),intent(out):: fluxOut_II(nT, nE)

  character(len=8) :: StringDate
  integer :: iTime1, iTime2, iSpec, NEL_
  integer :: i, j, k
  integer :: lE, rE, pE
  real(DP) :: y
  real(DP), allocatable :: flux_II(:,:),logFlux_II(:,:),logELan(:),logERam(:)

  ! Usual debug variables.
  logical :: DoTest, DoTestMe, IsInitialized=.false.
  character(len=*), parameter :: NameSub='get_geomlt_flux'
  !------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,'(a, i2.2,":",i2.2)')NameSub//': Getting flux at t= ', &
       TimeRamNow%iHour, TimeRamNow%iMinute

  ! Initialize both species on first call.
  if(.not.IsInitialized)then
     if(DoTest)write(*,*)NameSub//': Initializing fluxes...'
     call read_geomlt_file('prot')
     if(electrons) call read_geomlt_file('elec')
     timeOffset=3600.0*TimeRamStart%iHour + 60.0*TimeRamStart%iMinute + &
          TimeRamStart%iSecond + TimeRamStart%FracSecond
     IsInitialized=.true.
  end if

  ! Check date of file against current date.
  write(StringDate, '(i4.4,i2.2,i2.2)') TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay
  if(DoTest)write(*,*)'Nowdate, Filedate = ', Stringdate, '  ', StringFileDate
  ! If the day has changed, read in the new geomlt file
  if(StringDate.ne.StringFileDate)then
     call read_geomlt_file('prot')
     if(electrons) call read_geomlt_file('elec')
  end if

  ! Possible different sizes for two boundary input files
  ! These sizes are configurable in the PARAM file
  if (NameParticleIn .eq. 'prot') then
     NEL_  = NEL_prot
     iSpec = 1
  elseif (NameParticleIn.eq.'elec') then
     NEL_  = NEL
     iSpec = 2
  endif

  if(allocated(flux_II)) then
     deallocate(flux_II,logFlux_II,logELan,logERam)
  endif

  ! May need to add energy points on either end of the boundary files
  ! depending on if our desired energy range for the simulation lies
  ! outside of the files energy range.
  rE = 0
  lE = 0
  if (eGrid_SI(iSpec,1).gt.EkeV(1))     lE = 1
  if (eGrid_SI(iSpec,NEL_).lt.EkeV(nE)) rE = 1
  pE = rE + lE
  allocate(flux_II(0:NTL,NEL_+pE), logFlux_II(nT,NEL_+pE), logELan(NEL_+pE), logERam(nE))
  flux_II = 0.0; logFlux_II = 0.0; logELan = 0.0; logERam = 0.0

  ! Get closest times from geomlt file
  iTime1 = 0
  iTime2 = 0
  do i=1,NBD
     if (abs(TimeRamElapsed-tGrid_SI(iSpec,i)).le.1e-9) then
        iTime2 = i
        iTime1 = i
        exit  
     elseif (TimeRamElapsed.lt.tGrid_SI(iSpec,i)) then
        iTime2 = i
        iTime1 = max(1,i-1)
        exit
     endif
  enddo
  if ((iTime1.eq.0).and.(iTime2.eq.0)) then
     if (TimeRamElapsed.lt.tGrid_SI(iSpec,1)) then
        iTime1 = 1
        iTime2 = 1
     else
        iTime1 = NBD
        iTime2 = NBD
     endif
  endif

  ! If at end of full day of flux, set fluxLast for interpolation.
  if (iTime1==NBD) fluxLast_SII(iSpec,:,:)=flux_SIII(iSpec,NBD,:,:)

  ! If needed, interpolate between file time steps
  if (iTime1.eq.iTime2) then
    flux_II(1:24,1+lE:NEL_+le) = flux_SIII(iSpec,iTime1,:,1:NEL_)
  else
     do j=1,NTL-1
        do k=1,NEL_
           CALL GSL_Interpolation_1D(tGrid_SI(iSpec,:),flux_SIII(iSpec,:,j,k), &
                                     TimeRamElapsed,flux_II(j,k+lE),GSLerr)
        enddo
     enddo
  endif

  ! If needed, extrapolate to lower and higher energies 
  ! (for now just set to 10^8 and 0.1 for testing)
  if (lE.eq.1) then
     do j=1,NTL-1
        flux_II(j,1) = 10.**8
     enddo
  endif

  if (rE.eq.1) then
     do j=1,NTL-1
        flux_II(j,NEL_+pE) = 0.1
     enddo
  endif

  ! Create array that is periodic in L.
  flux_II(0,:)  = flux_II(24,:)
  flux_II(25,:) = flux_II(1,:)

  if(DoTest)then
     write(*,*)NameSub//'Raw flux at E-min='
     write(*,'(3(ES12.5,1x))')flux_II(:,1)
  end if

  ! Interpolate in LT-space, get log for E-interpolation...
  do j=1,nT
     do k=1,NEL_+pE
        call GSL_Interpolation_1D(lGrid_SI(iSpec,:),flux_II(:,k), MLT(j), y, GSLerr)
        logFlux_II(j,k)=log10(y)
     end do
  end do

  if(DoTest)then
     write(*,*)NameSub//'Log10-LT-interpolated flux at E-min='
     write(*,'(3(ES12.5,1x))')LogFlux_II(:,1)
  end if

  ! Place energy grids into log space.
  logELan(1+lE:NEL_+le) = log10(eGrid_SI(iSpec,1:NEL_))
  if (lE.eq.1) logELan(1)       = log10(0.1000)
  if (rE.eq.1) logELan(NEL_+pE) = log10(1000.00)
  logERam = log10(Ekev)

  ! Interpolate/Extrapolate in energy space; return to normal units.
  ! Interpolation in Log space isn't really needed with the GSL interpolation
  ! routines but is left over from the previous interpolation routines. It
  ! shouldn't mess anything up so for now we will leave it this way. -ME
  do j = 1, nT
     do k = 1, nE
        CALL GSL_Interpolation_1D(logELan, logFlux_II(j,:), logERam(k), y, GSLerr)
        fluxOut_II(j,k)=10.**y
        if (NameParticleIn.eq.'elec') then
           if (fluxOut_II(j,k) > ElectronFluxCap) fluxOut_II(j,k) = ElectronFluxCap
        end if
        if (NameParticleIn.eq.'prot') then
           if (fluxOut_II(j,k) > ProtonFluxCap) fluxOut_II(j,k) = ProtonFluxCap
        end if
     end do
  end do

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

!==================================================================================================
! *************************************************************************
!                                GEOSB
!                       Boundary conditions set up
!**************************************************************************
  SUBROUTINE GEOSB(S)

!!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: IsComponent, boundary, BoundaryFiles, WriteBoundary, &
                               FixedComposition, verbose
    use ModRamGrids,     ONLY: NTL, NEL, NT, NE, NR
    use ModRamTiming,    ONLY: TimeRamElapsed
    use ModRamVariables, ONLY: FFACTOR, UPA, Kp, F107, species
    use ModRamCouple,    ONLY: FluxBats_IIS

!!!! Module Subroutines and Functions
    use ModRamIO, ONLY: write_dsbnd

    implicit none
    integer, intent(in) :: S

    integer  :: ij, ik, l, u
    real(DP) :: bexp, ahe0, ahe1, gexp, T
    real(DP), ALLOCATABLE :: RELAN(:,:), FLAN(:,:), FluxLanl(:,:)

    ALLOCATE(RELAN(NTL,NEL),FLAN(NT,NEL),FluxLanl(nT,nE))
    RELAN = 0.0; FLAN = 0.0; FluxLanl = 0.0
    T = TimeRamElapsed

    ! Zero out FGEOS at this species.
    FGEOS(S,:,:,:)=0.

    ! Read LANL flux (1/cm2/s/sr/keV) assumed isotropic
    IF (boundary.EQ.'LANL') THEN
      ! LANL interpolated flux files.
      select case(species(S)%s_name)
      case ("Electron")
        call get_geomlt_flux('elec', FluxLanl)
      case ("Hydrogen", "OxygenP1", "HeliumP1", "Nitrogen")
        ! We assume a fraction of the oxygen is actually nitrogen
        ! The specific percentage is configurable in the PARAM file
        ! By default the nitrogen fraction is assumed to be zero
        call get_geomlt_flux('prot', FluxLanl)
      case default
        FluxLanl = 0._dp
      end select
      FluxLanl(1,:) = FluxLanl(nT,:)
      ! Adjust flux for composition
      ! composition is either fixed at the beginning of a run or calculated in
      ! ModRamIndices based on the Young et al. composition model
      FluxLanl=FluxLanl*species(S)%s_comp
      do ik=1,NE
        do ij=1,nT
          u = int(upa(nR)-1,kind=4)
          do L=2,u
            FGEOS(s,iJ,iK,L)=FluxLanl(iJ,iK) * FFACTOR(S,NR,IK,L)
          end do
        end do
      end do

    ELSEIF (boundary .EQ. 'SWMF') THEN
      !! Read SWMF flux (1/cm2/s/sr/keV) assumed isotropic
      do iK=1, nE
        do iJ=1,nT
          u = int(UPA(NR)-1,kind=4)
          do L=2,u
            FGEOS(s,iJ,iK,L)=FluxBats_IIS(iK,iJ,s)*FFACTOR(S,NR,IK,L)
          end do
        end do
      end do
    END IF

    if (verbose) then
       write(*,'(1x,a,2E10.2)') 'Max, Min '//species(S)%s_code//' boundary flux = ', &
             maxval(FGEOS(s,:,:,:)), minval(FGEOS(s,:,:,:), MASK=FGEOS(s,:,:,:).GE.0)
    endif

    ! Write interpolated fluxes to file.
    if (WriteBoundary) call write_dsbnd(S)

    DEALLOCATE(RELAN,FLAN,FluxLanl)
    RETURN
  END SUBROUTINE GEOSB

END MODULE ModRamBoundary
