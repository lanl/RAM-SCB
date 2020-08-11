!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamEField
! Contains subroutines related to getting the electric field for RAM

  implicit none

  contains

!==============================================================================
subroutine get_electric_field

  use ModRamMain,      ONLY: Real8_
  use ModRamTiming,    ONLY: TimeRamNow, DtEfi, TimeRamElapsed
  use ModRamConst,     ONLY: RE
  use ModRamParams,    ONLY: electric, IsComponent
  use ModRamGrids,     ONLY: NR, NT
  use ModRamVariables, ONLY: KP, PHI, LZ, PHIOFS, VT, VTOL, VTN, TOLV
  use ModRamIO,        ONLY: write_prefix

  use ModTimeConvert, ONLY: TimeType, time_real_to_int


  implicit none

  real(kind=Real8_) :: AVS
  integer :: I, J
  type(TimeType) :: TimeNext
  character(len=214) :: NameEfile

  ! Set "new" values of indices and E-field to "old" values. 
  VTOL = VTN
  write(*,*) ''
  call write_prefix
  write(*,'(a,F10.2)')'Updating E field at time: ', TimeRamElapsed

  ! Open this file and save contents.
  if(electric .ne. 'VOLS') call ram_get_electric(vtn)

  ! Set time of "old" file for interpolation:
  TOLV=TimeRamElapsed

  ! No interpolation for IESC/RSCE/WESC/W5SC  cases.
  if (electric=='IESC' .or. electric == 'RSCE' &
    .OR. electric=='WESC' .or. electric=='W5SC') then
     VTOL = VTN
  endif

  DO I=1,NR+1
     DO J=1,NT
        if (electric .ne. 'VOLS') then
           VT(I,J)=VTOL(I,J)+(VTN(I,J)-VTOL(I,J))*(TimeRamElapsed-TOLV)/DtEfi
        else
           AVS=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE  ! Voll-Stern parameter
           VT(I,J)=AVS*(LZ(I)*RE)**2*SIN(PHI(J)-PHIOFS) ! [V]
        endif
     ENDDO
  ENDDO

end subroutine get_electric_field

!============================================================================
subroutine ram_get_electric(EOut_II)
  ! Open E-field file "NextFile", collect E-field and indices, update
  ! "NextFile" and return all of the above to the caller.

  use ModRamMain,      ONLY: DP, PathRamOut
  Use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamNow
  use ModRamParams,    ONLY: electric, IsComponent, verbose, WritePotential
  use ModRamGrids,     ONLY: NR, NT, RadiusMax, RadiusMin
  use ModRamVariables, ONLY: PHIOFS, Kp, F107, VT
  use ModRamFunctions, ONLY: RamFileName
  use ModRamGSL,       ONLY: GSL_Interpolation_2D
  use ModScbMain,      ONLY: prefixOut
  use ModScbGrids,     ONLY: nzeta
  use ModScbVariables, ONLY: PhiIono, radRaw, azimRaw, x, y, nThetaEquator

  use ModRamFunctions, ONLY: RamFileName

  use ModIOUnit,      ONLY: UNITTMP_

  use nrtype, ONLY: pi_d

  implicit none

  ! Arguments
  real(DP), intent(inout) :: EOut_II(:,:)
  
  integer :: iError, i, j, k, jw, GSLerr
  real(DP) :: day, tth, ap, rsun, RRL, PH, wep(49), radOut, KpOut, F107Out
  REAL(DP), ALLOCATABLE :: Epot_Cart(:,:), xo(:,:), yo(:,:)
  character(len=200) :: StringHeader, NameFileOut
  character(len=*), parameter :: NameSub = 'ram_get_electric'
  !--------------------------------------------------------------------------
  ALLOCATE(Epot_Cart(nR+1,nT), xo(nR+1,nT), yo(nR+1,nT))
  Epot_Cart = 0.0; xo = 0.0; yo = 0.0

  radOut = RadiusMax+0.25

  EOut_II  = 0.0
  Epot_Cart = 0.0
  DO j = 0,nR
     radRaw(j) = RadiusMin + (radOut-RadiusMin) * REAL(j,DP)/REAL(nR,DP)
  END DO
  DO k = 1,nT
     azimRaw(k) = 24.0 * REAL(k-1,DP)/REAL(nT-1,DP)
  END DO

  ! Get potential on ionospheric grid
  CALL ionospheric_potential

  ! Interpolate from equatorially mapped ionospheric grid to radial grid
  DO i = 1,nR+1
     DO j = 1,nT
        xo(i,j) = radRaw(i-1) * COS(azimRaw(j)*2._dp*pi_d/24._dp - pi_d)
        yo(i,j) = radRaw(i-1) * SIN(azimRaw(j)*2._dp*pi_d/24._dp - pi_d)
     ENDDO
  ENDDO
  CALL GSL_Interpolation_2D(x(nThetaEquator,:,2:nzeta),y(nThetaEquator,:,2:nzeta), &
                            PhiIono(:,2:nzeta), xo(:,:), yo(:,:), Epot_Cart(:,:), GSLerr)


  ! Make sure the potential values are reasonable
  if ((maxval(abs(Epot_Cart)).gt.150000.0).or.(isnan(sum(Epot_Cart)))) then
     if (verbose) write(*,*) 'Bad Electric Field data, using previous time step'
     EOut_II = VT
  else
     EOut_II = Epot_Cart
  endif
  if (verbose) write(*,'(1x,a,2F10.2)') 'EOut_II max and min', maxval(EOut_II), minval(EOut_II)
  if (verbose) write(*,*) ''

  if (WritePotential) then
     open(UNITTMP_,FILE=trim(PathRamOut)//RamFileName('potential','dat',TimeRamNow))
     write(UNITTMP_,*) nR+1,nT
     do j = 1, nT
      do i = 1, nR+1
       write(UNITTMP_,*) xo(i,j), yo(i,j), EOut_II(i,j)
      enddo
     enddo
     close(UNITTMP_)
  endif

  DEALLOCATE(Epot_Cart, xo, yo)
  RETURN

end subroutine ram_get_electric

!==================================================================================================
SUBROUTINE ionospheric_potential
  !!!! Module Variables
  use ModRamMain,      ONLY: PathRamOut
  use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed
  use ModRamParams,    ONLY: IsComponent, electric, UseSWMFFile, NameOmniFile, &
                             WritePotential
  use ModRamConst,     ONLY: RE
  use ModRamVariables, ONLY: Kp
  use ModScbMain,      ONLY: iConvE
  use ModScbGrids,     ONLY: npsi, nzeta, nzetap
  use ModScbVariables, ONLY: phiiono, x, y, z, r0Start, dPhiIonodAlpha, &
                             dPhiIonodBeta, f, fzet, zetaVal, rhoVal, tilt, &
                             nThetaEquator
  use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
  use ModSceVariables, ONLY: IONO_NORTH_Phi, IONO_NORTH_JR, PhiIono_Weimer
  use ModRamCouple,    ONLY: SWMFIonoPot_II, nIeTheta, nIePhi
  !!!! Module Subroutine/Functions
  use ModRamGSL,    ONLY: GSL_Interpolation_2D, GSL_Derivs
  use ModRamFunctions, ONLY: RamFileName
  use ModScbIO,     ONLY: Write_Ionospheric_Potential
  use ModSceRun,    ONLY: sce_run
  !!!! Share Modules
  use ModTimeConvert, ONLY: n_day_of_year
  use ModIOUnit, ONLY: UNITTMP_
  !!!! Extra Modules
  use w05
  !!!! NR Modules
  use nrtype, ONLY: DP, pi_d


  implicit none

  integer :: doy, GSLerr, i, j, k, ierr, iYear_l, &
             iDoy_l, iHour_l, iMin_l, iLines, isec_l, imsec_l, imonth_l, iday_l, &
             AL_l, SymH_l, sgn
  REAL(DP) :: radius, angle, bzimf_l, bndylat, byimf_l, pdyn_l, Nk_l, &
              Vk_l, bTot_l, bximf_l, vx_l, vy_l, vz_l, t_l, AVS, VT
  REAL(DP), ALLOCATABLE :: colat(:), lon(:), phiIonoRaw(:,:), dPhiIonodRho(:,:), &
                           dPhiIonodZeta(:,:), colatGrid(:,:), lonGrid(:,:), latGrid(:,:)
  REAL(DP), EXTERNAL :: EpotVal, BoundaryLat
  CHARACTER(LEN = 100) :: header
  LOGICAL :: UseAL

  !================================================================================================
  ALLOCATE(dPhiIonodRho(npsi, nzeta+1), dPhiIonodZeta(npsi, nzeta+1), &
           colatGrid(npsi,nzeta+1), lonGrid(npsi,nzeta+1), latGrid(npsi,nzeta+1))
  dPhiIonodRho = 0.0; dPhiIonodZeta = 0.0; colatGrid=0.0; lonGrid = 0.0; latGrid = 0.0

  ! Set latitude/colatitude and longitude grids
  DO k = 2, nzeta
     DO j = 1, npsi
        radius = SQRT((x(1,j,k))**2 + y(1,j,k)**2)
        angle = ATAN2(y(1,j,k), x(1,j,k)) ! Angle from noon
        IF (angle < 0.) angle = angle + 2.*pi_d
        colatGrid(j,k) = 0.5_dp*pi_d - ASIN(-z(1,j,k)/r0Start)
        latGrid(j,k) = 0.5 * pi_d - colatGrid(j,k) ! In radians
        lonGrid(j,k) = angle
     END DO
  END DO

  ! Periodicity
  colatGrid(:,nzetap) = colatGrid(:,2)
  lonGrid(:,nzetap) = lonGrid(:,2) + 2.*pi_d
  latGrid(:,nzetap) = latGrid(:,2)
  colatGrid(:,1) = colatGrid(:,nzeta)
  lonGrid(:,1) = lonGrid(:,nzeta) - 2.*pi_d
  latGrid(:,1) = latGrid(:,nzeta)

  SELECT CASE (electric)

  CASE('IESC') ! interpolate SWMF iono potentials on 3Dcode ionospheric grid
     !call con_stop('IESC SWMF electric field not currently working')
     ! If coupled to SWMF, we don't need to open any files and crap.
     ! Initialize arrays and fill with correct values.
     ! IM_wrapper ensures that, despite the IE grid, we always
     ! have values from and including 10 to 80 degrees colat.
     if (.not. allocated(PhiIonoRaw)) then
        allocate(PhiIonoRaw(nIeTheta, nIePhi))
        PhiIonoRaw = 0.0
     endif
     if(.not. allocated(colat) .and. .not. allocated(lon)) then
        allocate(colat(nIeTheta), lon(nIePhi))
        colat = 0.0
        lon   = 0.0
        do i=2, nIePhi ! Longitude goes from 0 to 360.
           lon(i) = lon(i-1) + 2.0_dp*pi_d/real(nIePhi-1)
        end do
        do i=2, nIeTheta ! Colat goes from 10 to 80.
           colat(i) = colat(i-1) + (7.0_dp/18.0_dp) *pi_d/real(nIeTheta-1)
        end do
        colat = colat + pi_d/18.0_dp
     end if
     ! Plug SWMF/IE potential into PhiIonoRaw; then we're all set.
     PhiIonoRaw = SwmfIonoPot_II

     CALL GSL_Interpolation_2D(colat, lon, PhiIonoRaw, colatGrid(1:npsi,2:nzetap), &
                            lonGrid(1:npsi,2:nzetap), PhiIono(1:npsi,2:nzetap), GSLerr)

  CASE('RSCE')
     ! Potential from self-consistent electric field calculated in the ionosphere (Yu. 2016 August)
     if(.not.allocated(PhiIonoRaw))then
        allocate(PhiIonoRaw(IONO_nTheta, IONO_nPsi))
        PhiIonoRaw = 0.0
     end if
     if(.not.allocated(colat).and..not.allocated(lon))then
        allocate(colat(IONO_nTheta), lon(IONO_nPsi))
        colat = 0.0
        lon   = 0.0
        do i=2, IONO_nPsi ! Longitude goes from 0 to 360.
           lon(i) = lon(i-1) + 2.0_dp*pi_d/real(IONO_nPsi-1)
        end do
        do i=2, IONO_nTheta ! Colat goes from 0 to 90.
           colat(i) = colat(i-1) +  0.5*pi_d/real(IONO_nTheta-1)
        end do
     end if
     ! Plug self-consisent/IE potential into PhiIonoRaw (only northern)
     CALL sce_run
     PhiIonoRaw = Iono_North_Phi

     CALL GSL_Interpolation_2D(colat, lon, PhiIonoRaw, colatGrid(1:npsi,2:nzetap), &
                               lonGrid(1:npsi,2:nzetap), PhiIono(1:npsi,2:nzetap), GSLerr)

     deallocate(PhiIonoRaw, colat, lon)

  CASE('WESC', 'W5SC') ! Potential by calling Weimer function
     OPEN(UNITTMP_, FILE=NameOmniFile, status = 'UNKNOWN', action = 'READ')
     ! Omni files used for the SWMF have a different format from those used in older RAM simulations
     if (UseSWMFFile) then
        UseAL = .false.
        Read_SWMF_file: DO
           read(UNITTMP_,*) iyear_l, imonth_l, iday_l, ihour_l, imin_l, isec_l, imsec_l, &
                            bximf_l, byimf_l, bzimf_l, vx_l, vy_l, vz_l, nk_l, t_l
           if ((TimeRamNow%iYear.eq.iyear_l).and.(TimeRamNow%iMonth.eq.imonth_l).and. &
               (TimeRamNow%iDay.eq.iday_l).and. &
               (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
              bTot_l = SQRT(bximf_l**2 + byimf_l**2 + bzimf_l**2)
              vk_l   = SQRT(vx_l**2 + vy_l**2 + vz_l**2)
              EXIT Read_SWMF_file
           END IF
        END DO Read_SWMF_file
     else
        UseAL = .true.
        Read_OMNI_file: DO
           READ(UNITTMP_,*) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, &
                            bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
           doy = n_day_of_year(TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay)
           if ((TimeRamNow%iYear.eq.iyear_l).and.(doy.eq.idoy_l).and. &
               (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
              EXIT Read_OMNI_file
           END IF
        END DO Read_OMNI_file
     endif
     CLOSE(UNITTMP_)

     if (bzimf_l.lt.-20) bzimf_l = -20
     if (bzimf_l.gt.20) bzimf_l = 20
     angle = 180._dp/pi_d*ACOS(bzimf_l/bTot_l)
     latGrid = 180._dp/pi_d * latGrid ! Transform latitude to degrees for Weimer EpotVal function
     lonGrid = 12._dp/pi_d * lonGrid ! Transform lonGrid to hours (MLT) for Weimer EpotVal function
     DO k = 2, nzeta
        DO j = 1, npsi
           IF (lonGrid(j,k) <= 12._dp) THEN
              lonGrid(j,k) = lonGrid(j,k) + 12._dp
           ELSE
              lonGrid(j,k) = lonGrid(j,k) - 12._dp
           END IF
        END DO
     END DO

     ! Choose which Weimer model to use 
     if (electric == 'WESC') then
        write(*,*) 'Year ', ' DOY ', ' Hour ', ' Min ', ' angle ', '  bT  ', ' tilt ', '    vT    ', '  Nd  '
        WRITE(*,'(1x,I4,2x,I3,2x,I2,4x,I2,3x,3F6.2,F10.2,F6.2)') &
              iYear_l, iDoy_l, iHour_l, iMin_l, angle, bTot_l, tilt, Vk_l, Nk_l

        CALL SetModel(angle,bTot_l,tilt,Vk_l,Nk_l,REAL(AL_l,dp),UseAL)
        DO k = 2, nzeta
           DO j = 1, npsi
              PhiIono(j,k) = 1.E3 * EpotVal(latGrid(j,k), lonGrid(j,k)) !!! Achtung RAM needs it in volts !!!
              bndylat = BoundaryLat(lonGrid(j,k))
           END DO
        END DO
     elseif (electric == 'W5SC') then
        write(*,*) 'Year ', ' DOY ', ' Hour ', ' Min ', '  bY  ', '  bZ  ', ' tilt ', '    vT    ', '  Nd  '
        WRITE(*,'(1x,I4,2x,I3,2x,I2,4x,I2,1x,3F6.2,F10.2,F6.2)') &
              iYear_l, iDoy_l, iHour_l, iMin_l, byimf_l, bzimf_l, tilt, Vk_l, Nk_l

        CALL SetModel05(byimf_l, bzimf_l, tilt, Vk_l, Nk_l)
        AVS=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE  ! Voll-Stern parameter
        DO k = 2, nzeta
           DO j = 1, npsi
              radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
              angle = ATAN2(y(nThetaEquator,j,k), x(nThetaEquator,j,k))
              VT = AVS*(radius*Re)**2*SIN(angle+pi_d)
              call EpotVal05(latGrid(j,k), lonGrid(j,k), 0.0_dp, PhiIono(j,k))
              PhiIono(j,k) = PhiIono(j,k) * 1.0e3   ! ram needs it in Volts
              ! If the weimer model gives us 0 (because we are outside of its
              ! valid domain) we set the the potential to the Voll-Stern
              ! potential
              if (abs(PhiIono(j,k)) < 0.001) then
               PhiIono(j,k) = VT
              elseif (abs(VT) > abs(PhiIono(j,k))) then
               sgn = abs(PhiIono(j,k))/PhiIono(j,k)
               PhiIono(j,k) = sgn*abs(VT)
              endif
           END DO
        END DO
     endif

  CASE DEFAULT
     DEALLOCATE(dPhiIonodRho, dPhiIonodZeta, colatGrid, lonGrid,latGrid)
     call con_stop('Unrecognized electric field parameter')

  END SELECT

  ! Add corotation potential - only if RAM takes all E-field info from SCB
  PhiIono(:,:) = PhiIono(:,:) - 2._dp*pi_d*0.31_dp*6.4**2 * 1.E3_dp &
                               /(24._dp*36._dp*SQRT(x(nThetaEquator,:,:)**2+y(nThetaEquator,:,:)**2))

  ! Azimuthal periodicity
  PhiIono(:,nzetap) = PhiIono(:,2)
  PhiIono(:,1) = PhiIono(:,nzeta)

  CALL GSL_Derivs(rhoVal, zetaVal(2:nzeta), PhiIono(:,2:nzeta), &
                  dPhiIonodRho(:,2:nzeta), dPhiIonodZeta(:,2:nzeta), GSLerr)

  dPhiIonodRho(:,nzetap) = dPhiIonodRho(:,2)
  dPhiIonodRho(:,1) = dPhiIonodRho(:,nzeta)
  dPhiIonodZeta(:,nzetap) = dPhiIonodZeta(:,2)
  dPhiIonodZeta(:,1) = dPhiIonodZeta(:,nzeta)

  DO j = 1, npsi
     dPhiIonodAlpha(j,1:nzetap) = dPhiIonodRho(j,1:nzetap) / f(j)
  END DO
  DO k = 1, nzetap
     dPhiIonodBeta(1:npsi, k) = dPhiIonodZeta(1:npsi,k) / fzet(k)
  END DO

  if (WritePotential) call Write_Ionospheric_Potential

  DEALLOCATE(dPhiIonodRho, dPhiIonodZeta, colatGrid, lonGrid,latGrid)
  RETURN

END SUBROUTINE ionospheric_potential

!============================================================================

END MODULE ModRamEField
