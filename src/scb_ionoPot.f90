SUBROUTINE ionospheric_potential
  !!!! Module Variables
  use ModRamMain,      ONLY: PathSwmfOut
  use ModRamTiming,    ONLY: TimeRamNow
  use ModRamParams,    ONLY: IsComponent, electric
  use ModRamCouple,    ONLY: SwmfIonoPot_II, nIePhi, nIeTheta
  use ModRamIndices,   ONLY: NameOmniFile
  use ModScbMain,      ONLY: iConvE
  use ModScbGrids,     ONLY: npsi, nzeta, nthe, nzetap
  use ModScbVariables, ONLY: phiiono, x, y, z, r0Start, dPhiIonodAlpha, &
                             dPhiIonodBeta, f, fzet, zetaVal, rhoVal, tilt
  !!!! Module Subroutine/Functions
  use ModScbSpline, ONLY: Spline_2D_derivs, Spline_2D_periodic
  !!!! Share Modules
  use ModTimeConvert, ONLY: n_day_of_year
  use ModIOUnit, ONLY: UNITTMP_
  !!!! Extra Modules
  use w05
  !!!! NR Modules
  use nrtype, ONLY: DP

  IMPLICIT NONE

  integer :: doy
  INTEGER :: i, ierralloc, j, j1, k1, k, ierr, ierrDom, idealerr
  REAL(DP) :: dPhiIonodRho(npsi, nzeta+1), dPhiIonodZeta(npsi, nzeta+1), &
       colatGrid(npsi,nzeta+1), lonGrid(npsi,nzeta+1), latGrid(npsi,nzeta+1)
  REAL(DP) :: radius, angle, lineData(3)
  REAL(DP) :: thangle, zangle, byimf, bzimf_l
  REAL(DP) :: bt, swvel, swden, alindex, bndylat
  REAL(DP), PARAMETER :: tiny = 1.E-6_dp
  INTEGER :: ier, iCount_neighbor, iDomain
  !C  INTEGER, PARAMETER :: mlat_range = 81, mlon_range = 181
  INTEGER, PARAMETER :: mlat_range = 36, mlon_range = 91
  INTEGER :: iTimeArr(0:24), dstArr(0:24)
  INTEGER :: iYear_l, iDoy_l, iHour_l, iMin_l, iLines
  REAL(DP) :: bzimfArr_l, byimf_l, pdyn_l, Nk_l, Vk_l, bTot_l
  INTEGER :: AL_l, SymH_l
  REAL(DP), ALLOCATABLE :: colat(:), lon(:), phiIonoRaw(:,:)
  REAL(DP), EXTERNAL :: EpotVal, BoundaryLat
  CHARACTER(LEN = 3) :: zeroChars
  CHARACTER(LEN = 15) :: StringDateTime
  CHARACTER(LEN = 100) :: header
  LOGICAL :: UseAL

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

  CASE('SWMF') ! interpolate SWMF iono potentials on 3Dcode ionospheric grid
     if(IsComponent) then
        ! If coupled to SWMF, we don't need to open any files and crap.
        ! Initialize arrays and fill with correct values.
        ! IM_wrapper ensures that, despite the IE grid, we always
        ! have values from and including 10 to 80 degrees colat.
        if(.not. allocated(PhiIonoRaw)) &
             allocate(PhiIonoRaw(nIeTheta, nIePhi))
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
        
     else
        ! We now return you to your regularly scheduled 3DEQ code.
        IF (.NOT. ALLOCATED(PhiIonoRaw)) &
             ALLOCATE(PhiIonoRaw(mlat_range,mlon_range), STAT = ierralloc)
        IF (.NOT. ALLOCATED(colat) .AND. .NOT. ALLOCATED(lon)) &
             ALLOCATE(colat(mlat_range), lon(mlon_range), STAT = ierralloc)
        
        ! Use current time to build IE file name.
        write(StringDateTime, '(i4.4, 2i2.2, "_",3i2.2)') &
             TimeRamNow%iYear, TimeRamNow%iMonth,  TimeRamNow%iDay, &
             TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond
        ! Open IE file for reading.
        write(*,*)'Reading iono data at ', &
             trim(PathSwmfOut)//'/'//'IE_'//StringDateTime//'_SCB.in'
        OPEN(UNITTMP_, file = trim(PathSwmfOut)//'/'//'IE_'//StringDateTime//'_SCB.in', &
             action = 'READ', status='OLD')

        DO i = 1, 9
           READ(UNITTMP_, '(A)') HEADER
        END DO
        DO j = 1, mlon_range
           DO i = 1, mlat_range
              READ(UNITTMP_,*) lineData
              colat(i) = lineData(1)*pi_d/180._dp
              lon(j) = lineData(2)*pi_d/180._dp
              PhiIonoRaw(i,j) = lineData(3)
           END DO
        END DO
        CLOSE(UNITTMP_)

     end if

     ! Now there should be no difference between what 
     ! we do for coupled and non-coupled cases.

     !   lon(mlon_range) = 2._dp*pi_d + lon(1)
     !   phiIonoRaw(:,mlon_range) = phiIonoRaw(:,1)

     CALL Spline_2D_periodic(colat, lon, PhiIonoRaw, colatGrid(1:npsi,2:nzetap), lonGrid(1:npsi,2:nzetap), &
          PhiIono(1:npsi,2:nzetap), iDomain)

     IF (iDomain > 0) THEN
        PRINT*, 'Stop; problem with PhiIono domain; iDomain = ', iDomain
        STOP
     END IF
     
  CASE('WESC') ! Potential by calling Weimer function 
     UseAL = .TRUE.
     !C OPEN(UNITTMP_, FILE='./omni_5min_Sep2005.txt', status = 'UNKNOWN', action = 'READ')
     OPEN(UNITTMP_, FILE=NameOmniFile, status = 'UNKNOWN', action = 'READ')
     i = 0
     ierr = 0
     DO WHILE (ierr == 0) 
        READ(UNITTMP_, *, IOSTAT=ierr) header
        i = i+1
     END DO
     iLines = i-1
     !C PRINT*, 'IP: the file has ', iLines, ' lines.'
     
     ! Rewind file
     REWIND(UNITTMP_)

     PRINT*, 'IP: Year, Day, Hour, Min, Bt, By, Bz, V, N, Pdyn, AL, SYMH'
     Read_data_from_file: DO i = 1, iLines
        READ(UNITTMP_,*) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
        doy = n_day_of_year(TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay)
        if ((TimeRamNow%iYear.eq.iyear_l).and.(doy.eq.idoy_l).and. &
            (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
           WRITE(*,*) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l           
           EXIT Read_data_from_file
        END IF
        ! ACHTUNG this assumes the same starting file in omni file, and same cadence of E-field update as the B-field update in RAM-SCB !!! 
     END DO Read_data_from_file
     CLOSE(UNITTMP_)


     ! 101  FORMAT(A11, 1X, I2, A3, 1X, F5.1, 1X, F5.1, 1X, F5.1, 1X, F5.0, 1X, F5.2, 1X, I5)
101  FORMAT(2I4, 2I3, 3F8.2, F8.1, F7.2, F6.2, 2I6)
     ! 1 Year                          I4
     ! 2 Day                           I4        
     ! 3 Hour                          I3        
     ! 4 Minute                        I3        
     ! 5 Field magnitude average, nT   F8.2      
     ! 6 BY, nT (GSM)                  F8.2      
     ! 7 BZ, nT (GSM)                  F8.2      
     ! 8 Speed, km/s                   F8.1      
     ! 9 Proton Density, n/cc          F7.2      
     ! 10 Flow pressure, nPa            F6.2      
     ! 11 AL-index, nT                  I6        
     ! 12 SYM/H, nT                     I6        

     !C PRINT*, 'iono_potential: iHour_l = ', iHour_l

     !C Bt = SQRT(byimf**2 + bzimf_l**2)
     angle = 180._dp/pi_d*ACOS(bzimf_l/bTot_l)
     WRITE(*,'(A,4F12.2,I6)') 'IP: Bt, angle, SW Vel., SW Dens., AL = ', bTot_l, angle, Vk_l, Nk_l, AL_l
     CALL SetModel(angle,bTot_l,tilt,Vk_l,Nk_l,REAL(AL_l,dp),UseAL); CALL FLUSH(6) 
     latGrid = 180._dp/pi_d * latGrid ! Transform latitude to degrees for Weimer EpotVal function
     lonGrid = 12._dp/pi_d * lonGrid ! Transform lonGrid to hours (MLT) for Weimer EpotVal function
     ! Translate to correct MLT
     DO k = 2, nzeta
        DO j = 1, npsi
           IF (lonGrid(j,k) <= 12._dp) THEN
              lonGrid(j,k) = lonGrid(j,k) + 12._dp
           ELSE
              lonGrid(j,k) = lonGrid(j,k) - 12._dp
           END IF
        END DO
     END DO

     Weimer_potential_loop: DO k = 2, nzeta
        DO j = 1, npsi
           IF (ABS(pdyn_l-99.99_dp) < 1E-3_dp) THEN
              PRINT*, 'IP: not calling Weimer model, bad SW data'
              PRINT*, ' '
              EXIT Weimer_potential_loop ! Bad SW data, do not call Weimer model
           END IF
           PhiIono(j,k) = 1.E3 * EpotVal(latGrid(j,k), lonGrid(j,k)) !!! Achtung RAM needs it in volts !!!
           bndylat = BoundaryLat(lonGrid(j,k))
    !C       PRINT*, 'j, k, lat, mlt, bndylat, Phi_Weimer(j,k) = ', j, k, REAL(latGrid(j,k),sp), &
    !C           REAL(lonGrid(j,k),sp), real(bndylat), 1.E-3*REAL(PhiIono(j,k),sp)
           ! Add corotation potential - only if RAM takes all E-field info from the equilibrium code 
           ! PhiIono(j,k) = PhiIono(j,k) - 2._dp*pi_d*0.31_dp*6.4**2 * 1.E3_dp/(24._dp*36._dp*SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2))
        END DO
     END DO Weimer_potential_loop

  CASE('W5SC') ! Potential by calling Weimer 2005 function 'W5SC'
     open(Unittmp_, file=NameOmniFile, status = 'unknown',action='read')
     i = 0
     ierr = 0
     do while (ierr == 0)
        read(unittmp_, *, iostat=ierr) header
        i = i + 1
     end do
     ilines = i-1

     rewind(unittmp_)
     print*, 'IP: year, day, hour, min, bt, by, bz, v, n, pdyn, al. symh'
     read_data_from_file1: do i=1,ilines
        read(unittmp_,*) iyear_l,idoy_l,ihour_l, imin_l, btot_l, byimf_l, bzimf_l, vk_l, nk_l, pdyn_l, al_l, symh_l
        doy = n_day_of_year(TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay)
        if ((TimeRamNow%iYear.eq.iyear_l).and.(doy.eq.idoy_l).and. &
            (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
           write(*,*)iyear_l,idoy_l,ihour_l, imin_l, btot_l, byimf_l, bzimf_l, vk_l, nk_l, pdyn_l, al_l, symh_l
           exit read_data_from_file1
        end if
     end do read_data_from_file1
     close(unittmp_)
     angle = 180._dp/pi_d*ACOS(bzimf_l/bTot_l)
     WRITE(*,'(A,4F12.2,I6)') 'IP: Bt, angle, SW Vel., SW Dens., AL = ', bTot_l, angle, Vk_l, Nk_l, AL_l
     CALL SetModel05(byimf_l, bzimf_l, tilt, Vk_l, Nk_l)
     CALL FLUSH(6) 
     latGrid = 180._dp/pi_d * latGrid ! Transform latitude to degrees for Weimer EpotVal function
     lonGrid = 12._dp/pi_d * lonGrid ! Transform lonGrid to hours (MLT) for Weimer EpotVal function
     ! Translate to correct MLT
     DO k = 2, nzeta
        DO j = 1, npsi
           IF (lonGrid(j,k) <= 12._dp) THEN
              lonGrid(j,k) = lonGrid(j,k) + 12._dp
           ELSE
              lonGrid(j,k) = lonGrid(j,k) - 12._dp
           END IF
        END DO
     END DO

     Weimer_potential_loop1: DO k = 2, nzeta
        DO j = 1, npsi
           IF (ABS(pdyn_l-99.99_dp) < 1E-3_dp) THEN
              PRINT*, 'IP: not calling Weimer model, bad SW data'
              PRINT*, ' '
              EXIT Weimer_potential_loop1 ! Bad SW data, do not call Weimer model
           END IF
           call EpotVal05(latGrid(j,k), lonGrid(j,k), 0.0_dp, PhiIono(j,k)) 
           PhiIono(j,k) = PhiIono(j,k) * 1.0e3   ! ram needs it in Volts
           ! Add corotation potential - only if RAM takes all E-field info from the equilibrium code 
           ! PhiIono(j,k) = PhiIono(j,k) - 2._dp*pi_d*0.31_dp*6.4**2 * 1.E3_dp/(24._dp*36._dp*SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2))
        END DO

     END DO Weimer_potential_loop1


  CASE DEFAULT
     STOP 'ionospheric_potential: problem.'

  END SELECT

  ! Azimuthal periodicity
  PhiIono(:,nzetap) = PhiIono(:,2)
  PhiIono(:,1) = PhiIono(:,nzeta)

  IF (iConvE /= 1) RETURN

  CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), PhiIono(:,1:nzeta), dPhiIonodRho(:,1:nzeta), dPhiIonodZeta(:,1:nzeta))

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

  RETURN

END SUBROUTINE ionospheric_potential



