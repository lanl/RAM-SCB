SUBROUTINE TM03_DensPres(nbc, rad_max, NT, local_times, densPres)

  ! Tsyganenko-Mukai statistical plasma sheet model - gives out density and pressure
  ! from Tsyganenko, N. and T. Mukai, Tail plasma sheet models derived from Geotail particle data
  ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 108, NO. A3, 1136, doi:10.1029/2002JA009707, 2003

  ! nbc is passed from ram, the index of the boundary update

  USE ModRamMain, ONLY : Real8_, PI
  USE ModConst, ONLY : cProtonMass
  USE ModIoUnit, ONLY: UNITTMP_
  USE ModRamIndices, ONLY : NameOmniFile

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbc, NT
  REAL(Real8_), INTENT(IN) :: rad_max
  REAL(Real8_), INTENT(IN) :: local_times(NT)

  REAL(Real8_), INTENT(IN OUT) :: densPres(2,NT)  ! Should be passed as size (2,NT)

  REAL(Real8_), PARAMETER :: A1 = 0.057, A2 = 0.524, A3 = 0.0908, A4 = 0.527, &
       A5 = 0.078, A6 = -4.422, A7 = -1.533, A8 = -1.217, A9 = 2.54, &
       A10 = 0.32, A11 = 0.754, A12 = 1.048, A13 = -0.074, A14 = 1.015 
  REAL(Real8_), PARAMETER :: A1_d = -0.159, A2_d = 0.608, A3_d = 0.5055, A4_d = 0.0796, &
       A5_d = 0.2746, A6_d = 0.0361, A7_d = -0.0342, A8_d = -0.7935, A9_d = 1.162, &
       A10_d = 0.4756, A11_d = 0.7117
  REAL, PARAMETER :: cMass2Num = 1.0E6 * cProtonMass
  REAL(Real8_) :: Bperp, rhoStar, rhoStar10, F, FStar, theta, &
       NStar, PStar, VStar, BSStar, BNStar
  REAL(Real8_), ALLOCATABLE :: phi(:)
  INTEGER :: iYear_l, iDoy_l, iHour_l, iMin_l, i, ierr, iLines
  REAL(Real8_) :: bTot_l, byimf_l, bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
  CHARACTER(LEN = 100) :: header

  IF (.NOT.ALLOCATED(phi)) ALLOCATE(phi(NT), stat=ierr)

  rhoStar = 0.1 * rad_max 

  ! Obtains solar wind quantities from omni text file
  OPEN(UNITTMP_, FILE=NameOmniFile, status = 'UNKNOWN', action = 'READ')
  i = 0
  ierr = 0
  DO WHILE (ierr == 0) 
     READ(UNITTMP_, *, IOSTAT=ierr) header
     i = i+1
  END DO
  iLines = i-1
  PRINT*, 'TM03_DensPres: the file has ', iLines, ' lines.'

  ! Rewind file
  REWIND(UNITTMP_)

  PRINT*, 'TM03: NT, nbc = ', NT, nbc
  PRINT*, 'TM03: Year, Day, Hour, Min, Bt, By, Bz, V, N, Pdyn, AL, SYMH'
  Read_data_from_file: DO i = 0, iLines-1
     READ(UNITTMP_, *) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
     IF (i == nbc) THEN
        WRITE(*, *) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
        EXIT Read_data_from_file !!! ACHTUNG this assumes omni file has same data cadence (5 min.) as boundary update !!!
     END IF
  END DO Read_data_from_file
  CLOSE(UNITTMP_)
101 FORMAT(2I4, 2I3, 3F8.2, F8.1, F7.2, F6.2, 2I6)

  NStar = 0.1 * Nk_l ! Normalization
  PStar = pdyn_l/3.
  VStar = Vk_l/500.
  BSStar = MAX(-bzimf_l, 0._Real8_) / 5. ! Southward component of IMF Bz
  BNStar = MAX(bzimf_l, 0._Real8_) / 5. ! Northward component of IMF Bz

  theta = ATAN2(byimf_l,bzimf_l) ! IMF clock angle
  IF (theta < 0.0) theta = theta + 2*PI
  Bperp = SQRT(byimf_l**2 + bzimf_l**2) ! BSW perpendicular to Sun-Earth axis
  F = ABS(Bperp)*SQRT(SIN(0.5*theta))
  FStar = 0.2*F

  ! Transform local times into azimuthal angle, - atan(Y/X) ! TM03 takes angle positive in the dusk sector
  phi = PI/12. * local_times
  phi = 2*PI - phi ! TM03 dawn-dusk symmetric (only sin**2 terms)

  !C phi = ATAN2(xEqGsm, yEqGsm) ! azimuthal angle, atan(y/x)

  ! print*, 'TM03: size(densPres) = ', size(densPres)
  PRINT*, 'TM03: pdyn, theta, F, NStar, PStar, VStar, BSStar, BNStar = ', REAL(pdyn_l), REAL(theta), REAL(F), &
       REAL(NStar), REAL(PStar), REAL(VStar), REAL(BSStar), REAL(BNStar)

  ! Plasma sheet density
  densPres(1,:) = (A1_d + A2_d*NStar**A10_d + A3_d*BNStar + A4_d*VStar*BSStar) * rhoStar**A8_d + &
       (A5_d*NStar**A11_d + A6_d*BNStar + A7_d*VStar*BSStar) * rhoStar**A9_d * (SIN(phi))**2 
  !C print*, 'TM03: density at 9 RE in cm-3 = ', densPres(1,:)
  PRINT*, 'TM03: Max. density at 9 RE in cm-3 = ', MAXVAL(densPres(1,:))
  densPres(1,:) = cMass2Num * densPres(1,:) ! Transformation to mass density (kg/m^3)

  ! Plasma sheet pressure
  densPres(2,:) = A1*rhoStar**A6 + A2*PStar**A11*rhoStar**A7 + A3*FStar**A12*rhoStar**A8 + &
       (A4*PStar**A13*EXP(-A9*rhoStar) + A5*FStar**A14*EXP(-A10*rhoStar))*(SIN(phi))**2
  !C print*, 'TM03: pressure at 9 RE in nPa = ', densPres(2,:)
  PRINT*, 'TM03: Max. pressure at 9 RE in nPa = ', MAXVAL(densPres(2,:))
  densPres(2,:) = 1.E-9 * densPres(2,:) ! Transformation to Pascals

  IF (ALLOCATED(phi)) DEALLOCATE(phi, stat=ierr)


  RETURN 

END SUBROUTINE TM03_DensPres


