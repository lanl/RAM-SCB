!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Takes magnetic boundary conditions and initial guess for X, Y, Z 
!  from file obtained by (empirical or MHD) magnetic field model tracing
!  Copyright (c) 2016, Los Alamos National Security, LLC
!  All rights reserved.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

SUBROUTINE Computational_domain(i_comp_domain)
  USE ModRamMain, ONLY: PathScbIn, IsComponent
  USE nrtype
  USE Module1
  USE ezcdf
  USE ModIoUnit, ONLY: UNITTMP_
  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: i_comp_domain

  INTEGER                 :: i, j, k, ierr, ilat, ilon, ifld, inBlank
  INTEGER                 :: ntheTsyg, npsiTsyg, nzetaTsyg
  INTEGER                 :: netcdfId, netcdfid2
  INTEGER, DIMENSION(3)   :: dimlens 
  INTEGER, DIMENSION(3)   :: dimlens2
  REAL(DP)                :: er, error, xpsi_in, xpsi_out
  CHARACTER(LEN = 10)      :: statusBound
  CHARACTER(LEN = 8)      :: pdyn_char, byimf_char, bzimf_char, dst_char, r00_char, xpsiin_char, &
       xpsiout_char, xpsiin_char_fl, xpsiout_char_fl, xpsiin_char_ce, xpsiout_char_ce, &
       constz_char, consttheta_char, nthe_char, npsi_char, nzeta_char
  CHARACTER               :: xtype, xtype2, xtype3, xtype4
  CHARACTER*6            :: kpChar
  CHARACTER*10           :: timeChar
  CHARACTER*4            :: KpTrunc
  CHARACTER*1             :: KpFl, KpCe
  REAL(SP)                :: tiltSP = 0.0
  REAL(DP)                :: tiltDP = 0.0
  REAL(DP)                :: KpReal, ratioFl=1
  REAL(SP), ALLOCATABLE :: xFl(:,:,:), yFl(:,:,:), zFl(:,:,:), xCe(:,:,:), yCe(:,:,:), zCe(:,:,:) ! SP for SWMF
  REAL(DP), ALLOCATABLE :: xFlDbl(:,:,:), yFlDbl(:,:,:), zFlDbl(:,:,:), &
       xCeDbl(:,:,:), yCeDbl(:,:,:), zCeDbl(:,:,:) ! DP for empirical models
  REAL(DP) :: xpsiin_fl, xpsiin_ce, xpsiout_fl, xpsiout_ce
  INTEGER             :: iKpFl, iKpCe
  CHARACTER*300        :: fileNameTsyga, fileNameTsyga2
  CHARACTER*2 :: minuteInteger
  CHARACTER*4 :: ST3

  ALLOCATE(xFl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
  ALLOCATE(yFl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
  ALLOCATE(zFl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
  ALLOCATE(xFlDbl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
  ALLOCATE(yFlDbl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
  ALLOCATE(zFlDbl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)

  IF (iPressureChoice == 5 .OR. iPressureChoice == 6 .OR. iPressureChoice == 7 .OR. &
       iPressureChoice == 8 .OR. iPressureChoice == 11) THEN ! RAM coupling case, or Roeder pressures
     SELECT CASE (i_comp_domain) 
     CASE(1)
        !C CALL  cdf_open (netcdfId, 't96_config.cdf', 'r')
        fileNameTsyga = TRIM(ADJUSTL('t96_2007_324_'//TRIM(ADJUSTL(Hour))//'00.cdf'))
        inBlank = INDEX(fileNameTsyga, ' ') - 1
        CALL cdf_open(netcdfId, fileNameTsyga(1:inBlank), 'r')
        PRINT*, 'Calling magnetic boundary file: ', fileNameTsyga(1:inBlank)
        blendGlobal = blendInitial
     CASE(2)  ! TS04 Boundary conditions
        ALLOCATE(xCe(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
        ALLOCATE(yCe(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
        ALLOCATE(zCe(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
        !C IF (rank == 0) PRINT*, 'Hour is ', HourDecimal; CALL FLUSH(6)
        !   ratioFl = 1. - (KpReal - iKpFl) ! The closest to Fl, the higher (toward 1) the ratio must be
        !   ratioFl = 1. - (rHour-FLOOR(rHour)) 
        ratioFl = 1. ! T04S boundaries not implemented with linear interpolation yet
        ! WRITE(KpFl, '(I1)') iKpFl
        ! WRITE(KpCe, '(I1)') iKpCe

        WRITE(minuteInteger, '(I2)') MOD(iST3,12)*5+2
        IF (MOD(iST3,12)*5+2 < 10) minuteInteger(1:1)= '0'
        fileNameTsyga = TRIM(ADJUSTL(PathScbIn))//&
             '/t04s_config_'//event//'_'//Day//'_'//HourInteger//'_'//minuteInteger//'.cdf'

        IF (rank == 0) THEN
           !       print*, 'iHour, Hour, Hour+1 = ', iHour, HourInteger, ' ', HourIntegerP1
           inBlank = INDEX(fileNameTsyga, ' ') - 1
           PRINT*, 'Calling magnetic boundary file: ', fileNameTsyga(1:inBlank) ! , fileNameTsyga2
           !           WRITE(*, '(A, F10.2)') 'ratioFl = ', ratioFl
        END IF
        blendGlobalInitial = blendInitial
        blendGlobal = blendInitial
        !C        IF (blendInitial >= 0.05_dp) blendGlobal = MIN(0.5_dp, 0.01+(blendGlobalInitial-0.05) * 18./ REAL(iKpFl**2+iKpCe**2,dp)) ! decrease at high Kp
        !C        print*, 'Computational_domain: KpFl, KpCe, blendGlobal = ', iKpFl, iKpCe, blendGlobal
        CALL  cdf_open (netcdfId, TRIM(fileNameTsyga), 'r')

     CASE(20)  ! SWMF-based boundary conditions, 5-min resolution
        if(IsComponent)then
           call build_scb_init
           fileNameTsyga = 'SWMF_config_test.cdf'
        else
           WRITE(minuteInteger, '(I2)') MOD(iST3,12)*5 + 0
           IF (MOD(iST3,12)*5+0 < 10) minuteInteger(1:1) = '0'
           WRITE(ST3, '(I4)') iST3
           IF (iST3 < 10) ST3(3:3) = '0'
           IF (iST3 < 100) ST3(2:2) = '0'
           IF (iST3 < 1000) ST3(1:1) = '0'
           !C        PRINT*, 'ST3 = ', ST3
           fileNameTsyga = TRIM(ADJUSTL(prefixIn))//'SWMF_config_'//event//'_'//ST3//'.cdf'
        end if
        inBlank = INDEX(fileNameTsyga, ' ') - 1
        IF (rank == 0) PRINT*, 'Calling magnetic boundary file: ', fileNameTsyga(1:inBlank) ! , fileNameTsyga2
        blendGlobalInitial = blendInitial
        blendGlobal = blendInitial
        CALL  cdf_open (netcdfId, fileNameTsyga(1:inBlank), 'r')
     CASE(3)
        ALLOCATE(xCeDbl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
        ALLOCATE(yCeDbl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
        ALLOCATE(zCeDbl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
        ! To improve smoothness takes an average boundary from blending a "floor" T89 with a "ceiling" T89 field, based on the 
        ! non-integer Kp at that particular time(hour)
        IF (iPressureChoice /=7 ) THEN ! Vania pressure file
           OPEN(UNITTMP_, file=TRIM(fileNamePressure), status = 'old')
           READ(UNITTMP_, '(A40, 2X, A6)') timeChar, kpChar
           CLOSE(UNITTMP_)
           KpTrunc = kpChar(4:6)
           READ(KpTrunc, '(F3.1)') KpReal
           iKpFl = FLOOR(KpReal)
           iKpCe = CEILING(kpReal)
           ratioFl = 1. - (KpReal - iKpFl) ! The closest to Fl, the higher (toward 1) the ratio must be
           !C IF (iKpFl < 1) iKpFl = 1  ! Not needed - T89c has KP = 0 case
           IF (iKpFl > 6) iKpFl = 6
           IF (iKpCe > 6) iKpCe = 6
        END IF

        !C   ratioFl = 1 ! For GEM Challenge, only use 1 T89-based file at a time

        IF (iPressureChoice == 7) THEN ! Roeder, statistical (take Kp = 1)
           iKpFl = 1
           iKpCe = 1
           ratioFl = 1
        END IF

        WRITE(KpFl, '(I1)') iKpFl
        WRITE(KpCe, '(I1)') iKpCe

        fileNameTsyga = TRIM(ADJUSTL(PathScbIn))//'/'//'t89_config_KP'//KpFl//'.cdf'
        fileNameTsyga2 = TRIM(ADJUSTL(PathScbIn))//'/'//'t89_config_KP'//KpCe//'.cdf'
        inBlank = INDEX(fileNameTsyga, ' ') - 1
        IF (rank == 0) THEN
           PRINT*, 'Will call files: ', fileNameTsyga(1:inBlank), ' and ', fileNameTsyga2(1:inBlank)
           WRITE(*, '(A, 2 F10.2)') 'KpReal, ratioFl = ', KpReal, ratioFl
        END IF
        blendGlobalInitial = blendInitial
        !C        blendGlobal = MIN(0.5_dp, 0.03+(blendGlobalInitial-0.05) * 18./ REAL(iKpFl**2+iKpCe**2,dp)) ! decrease at high Kp
        blendGlobal = blendGlobalInitial
        CALL  cdf_open (netcdfId, TRIM(fileNameTsyga), 'r')
        CALL  cdf_open (netcdfId2, TRIM(fileNameTsyga2), 'r')
     CASE(4)
        CALL  cdf_open (netcdfId, trim(PathScbIn)//'dipole_config.cdf', 'r')
        blendGlobalInitial = blendInitial
        blendGlobal = blendGlobalInitial
     END SELECT
  ELSE
     SELECT CASE (i_comp_domain) 
     CASE(1)
        CALL  cdf_open (netcdfId, 't96_config.cdf', 'r')
        blendGlobal = blendInitial  
     CASE(2)
        CALL  cdf_open (netcdfId, 't01_storm_config.cdf', 'r')
     CASE(3)
        CALL  cdf_open (netcdfId, TRIM(ADJUSTL(prefixIn))//'t89_config_KP2.cdf', 'r')
        blendGlobal = blendInitial
     CASE(4)
        CALL  cdf_open (netcdfId, 'dipole_config.cdf', 'r')
        blendGlobalInitial = blendInitial
        blendGlobal = blendGlobalInitial
     END SELECT
  END IF

  !  IF (rank == 0) PRINT*, 'External file has been called.'; CALL FLUSH(6)

  CALL  cdf_inquire (netcdfId, 'statusBC', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'pdyn', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'byimf', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'bzimf', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'DST', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'xpsi_in', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'xpsi_out', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'rStart', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'Constz', dimlens, xtype)
  CALL  cdf_inquire(netcdfId, 'Consttheta', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'ntheta', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'npsi', dimlens, xtype)
  CALL  cdf_inquire (netcdfId, 'nzeta', dimlens, xtype)
  CALL  cdf_inquire(netcdfId, 'tilt', dimlens, xtype4)
  CALL  cdf_inquire(netcdfId, 'wTsyg', dimlens, xtype3)
  CALL  cdf_inquire (netcdfId, 'xTsyg', dimlens2, xtype2)
  CALL  cdf_inquire (netcdfId, 'yTsyg', dimlens2, xtype2)
  CALL  cdf_inquire (netcdfId, 'zTsyg', dimlens2, xtype2)

  !C  CALL  cdf_read (netcdfId, 'statusBC', statusBound) 
  IF (i_comp_domain /= 20) THEN
     CALL  cdf_read (netcdfId, 'pdyn', pdyn_char)
     CALL  cdf_read (netcdfId, 'byimf', byimf_char)
     CALL  cdf_read (netcdfId, 'bzimf', bzimf_char)
     CALL  cdf_read (netcdfId, 'DST', dst_char)
     CALL  cdf_read (netcdfId, 'wTsyg', wTsyg)
  END IF
  CALL  cdf_read (netcdfId, 'xpsi_in', xpsiin_char_fl)
  CALL  cdf_read (netcdfId, 'xpsi_out', xpsiout_char_fl)
  CALL  cdf_read (netcdfId, 'rStart', r00_char)
  CALL  cdf_read (netcdfId, 'Constz', constz_char)
  CALL  cdf_read (netcdfId, 'Consttheta', consttheta_char)
  CALL  cdf_read (netcdfId, 'ntheta', nthe_char)
  CALL  cdf_read (netcdfId, 'npsi', npsi_char)
  CALL  cdf_read (netcdfId, 'nzeta', nzeta_char)
  IF (i_comp_domain == 20) THEN
     CALL  cdf_read (netcdfId, 'tilt', tiltSP) 
  ELSE
     CALL cdf_read(netcdfId, 'tilt', tiltDP)
  END IF
  PRINT*, 'tilt=',tiltDP

  IF (iCompDomain == 20) THEN ! SP
     CALL  cdf_read (netcdfId, 'xTsyg', xFl)
     CALL  cdf_read (netcdfId, 'yTsyg', yFl)
     CALL  cdf_read (netcdfId, 'zTsyg', zFl)
  ELSE
     CALL  cdf_read (netcdfId, 'xTsyg', xFlDbl)
     CALL  cdf_read (netcdfId, 'yTsyg', yFlDbl)
     CALL  cdf_read (netcdfId, 'zTsyg', zFlDbl)
  END IF

  IF ((i_comp_domain == 3) .AND. isotropy == 0) THEN
     CALL  cdf_read (netcdfId2, 'xTsyg', xCeDbl)
     CALL  cdf_read (netcdfId2, 'yTsyg', yCeDbl)
     CALL  cdf_read (netcdfId2, 'zTsyg', zCeDbl)
     CALL  cdf_read (netcdfId2, 'xpsi_in', xpsiin_char_ce)
     CALL  cdf_read (netcdfId2, 'xpsi_out', xpsiout_char_ce)
     CALL cdf_close(netcdfId2)
  END IF

  CALL cdf_close(netcdfId)
  ! transform characters to variables using internal file read
  IF (i_comp_domain /= 20) THEN
     READ(pdyn_char, '(SP, F6.2)') p_dyn
     READ(byimf_char, '(SP, F6.2)') by_imf
     READ(bzimf_char, '(SP, F6.2)') bz_imf
     READ(dst_char, '(SP, F6.1)') dst_global
     pdynGlobal = p_dyn
     byimfGlobal = by_imf
     bzimfGlobal = bz_imf
  END IF
  READ(xpsiin_char_fl, '(SP, F6.2)') xpsiin_fl
  READ(xpsiout_char_fl, '(SP, F6.2)') xpsiout_fl
  IF (i_comp_domain == 3 .AND. isotropy==0) THEN
     READ(xpsiin_char_ce, '(SP, F6.2)') xpsiin_ce
     READ(xpsiout_char_ce, '(SP, F6.2)') xpsiout_ce
  END IF
  READ(r00_char, '(SP, F6.2)') r0Start
  READ(constz_char, '(SP, F6.2)') constz
  READ(consttheta_char, '(SP, F6.2)') constTheta
  READ(nthe_char, '(I4)') ntheTsyg
  READ(npsi_char, '(I4)') npsiTsyg
  READ(nzeta_char, '(I4)') nzetaTsyg

  IF ((iPressureChoice==5 .OR. iPressureChoice==6 .OR. iPressureChoice == 11) .AND. &
       (i_comp_domain == 3)) THEN ! If coupled with RAM
     x = ratioFl*xFlDbl + (1._dp - ratioFl)*xCeDbl
     y = ratioFl*yFlDbl + (1._dp - ratioFl)*yCeDbl
     z = ratioFl*zFlDbl + (1._dp - ratioFl)*zCeDbl
     xpsiin = ratioFl*xpsiin_fl + (1._dp - ratioFl)*xpsiin_ce
     xpsiout = ratioFl*xpsiout_fl + (1._dp - ratioFl)*xpsiout_ce
  ELSE IF (i_comp_domain==2 .OR. i_comp_domain==4) THEN ! 5-min. cadence TS04, no interpolations; or dipole (double prec.) 
     x = ratioFl*xFlDbl 
     y = ratioFl*yFlDbl 
     z = ratioFl*zFlDbl 
     xpsiin = ratioFl*xpsiin_fl 
     xpsiout = ratioFl*xpsiout_fl 
  ELSE
     x = xFl
     y = yFl
     z = zFl
     xpsiin = xpsiin_fl
     xpsiout = xpsiout_fl
  END IF

  IF (ntheTsyg /= nthe .OR. npsiTsyg /= npsi .OR. nzetaTsyg /= nzeta) THEN
     PRINT*, 'ntheTsyg, nthe = ', ntheTsyg, nthe
     PRINT*, 'npsiTsyg, npsi = ', npsiTsyg, npsi
     PRINT*, 'nzetaTsyg, nzeta = ', nzetaTsyg, nzeta
     !   STOP 'Number of grid points does not correspond to data file.'
  END IF

  IF (rank == 0) THEN
     !C IF (rank == 0 .AND. iHourChoice == 1) THEN
     ! Write to screen values of parameters
     !  WRITE(*, '(A12, A11)') 'statusBC = ', TRIM(ADJUSTL(statusBound))
     IF (i_comp_domain == 20) THEN
        WRITE(*, '(A12, F10.2)') 'tilt = ', tiltSP
        tilt = REAL(tiltSP, DP) ! tilt is defined as DP in Module1
     ELSE
        tilt = tiltDP
     END IF
     IF(i_comp_domain /= 3 .AND. i_comp_domain /= 20) THEN
        WRITE(*, '(A12, F11.3)') 'p_dyn = ', p_dyn
        WRITE(*, '(A12, F11.3)') 'by_imf = ', by_imf
        WRITE(*, '(A12, F11.3)') 'bz_imf = ', bz_imf
        WRITE(*, '(A12, F11.3)') 'dst = ', dst_global
        IF (i_comp_domain == 2 .OR. i_comp_domain == 21) WRITE(*, '(A12, 6F9.2)') 'W(1:6) = ', wTsyg(1:6)
     END IF
     IF (i_comp_domain == 3) THEN
        WRITE(*, '(A35, 3F11.3)') 'xpsiinFl, xpsiinCe, xpsiin = ', xpsiin_fl, xpsiin_ce, xpsiin
        WRITE(*, '(A35, 3F11.3)') 'xpsioutFl, xpsioutCe, xpsiout = ', xpsiout_fl, xpsiout_ce, xpsiout
     END IF
     WRITE(*, '(A12, F11.3)') 'r00 = ', r0Start
     WRITE(*, '(A12, F11.3)') 'constZ = ', constz
     ! constTheta = 0._dp ! Temporary here
     WRITE(*, '(A12, F11.3)') 'constTheta = ', constTheta
     PRINT*, ' '

     !     WRITE(*, *) 'Number of grid points:'
     !     WRITE(*, '(A8, I3)') 'nthe = ', ntheTsyg
     !     WRITE(*, '(A8, I3)') 'npsi = ', npsiTsyg
     !     WRITE(*, '(A8, I3)') 'nzeta = ', nzetaTsyg

     PRINT*, 'Domain starts from a sphere of radius ', INT(r0Start), ' R_E'
     PRINT*, 'MIN, MAX x: ', MINVAL(x), MAXVAL(x)
     PRINT*, ' '
  END IF

  IF (ALLOCATED(xFl)) DEALLOCATE(xFl, STAT = ierr)
  IF (ALLOCATED(yFl)) DEALLOCATE(yFl, STAT = ierr)
  IF (ALLOCATED(zFl)) DEALLOCATE(zFl, STAT = ierr)
  IF (ALLOCATED(xCe)) DEALLOCATE(xCe, STAT = ierr)
  IF (ALLOCATED(yCe)) DEALLOCATE(yCe, STAT = ierr)
  IF (ALLOCATED(zCe)) DEALLOCATE(zCe, STAT = ierr)
  IF (ALLOCATED(xFlDbl)) DEALLOCATE(xFlDbl, STAT = ierr)
  IF (ALLOCATED(yFlDbl)) DEALLOCATE(yFlDbl, STAT = ierr)
  IF (ALLOCATED(zFlDbl)) DEALLOCATE(zFlDbl, STAT = ierr)
  IF (ALLOCATED(xCeDbl)) DEALLOCATE(xCeDbl, STAT = ierr)
  IF (ALLOCATED(yCeDbl)) DEALLOCATE(yCeDbl, STAT = ierr)
  IF (ALLOCATED(zCeDbl)) DEALLOCATE(zCeDbl, STAT = ierr)

  tsygcorrect = 0  ! If = 1, forces N-S symmetry, taking southern mapping data

  IF (tsygcorrect /= 1) THEN
     x(:,:,1) = x(:,:,nzeta)
     y(:,:,1) = y(:,:,nzeta)
     z(:,:,1) = z(:,:,nzeta)
     x(:,:,nzeta+1) = x(:,:,2)
     y(:,:,nzeta+1) = y(:,:,2)
     z(:,:,nzeta+1) = z(:,:,2)
     RETURN
  END IF

  DO j = 1, npsi
     DO k = 2, nzeta
        DO i = 1, nthe/2
           x(nthe + 1 - i,j,k) = x(i, j, k)
           y(nthe + 1 - i,j,k) = y(i, j, k)
           z(nthe + 1 - i,j,k)= - z(i, j, k)  ! Force N-S symmetry on Tsyganenko data, using south data 
           ! (the tracing can lead to slight asymmetry, even with zero tilt)
        END DO
        z(nthe/2 + 1, j, k) = 0.0_dp
     END DO
  END DO
  x(:,:,1) = x(:,:,nzeta)
  y(:,:,1) = y(:,:,nzeta)
  z(:,:,1) = z(:,:,nzeta)
  x(:,:,nzeta+1) = x(:,:,2)
  y(:,:,nzeta+1) = y(:,:,2)
  z(:,:,nzeta+1) = z(:,:,2)


  RETURN

END SUBROUTINE Computational_domain
