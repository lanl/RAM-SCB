MODULE ModScbIO
  
  use ezcdf
  
  implicit none
  save
  
  contains

!=============================================================================!
!============================= INPUT ROUTINES ================================!
!=============================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Takes magnetic boundary conditions and initial guess for X, Y, Z 
!  from file obtained by (empirical or MHD) magnetic field model tracing
!  Copyright (c) 2016, Los Alamos National Security, LLC
!  All rights reserved.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  SUBROUTINE computational_domain
    !!!! Module Variables  
    USE ModRamVariables, ONLY: Kp
    use ModRamParams,    ONLY: IsComponent, NameBoundMag, boundary
    USE ModScbMain,      ONLY: PathScbIn, blendInitial, tsygcorrect
    USE ModScbParams,    ONLY: decreaseConvAlphaMin, decreaseConvPsiMin, &
                               decreaseConvAlphaMax, decreaseConvPsiMax
    USE ModScbGrids,     ONLY: npsi, nthe, nzeta
    USE ModScbVariables, ONLY: by_imf, bz_imf, dst_global, p_dyn, wTsyg, tilt, constZ, &
                               constTheta, xpsiin, xpsiout, r0Start, byimfglobal, &
                               bzimfglobal, pdynglobal, blendGlobal, blendGlobalInitial, &
                               x, y, z
    !!!! Module Subroutines/Functions
    use ModScbCouple, ONLY: build_scb_init 
    !!!! Share Modules
    USE ModIoUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, SP
  
    IMPLICIT NONE
  
    INTEGER :: i, j, k, ierr, ilat, ilon, ifld, inBlank, iKpFl, iKpCe, &
               ntheTsyg, npsiTsyg, nzetaTsyg, netcdfId, netcdfId2
  
    INTEGER, DIMENSION(3) :: dimlens, dimlens2
  
    REAL(SP) :: tiltSP = 0.0
    REAL(DP) :: tiltDP = 0.0
    REAL(DP) :: er, error, xpsi_in, xpsi_out, KpReal, xpsiin_fl, &
                xpsiin_ce, xpsiout_fl, xpsiout_ce
  
    REAL(SP), ALLOCATABLE :: xFl(:,:,:), yFl(:,:,:), zFl(:,:,:), &
                             xCe(:,:,:), yCe(:,:,:), zCe(:,:,:) ! SP for SWMF
    REAL(DP), ALLOCATABLE :: xFlDbl(:,:,:), yFlDbl(:,:,:), zFlDbl(:,:,:), &
                             xCeDbl(:,:,:), yCeDbl(:,:,:), zCeDbl(:,:,:) ! DP for empirical models
  
    CHARACTER(LEN=300) :: fileNameTsyga, fileNameTsyga2
    CHARACTER(LEN=10)  :: statusBound, timeChar
    CHARACTER(LEN=8)   :: pdyn_char, byimf_char, bzimf_char, dst_char, r00_char, &
                          xpsiin_char, xpsiout_char, xpsiin_char_fl, xpsiout_char_fl, &
                          xpsiin_char_ce, xpsiout_char_ce, constz_char, consttheta_char, &
                          nthe_char, npsi_char, nzeta_char
    CHARACTER(LEN=6)   :: kpChar
    CHARACTER(LEN=4)   :: KpTrunc, ST3
    CHARACTER(LEN=2)   :: minuteInteger
    CHARACTER(LEN=1)   :: KpFl, KpCe
    CHARACTER          :: xtype, xtype2, xtype3, xtype4
  
    REAL(DP) :: ratioFl=1

    ALLOCATE(xFl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    ALLOCATE(yFl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
    ALLOCATE(zFl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
    ALLOCATE(xFlDbl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    ALLOCATE(yFlDbl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
    ALLOCATE(zFlDbl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
    ALLOCATE(xCeDbl(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    ALLOCATE(yCeDbl(SIZE(y,1), SIZE(y,2), SIZE(y,3)), STAT = ierr)
    ALLOCATE(zCeDbl(SIZE(z,1), SIZE(z,2), SIZE(z,3)), STAT = ierr)
  
!    IF (boundary == 'SWMF' .OR. boundary == 'LANL') THEN
       SELECT CASE (NameBoundMag)
       CASE('SWMF')  ! SWMF-based boundary conditions, 5-min resolution
          call build_scb_init
          fileNameTsyga = 'SWMF_config_test.cdf'
          inBlank = INDEX(fileNameTsyga, ' ') - 1
          PRINT*, 'Calling magnetic boundary file: ', fileNameTsyga(1:inBlank)
          blendGlobalInitial = blendInitial
          blendGlobal = blendInitial
          CALL  cdf_open (netcdfId, fileNameTsyga(1:inBlank), 'r')
       CASE('DIPL')
          CALL  cdf_open (netcdfId, trim(PathScbIn)//'dipole_config.cdf', 'r')
          blendGlobalInitial = blendInitial
          blendGlobal = blendGlobalInitial
       CASE('DIPS')
          CALL  cdf_open (netcdfId, trim(PathScbIn)//'dipole_config.cdf', 'r')
          blendGlobalInitial = blendInitial
          blendGlobal = blendGlobalInitial
       CASE('T89C')
          blendGlobalInitial = blendInitial
          blendGlobal = blendGlobalInitial
          KpReal = Kp
          iKpFl = min(FLOOR(KpReal),6)
          iKpCe = min(CEILING(KpReal),6)
          ratioFl = 1. - (KpReal - iKpFl)
          WRITE(KpFl, '(I1)') iKpFl
          WRITE(KpCe, '(I1)') iKpCe
          fileNameTsyga = TRIM(ADJUSTL(PathScbIn))//'/t89_config_KP'//KpFl//'.cdf'
          fileNameTsyga2 = TRIM(ADJUSTL(PathScbIn))//'/t89_config_KP'//KpCe//'.cdf'
          inBlank = INDEX(fileNameTsyga, ' ') - 1
          PRINT*, 'Will call files: ', fileNameTsyga(1:inBlank), ' and ', fileNameTsyga2(1:inBlank)
          WRITE(*, '(A, 2 F10.2)') 'KpReal, ratioFl = ', KpReal, ratioFl
          call cdf_open (netcdfId, TRIM(fileNameTsyga), 'r')
          call cdf_open (netcdfId2, TRIM(fileNameTsyga2), 'r')
       END SELECT
!    END IF
  
    CALL cdf_inquire(netcdfId, 'statusBC',   dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'pdyn',       dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'byimf',      dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'bzimf',      dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'DST',        dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'xpsi_in',    dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'xpsi_out',   dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'rStart',     dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'Constz',     dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'Consttheta', dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'ntheta',     dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'npsi',       dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'nzeta',      dimlens,  xtype)
    CALL cdf_inquire(netcdfId, 'tilt',       dimlens,  xtype4)
    CALL cdf_inquire(netcdfId, 'wTsyg',      dimlens,  xtype3)
    CALL cdf_inquire(netcdfId, 'xTsyg',      dimlens2, xtype2)
    CALL cdf_inquire(netcdfId, 'yTsyg',      dimlens2, xtype2)
    CALL cdf_inquire(netcdfId, 'zTsyg',      dimlens2, xtype2)
  
    IF (NameBoundMag /= 'SWMF') THEN
       CALL cdf_read(netcdfId, 'pdyn',  pdyn_char)
       CALL cdf_read(netcdfId, 'byimf', byimf_char)
       CALL cdf_read(netcdfId, 'bzimf', bzimf_char)
       CALL cdf_read(netcdfId, 'DST',   dst_char)
       CALL cdf_read(netcdfId, 'wTsyg', wTsyg)
    END IF
    CALL cdf_read(netcdfId, 'xpsi_in',    xpsiin_char_fl)
    CALL cdf_read(netcdfId, 'xpsi_out',   xpsiout_char_fl)
    CALL cdf_read(netcdfId, 'rStart',     r00_char)
    CALL cdf_read(netcdfId, 'Constz',     constz_char)
    CALL cdf_read(netcdfId, 'Consttheta', consttheta_char)
    CALL cdf_read(netcdfId, 'ntheta',     nthe_char)
    CALL cdf_read(netcdfId, 'npsi',       npsi_char)
    CALL cdf_read(netcdfId, 'nzeta',      nzeta_char)
  
    IF (NameBoundMag == 'SWMF') THEN ! SP
       CALL cdf_read(netcdfId, 'xTsyg', xFl)
       CALL cdf_read(netcdfId, 'yTsyg', yFl)
       CALL cdf_read(netcdfId, 'zTsyg', zFl)
       CALL cdf_read(netcdfId, 'tilt',  tiltSP)
    ELSE
       CALL cdf_read(netcdfId, 'xTsyg', xFlDbl)
       CALL cdf_read(netcdfId, 'yTsyg', yFlDbl)
       CALL cdf_read(netcdfId, 'zTsyg', zFlDbl)
       CALL cdf_read(netcdfId, 'tilt',  tiltDP)
    END IF
  
    if (NameBoundMag == 'T89C') THEN
       CALL cdf_read(netcdfId2, 'xTsyg',    xCeDbl)
       CALL cdf_read(netcdfId2, 'yTsyg',    yCeDbl)
       CALL cdf_read(netcdfId2, 'zTsyg',    zCeDbl)
       CALL cdf_read(netcdfId2, 'xpsi_in',  xpsiin_char_ce)
       CALL cdf_read(netcdfId2, 'xpsi_out', xpsiout_char_ce)
       CALL cdf_close(netcdfId2)
    end if
  
    CALL cdf_close(netcdfId)
  
    ! transform characters to variables using internal file read
    IF (NameBoundMag /= 'SWMF') THEN
       READ(pdyn_char,  '(SP, F6.2)') p_dyn
       READ(byimf_char, '(SP, F6.2)') by_imf
       READ(bzimf_char, '(SP, F6.2)') bz_imf
       READ(dst_char,   '(SP, F6.1)') dst_global
       pdynGlobal  = p_dyn
       byimfGlobal = by_imf
       bzimfGlobal = bz_imf
    END IF
    READ(xpsiin_char_fl,  '(SP, F6.2)') xpsiin_fl
    READ(xpsiout_char_fl, '(SP, F6.2)') xpsiout_fl
    if (NameBoundMag == 'T89C') THEN
       READ(xpsiin_char_ce,  '(SP, F6.2)') xpsiin_ce
       READ(xpsiout_char_ce, '(SP, F6.2)') xpsiout_ce
    end if
    READ(r00_char,        '(SP, F6.2)') r0Start
    READ(constz_char,     '(SP, F6.2)') constz
    READ(consttheta_char, '(SP, F6.2)') constTheta
    READ(nthe_char,  '(I4)') ntheTsyg
    READ(npsi_char,  '(I4)') npsiTsyg
    READ(nzeta_char, '(I4)') nzetaTsyg
  
    IF (NameBoundMag=='DIPL' .OR. NameBoundMag=='DIPS') THEN ! 5-min. cadence TS04, no interpolations; or dipole (double prec.) 
       x = ratioFl*xFlDbl
       y = ratioFl*yFlDbl
       z = ratioFl*zFlDbl
       xpsiin  = ratioFl*xpsiin_fl
       xpsiout = ratioFl*xpsiout_fl
    ELSE IF (NameBoundMag == 'T89C') THEN
       x = ratioFl*xFlDbl + (1._dp - ratioFl)*xCeDbl
       y = ratioFl*yFlDbl + (1._dp - ratioFl)*yCeDbl
       z = ratioFl*zFlDbl + (1._dp - ratioFl)*zCeDbl
       xpsiin  = ratioFl*xpsiin_fl + (1._dp - ratioFl)*xpsiin_ce
       xpsiout = ratioFl*xpsiout_fl + (1._dp - ratioFl)*xpsiout_ce
    ELSE
       x = xFl
       y = yFl
       z = zFl
       xpsiin  = xpsiin_fl
       xpsiout = xpsiout_fl
    END IF
  
    IF (NameBoundMag == 'SWMF') THEN
       tilt = REAL(tiltSP, DP) ! tilt is defined as DP in Module1
    ELSE
       tilt = tiltDP
    END IF
       IF(NameBoundMag /= 'T89C' .AND. NameBoundMag /= 'SWMF') THEN
          WRITE(*, '(A12, F11.3)') 'p_dyn = ', p_dyn
          WRITE(*, '(A12, F11.3)') 'by_imf = ', by_imf
          WRITE(*, '(A12, F11.3)') 'bz_imf = ', bz_imf
          WRITE(*, '(A12, F11.3)') 'dst = ', dst_global
       END IF
       IF (NameBoundMag == 'T89C') THEN
          WRITE(*, '(A35, 3F11.3)') 'xpsiinFl, xpsiinCe, xpsiin = ', xpsiin_fl, xpsiin_ce, xpsiin
          WRITE(*, '(A35, 3F11.3)') 'xpsioutFl, xpsioutCe, xpsiout = ', xpsiout_fl, xpsiout_ce, xpsiout
       END IF
       WRITE(*, '(A12, F11.3)') 'r00 = ', r0Start
       WRITE(*, '(A12, F11.3)') 'constZ = ', constz
       ! constTheta = 0._dp ! Temporary here
       WRITE(*, '(A12, F11.3)') 'constTheta = ', constTheta
       PRINT*, ' '
  
       PRINT*, 'Domain starts from a sphere of radius ', INT(r0Start), ' R_E'
       PRINT*, 'MIN, MAX x: ', MINVAL(x), MAXVAL(x)
       PRINT*, ' '
  
    IF (ALLOCATED(xFl)) DEALLOCATE(xFl, STAT = ierr)
    IF (ALLOCATED(yFl)) DEALLOCATE(yFl, STAT = ierr)
    IF (ALLOCATED(zFl)) DEALLOCATE(zFl, STAT = ierr)
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
             z(nthe + 1 - i,j,k) = -z(i, j, k)  ! Force N-S symmetry on Tsyganenko data, using south data 
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

!=============================================================================!
!============================= OUTPUT ROUTINES ===============================!
!=============================================================================!
SUBROUTINE Write_ionospheric_potential

  use ModRamTiming, ONLY: TimeRamElapsed

  use ModScbMain,      ONLY: prefixOut
  use ModScbGrids,     ONLY: npsi, nzeta
  use ModScbVariables, ONLY: x, y, z, PhiIono, dPhiIonodAlpha, dPhiIonodBeta, &
                             alphaVal, psiVal, nThetaEquator, bnormal

  USE nrtype
  USE netcdf

  IMPLICIT NONE

  CHARACTER*500 :: filename

  INTEGER :: alphaid, betaid, timeid, alphavarid, betavarid, &
       xeqid, yeqid, xionoid, yionoid, phiionoid, dphiionodalphaid, &
       dphiionodbetaid, timevarid, ncid

  integer :: START(1), COUNT(1)
  integer :: START1D(2), COUNT1D(2)
  INTEGER :: START2D(3), COUNT2D(3) ! For 2-D arrays (+ time)

  INTEGER, SAVE :: iCALLIP = 0

  REAL :: time(1)

  time(1) = TimeRamElapsed

  START = (/iCALLIP+1/)
  COUNT = (/1/)

  START1D = (/1,iCALLIP+1/)
  START2D = (/1,1,iCALLIP+1/)


  fileName = TRIM(ADJUSTL(prefixOut))//'ionospheric_potential.nc'

  First_time_call : IF(iCALLIP == 0) THEN
     CALL check (nf90_create(filename, nf90_clobber, ncid))

     ! Define dimensions
     CALL check(nf90_def_dim(ncid, 'alpha', npsi, alphaid))
     CALL check(nf90_def_dim(ncid, 'beta', nzeta, betaid))
     CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, timeid))

     ! Define variables
     CALL check(nf90_def_var(ncid, 'alpha', nf90_float,(/alphaid,timeid/),alphavarid))
     CALL check(nf90_put_att(ncid,alphavarid,'title','Magnetic flux-like Euler potential'))

     CALL check(nf90_def_var(ncid, 'beta', nf90_float,(/betaid,timeid/),betavarid))
     CALL check(nf90_put_att(ncid,betavarid,'title','Azimuthal angle-like Euler potential'))

     CALL check(nf90_def_var(ncid, 'time', nf90_float,timeid,timevarid))
     CALL check(nf90_put_att(ncid,timevarid,'title','Time'))

     CALL check(nf90_def_var(ncid, 'xEq', nf90_float, (/alphaid,betaid,timeid/),xeqid))
     CALL check(nf90_put_att(ncid,xeqid,'title','2D array of xEq locations '))

     CALL check(nf90_def_var(ncid, 'yEq', nf90_float, (/alphaid,betaid,timeid/),yeqid))
     CALL check(nf90_put_att(ncid,yeqid,'title','2D array of yEq locations '))

     CALL check(nf90_def_var(ncid, 'xIono', nf90_float, (/alphaid,betaid,timeid/),xionoid))
     CALL check(nf90_put_att(ncid,xionoid,'title','2D array of xIono locations '))

     CALL check(nf90_def_var(ncid, 'yIono', nf90_float, (/alphaid,betaid,timeid/),yionoid))
     CALL check(nf90_put_att(ncid,yionoid,'title','2D array of yIono locations '))

     CALL check(nf90_def_var(ncid, 'PhiIono', nf90_float, (/alphaid,betaid,timeid/),phiionoid))
     CALL check(nf90_put_att(ncid,phiionoid,'title','2D array of phiIono values'))

     CALL check(nf90_def_var(ncid, 'dPhiIonodAlpha', nf90_float, (/alphaid,betaid,timeid/),dphiionodalphaid))
     CALL check(nf90_put_att(ncid,dphiionodalphaid,'title','2D array of dPhi/dAlpha values'))

     CALL check(nf90_def_var(ncid, 'dPhiIonodBeta', nf90_float, (/alphaid,betaid,timeid/),dphiionodbetaid))
     CALL check(nf90_put_att(ncid,dphiionodbetaid,'title','2D array of dPhi/dBeta values'))

    ! End define mode
     CALL check(nf90_enddef(ncid))

  ELSE ! Open existing NetCDF file
     CALL check(nf90_open(filename, nf90_write, ncid))

     CALL check( nf90_inq_dimid(ncid, 'alpha', alphaid))
     CALL check( nf90_inq_dimid(ncid, 'beta', betaid))
     CALL check( nf90_inq_dimid(ncid, 'time', timeid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'alpha',alphavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'beta', betavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',timevarid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xEq', xeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yEq', yeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'xIono', xionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yIono', yionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'PhiIono', phiionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'dPhiIonodAlpha', dphiionodalphaid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'dPhiIonodBeta', dphiionodbetaid))

  END IF First_time_call

 ! Write mode - write at all times
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START1D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START1D))
  CALL check(nf90_put_var(ncid, xeqid, REAL(x(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yeqid, REAL(y(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, xionoid, REAL(x(1,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yionoid, REAL(y(1,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, phiionoid, REAL(PhiIono(1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, dphiionodalphaid, REAL(dPhiIonodAlpha(1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, dphiionodbetaid, REAL(dPhiIonodBeta(1:npsi,1:nzeta)),START2D))

  CALL check(nf90_put_var(ncid, timevarid, time, START, COUNT))

  CALL check(nf90_close(ncid))

  iCALLIP = iCALLIP+1

  RETURN

CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT*, 'STATUS = ', status
       PRINT *, TRIM(nf90_strerror(status))
       STOP 2
    END IF
  END SUBROUTINE check

END SUBROUTINE Write_ionospheric_potential

!==============================================================================
  ! Previously test_Convergence_anisotropic
  SUBROUTINE Write_convergence_anisotropic
  !!!! Module Variables
  use ModRamTiming,    ONLY: TimeRamNow
  use ModScbMain,      ONLY: prefixOut
  use ModScbParams,    ONLY: isotropy
  USE ModScbGrids,     ONLY: nthe, npsi, nzeta, dt, dr, dpPrime
  USE ModScbVariables, ONLY: thetaVal, rhoVal, zetaVal, x, y, z, &
                             jacobian, normDiff, normGradP, GradZetaSq, &
                             GradThetaGradZeta, GradRhoGradTheta, GradRhoSq, &
                             GradRhoGradZeta, ppar, pper, nThetaEquator, &
                             normJxB, f, fzet, nZetaMidnight, pnormal, &
                             dPPerdRho, dPPerdZeta, dPPerdTheta, bnormal, &
                             pjconst, dPdAlpha, dPdPsi, vecd, vec1, vec2, &
                             vec3, vec4, vec6, vec7, vec8, vec9, vecr, vecx, &
                             alfa, psi, fp, alphaVal, psiVal
  !!!! Module Subroutine/Function
  use ModRamFunctions, ONLY: RamFileName
  use ModScbEquation,  ONLY: metric, metrica, newk, newj
  use ModScbSpline,    ONLY: Spline_coord_derivs
  !!!! Share Modules
  USE ModIoUnit, ONLY: UNITTMP_
  !!!! NR Modules
  use nrtype,    ONLY: DP

  IMPLICIT NONE

  INTEGER :: i, j, k, id, ierr, idealerr
  CHARACTER(len=200) :: FileName

  REAL(DP) :: normDiffRel, volume, bf(nthe,npsi,nzeta+1), bsq(nthe,npsi,nzeta+1), &
              distance(nthe,npsi,nzeta+1)
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: derivXTheta, derivXRho, derivXZeta, &
       derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq, derivBsqTheta, derivBsqRho, derivBsqZeta, &
       derivNU1, derivNU2
  ! gradRhoSq, gradRhoGradZeta are global

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jGradRhoPartialTheta, derivjGradRhoPartialTheta, &
       jGradRhoPartialZeta, derivjGradRhoPartialZeta, jGradRho, jGradRhoFactor, jGradZetaPartialRho, &
       derivjGradZetaPartialRho, jGradZetaPartialTheta, derivjGradZetaPartialTheta, jGradZeta, &
       jGradZetaFactor, jGradThetaPartialRho, derivjGradThetaPartialRho, jGradThetaPartialZeta, &
       derivjGradThetaPartialZeta, jGradTheta, jGradThetaFactor, phiToroid, derivPhiRho, derivPhiZeta, &
       derivPhiTheta, derivDiffPTheta

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jCrossBUpRho, jCrossBUpZeta, jCrossBUpTheta, &
       derivjCrossBUpRho, derivjCrossBUpZeta, derivjCrossBUpTheta, jCrossBMinusGradPSq, &
       jCrossBMinusGradPMod, jCrossBSq, jCrossB, gradPSq, gradP

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jrrInt, jrr, jzzInt, jzz, jrtInt, jrt, jztInt, jzt, &
       rhoCompSq, zetaCompSq, thetaCompSq, curlJCrossBSq, curlJCrossB

  REAL(DP), DIMENSION(npsi) :: rtemp, xtemp
!  LOGICAL, EXTERNAL :: isnand ! Intrinsic for Portland Group Fortran

  !**********************************************************************************************************!

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), &
                           derivXTheta, derivXRho, derivXZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), &
                           derivYTheta, derivYRho, derivYZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), &
                           derivZTheta, derivZRho, derivZZeta)
  ! Now I have all the point derivatives

  ! Time to build the Jacobian
  jacobian = derivXRho * (derivYZeta * derivZTheta - derivYTheta * derivZZeta) &
           + derivXZeta * (derivYTheta * derivZRho - derivYRho * derivZTheta) &
           + derivXTheta * (derivYRho * derivZZeta - derivYZeta * derivZRho)

  gradRhoX = (derivYZeta * derivZTheta - derivYTheta * derivZZeta) / jacobian
  gradRhoY = (derivZZeta * derivXTheta - derivZTheta * derivXZeta) / jacobian
  gradRhoZ = (derivXZeta * derivYTheta - derivXTheta * derivYZeta) / jacobian

  gradZetaX = (derivYTheta * derivZRho - derivYRho * derivZTheta) / jacobian
  gradZetaY = (derivZTheta * derivXRho - derivZRho * derivXTheta) / jacobian
  gradZetaZ = (derivXTheta * derivYRho - derivXRho * derivYTheta) / jacobian

  gradThetaX = (derivYRho * derivZZeta - derivYZeta * derivZRho) / jacobian
  gradThetaY = (derivZRho * derivXZeta - derivZZeta * derivXRho) / jacobian
  gradThetaZ = (derivXRho * derivYZeta - derivXZeta * derivYRho) / jacobian

  gradRhoSq = gradRhoX**2 + gradRhoY**2 + gradRhoZ**2
  gradRhoGradZeta = gradRhoX * gradZetaX + gradRhoY * gradZetaY + gradRhoZ * gradZetaZ
  gradRhoGradTheta = gradRhoX * gradThetaX + gradRhoY * gradThetaY + gradRhoZ * gradThetaZ

  gradThetaSq = gradThetaX**2 + gradThetaY**2 + gradThetaZ**2
  gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ

  gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2

  ! B field squared is obtained now, then the B field bf; in the following, i, j, k can be 
  ! taken over the full domain because we have all the required quantities everywhere;
  ! thus, extrapolation for Bfield is not necessary anymore

  DO  k = 1,nzeta
     DO  j = 1,npsi
        DO  i = 1,nthe
           bsq(i,j,k) = (gradRhoSq(i,j,k) * gradZetaSq(i,j,k) - gradRhoGradZeta(i,j,k) **2) &
                      * (f(j) * fzet(k)) **2
           bf(i,j,k) = SQRT(bsq(i,j,k))
        END DO
     END DO
  END DO

  ! j dot gradRho
  DO j = 1, npsi
     DO k = 1, nzeta
        jGradRhoPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * &
             gradRhoGradZeta(:,j,k))
        jGradRhoPartialZeta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradZetaSq(:,j,k) - gradRhoGradZeta(:,j,k) **2)
     END DO
  END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialTheta, &
       & derivjGradRhoPartialTheta, derivNU1, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialZeta, &
       & derivNU1, derivNU2, derivjGradRhoPartialZeta)

  jGradRho = (derivjGradRhoPartialTheta + derivjGradRhoPartialZeta)/jacobian

  ! j dot gradZeta
  DO j = 1, npsi
     DO k = 1, nzeta
        jGradZetaPartialRho(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k) **2 - gradRhoSq(:,j,k) * &
             gradZetaSq(:,j,k))
        jGradZetaPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * &
             & gradZetaSq(:,j,k))
     END DO
  END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialRho, &
       & derivNU1, derivjGradZetaPartialRho, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialTheta, &
       & derivjGradZetaPartialTheta, derivNU1, derivNU2)

  jGradZeta = (derivjGradZetaPartialRho + derivjGradZetaPartialTheta)/jacobian

  ! j dot gradTheta
  DO j = 1,npsi
     DO k = 1,nzeta
        jGradThetaPartialRho(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhogradTheta(:,j,k) * gradRhogradZeta(:,j,k) - gradThetagradZeta(:,j,k) * gradRhoSq(:,j,k))
        jGradThetaPartialZeta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhogradTheta(:,j,k) * gradZetaSq(:,j,k) - gradRhogradZeta(:,j,k) * gradThetagradZeta(:,j,k))
     ENDDO
  ENDDO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradThetaPartialRho, &
       & derivNU1, derivjGradThetaPartialRho, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradThetaPartialZeta, &
       & derivNU1, derivNU2, derivjGradThetaPartialZeta)

  jGradTheta = (derivjGradThetaPartialRho + derivjGradThetaPartialZeta)/jacobian

  ! The one below is for calculating \partial Pper/\partial \theta and
  ! \partial(J(Pper-Ppar))/\partial \theta for the anisotropic case  
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jacobian(1:nthe,1:npsi,1:nzeta) * &
       (pper(1:nthe,1:npsi,1:nzeta)-ppar(1:nthe,1:npsi,1:nzeta)), derivDiffPTheta, derivNU1, derivNU2)


  !***************************************************************************************************
  ! Now we have jGradRho, jGradZeta
  ! Time to compute |j x B - div dot P|
  DO k = 1, nzeta
     DO j = 1, npsi
        jCrossBSq(:,j,k) = f(j)**2 * fzet(k)**2 &
                         * (gradRhoSq(:,j,k) * jGradZeta(:,j,k)**2 &
                         + gradZetaSq(:,j,k) * jGradRho(:,j,k)**2 &
                         - 2._dp * jGradZeta(:,j,k) * jGradRho(:,j,k) * gradRhoGradZeta(:,j,k))

        if (isotropy.eq.0) then
           gradPSq(:,j,k) = gradRhoSq(:,j,k)*dPperdRho(:,j,k)**2 &
                          + gradZetaSq(:,j,k)*dPperdZeta(:,j,k)**2 &
                          + gradThetaSq(:,j,k)*dPperdTheta(:,j,k)**2 &
                          + 2.*dPperdRho(:,j,k)*dPperdZeta(:,j,k)*gradRhoGradZeta(:,j,k) &
                          + 2.*dPperdRho(:,j,k)*dPperdTheta(:,j,k)*gradRhoGradTheta(:,j,k) &
                          + 2.*dPperdZeta(:,j,k)*dPperdTheta(:,j,k)*gradThetaGradzeta(:,j,k) &
                          + (derivDiffPTheta(:,j,k)/jacobian(:,j,k))**2 &
                          - 2.*dPperdTheta(:,j,k) * derivDiffPTheta(:,j,k)/jacobian(:,j,k)

           jCrossBMinusGradPSq(:,j,k) = gradRhoSq(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k))**2 + &
                gradZetaSq(:,j,k)*(f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k))**2 + gradThetaSq(:,j,k)*dPperdTheta(:,j,k)**2 + &
                (derivDiffPTheta(:,j,k)/jacobian(:,j,k))**2 - &
                2.*gradRhoGradZeta(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k)) * &
                (f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k)) - &
                2.*gradRhoGradTheta(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k))*dPperdTheta(:,j,k) + &
                2.*gradThetaGradzeta(:,j,k)*(f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k))*dPperdTheta(:,j,k) - &
                2.*dPperdTheta(:,j,k)*derivDiffPTheta(:,j,k)/jacobian(:,j,k)
        else
           gradPSq(:,j,k) = (fzet(k)**2 * dpdAlpha(1:nthe,j,k)**2 * gradZetaSq(:,j,k) + f(j)**2 * &
                dpdPsi(1:nthe,j,k)**2 * gradRhoSq(:,j,k) + 2._dp*dpdAlpha(1:nthe,j,k)*dpdPsi(1:nthe,j,k) * &
                f(j) * fzet(k) * gradRhoGradZeta(:,j,k))

           jCrossBMinusGradPSq(:,j,k) = f(j)**2 * fzet(k)**2 * (gradRhoSq(:,j,k) * (jGradZeta(:,j,k) - &
                1. / fzet(k) * dpdPsi(1:nthe,j,k))**2 + &
                gradZetaSq(:,j,k) * (jGradRho(:,j,k) + 1./f(j) * dpdAlpha(1:nthe,j,k))**2 - 2._dp * &
                (jGradZeta(:,j,k) - 1./fzet(k) * dpdPsi(1:nthe,j,k)) * (jGradRho(:,j,k) + 1./f(j)* &
                dpdAlpha(1:nthe,j,k)) * gradRhoGradZeta(:,j,k))
        endif

     END DO
  END DO

  distance = sqrt(x**2+y**2+z**2)
  jCrossB = SQRT(jCrossBSq)*(bnormal*pjconst/6.4)*distance(:,:,1:nzeta)
  gradP = SQRT(ABS(gradPSq))*(pnormal*2.0)/distance(:,:,1:nzeta)
  jCrossBMinusGradPMod = SQRT(ABS(jCrossBMinusGradPSq))

  ! Force balance quantities
  FileName = trim(prefixOut)//'Force_balance_equatorial'
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) npsi, nzeta
  DO j = 1, npsi
     DO k = 1, nzeta
        WRITE(UNITTMP_, *) x(nThetaEquator, j, k), y(nThetaEquator, j, k), bf(nThetaEquator, j, k)*bnormal, &
             jCrossB(nThetaEquator,j,k), jCrossBMinusGradPMod(nThetaEquator,j,k), &
             gradP(nThetaEquator,j,k), jCrossB(nThetaEquator,j,k)/gradP(nThetaEquator,j,k)
     END DO
  END DO
  CLOSE(UNITTMP_)

  FileName = trim(prefixOut)//'Force_balance_midnight'
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) npsi, nthe
  DO j = 1, npsi
     DO i = 1, nthe
        WRITE(UNITTMP_, *) x(i,j,nZetaMidnight), z(i,j,nZetaMidnight), &
             bf(i,j,nZetaMidnight)*bnormal, jCrossB(i,j,nZetaMidnight), &
             jCrossBMinusGradPMod(i,j,nZetaMidnight), gradP(i,j,nZetaMidnight), &
             jCrossBMinusGradPMod(i,j,nZetaMidnight)/jCrossB(i,j,nZetaMidnight)
     END DO
  END DO
  CLOSE(UNITTMP_)

  call metric
  call newj
  DO j = 2, npsi-1
     i = nThetaEquator
     k = nZetaMidnight
     rtemp(j) = - vecd(i,j,k)*psi(i,j,k) &
                + vec1(i,j,k)*psi(i-1,j-1,k) &
                + vec2(i,j,k)*psi(i,j-1,k) &
                + vec3(i,j,k)*psi(i+1,j-1,k) &
                + vec4(i,j,k)*psi(i-1,j,k) &
                + vec6(i,j,k)*psi(i+1,j,k) &
                + vec7(i,j,k)*psi(i-1,j+1,k) &
                + vec8(i,j,k)*psi(i,j+1,k) &
                + vec9(i,j,k)*psi(i+1,j+1,k) &
                - vecr(i,j,k)
  ENDDO
  rtemp(1) = 0.
  rtemp(npsi) = 0.

  call metrica
  call newk
  DO j = 1, npsi
     i = nThetaEquator
     k = nZetaMidnight
     xtemp(j) = - vecd(i,j,k)*alfa(i,j,k)  &
                + vec1(i,j,k)*alfa(i-1,j,k-1) &
                + vec2(i,j,k)*alfa(i,j,k-1) &
                + vec3(i,j,k)*alfa(i+1,j,k-1) &
                + vec4(i,j,k)*alfa(i-1,j,k)  &
                + vec6(i,j,k)*alfa(i+1,j,k)  &
                + vec7(i,j,k)*alfa(i-1,j,k+1) &
                + vec8(i,j,k)*alfa(i,j,k+1) &
                + vec9(i,j,k)*alfa(i+1,j,k+1) &
               - vecx(i,j,k)
  ENDDO

  FileName = trim(prefixOut)//'Force_balance_line'
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace') 
  WRITE(UNITTMP_, *) npsi
  DO j = 1, npsi
     WRITE(UNITTMP_, *) x(nThetaEquator,j,nZetaMidnight), y(nThetaEquator,j,nZetaMidnight), &
          jCrossB(nThetaEquator,j,nZetaMidnight), gradP(nThetaEquator,j,nZetaMidnight), &
          jCrossBMinusGradPMod(nThetaEquator,j,nZetaMidnight), rtemp(j), xtemp(j), &
          jacobian(nThetaEquator,j,nZetaMidnight), f(j), fzet(nZetaMidnight)
  END DO
  CLOSE(UNITTMP_)

  normDiff = 0.0_dp
  normDiffRel = 0.0_dp
  normJxB  = 0.0_dp
  normGradP = 0.0_dp
  volume   = 0.0_dp
  DO i = 2, nthe-1
     DO j = 2, npsi-1
        DO k = 2, nzeta
           !IF (2.*pper(i,j,k) > 1.E-1_dp*bsq(i,j,k)) THEN
              ! Only in regions with beta > 1E-2 
              ! (in regions of low plasma beta, the pressure does not change the magnetic field)
              normDiff = normDiff + jacobian(i,j,k) * dr * dpPrime * dt * jCrossBMinusGradPMod(i,j,k)
              normDiffRel = normDiffRel + jacobian(i,j,k) * dr * dpPrime * dt * jCrossBMinusGradPMod(i,j,k) / jCrossB(i,j,k)
              normJxB = normJxB + jacobian(i,j,k) * dr * dpPrime * dt * jCrossB(i,j,k)
              normGradP = normGradP + jacobian(i,j,k) * dr * dpPrime * dt * gradP(i,j,k)
              volume = volume + jacobian(i,j,k) * dr * dpPrime * dt
           !END IF
        END DO
     END DO
  END DO

  ! Normalize to total computational volume
  normDiff = normDiff/volume
  normJxB = normJxB/volume
  normGradP = normGradP/volume

  !  Norms of |jxB-grad P|,      |jxB|,      |gradP| 
  WRITE(*, *) normDiff, normJxB, normGradP

  RETURN

  END SUBROUTINE Write_convergence_anisotropic

END MODULE ModScbIO
