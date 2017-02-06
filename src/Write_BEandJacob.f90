SUBROUTINE Write_BEandJacob

  !****************************************************************************
  ! New Netcdf style (time included as an unlimited dimension)
  ! Author: Sorin Zaharia, Nov. 2009
  ! Electric field stuff not done yet
  !
  ! Copyright (c) 2016, Los Alamos National Security, LLC
  ! All rights reserved.
  !****************************************************************************

  use ModNumConst, ONLY: cRadToDeg
  USE Module1
  USE NETCDF
  use ModRamMain, ONLY: DoWriteB
  IMPLICIT NONE

  INTEGER :: i, j, k, ncid
  INTEGER :: chiid, alphaid, betaid, timeid, chivarid, alphavarid, betavarid, &
       xid, yid, zid, bxid, byid, bzid, jacid, tiltid, timevarid
  INTEGER, SAVE :: iCALL = 0
  INTEGER :: START(1), COUNT(1) ! For time
  INTEGER :: START1D(1)
  INTEGER :: START3D(4), COUNT3D(4) ! For 3-D arrays
  INTEGER :: START2D(2), COUNT2D(2) ! For 1-D arrays
  CHARACTER*500 :: filenameOut1, filenameOut2, filenameOut3

  REAL :: time(1), tilt_local(1)
  REAL(DP) :: VGSEX = -400._dp, VGSEY=0._dp, VGSEZ=0._dp ! Solar wind speed
  REAL(DP) :: AAA(10), SPS, CPS, BBB(22)
  REAL(DP) :: G_GEOPACK(105), H_GEOPACK(105), REC(105)
  COMMON /GEOPACK1/ AAA,SPS,CPS,BBB
  COMMON /GEOPACK2/ G_GEOPACK,H_GEOPACK,REC

  time(1) = rHour

  START = (/iCALL+1/)
  COUNT = (/1/)

  START1D = (/iCALL+1/)

  START2D = (/1,iCALL+1/)

  START3D = (/1,1,1,iCALL+1/)
  COUNT3D = (/nthe,npsi,nzeta,1/)

  ! Subtract internal field
  CALL RECALC_08(2005,243+iDay,iHour,iMin,0,VGSEX,VGSEY,VGSEZ) ! Just an example (for Sep. 2005 storm, DOY=243 for Aug.31,2005)
  tilt_local(1) = REAL(cRadToDeg*acos(CPS))
 
 ! Change tilt back to zero (bxintern, byintern, bzintern are WANTED in SM)
  SPS = 0._dp
  CPS = 1._dp
  ! tilt_local(1) = 0.  ! Uncomment this if you want zero tilt

  ! Subtract internal field (dipolar or IGRF)
  DO k = 1, nzeta
     DO j = 1, npsi
        DO i = 1, nthe
           CALL DIP_08(x(i,j,k),y(i,j,k),z(i,j,k),bxintern(i,j,k),byintern(i,j,k),bzintern(i,j,k)) ! Dipolar
      !     bX(i,j,k) = bX(i,j,k) - bxintern/bnormal
      !     bY(i,j,k) = bY(i,j,k) - byintern/bnormal
      !     bZ(i,j,k) = bZ(i,j,k) - bzintern/bnormal
        END DO
     END DO
  END DO

  ! Normalize to code coordinates
  bxintern = bxintern/bnormal
  byintern = byintern/bnormal
  bzintern = bzintern/bnormal

  IF (DoWriteB.eqv..false.) THEN   !VJ  do not write output
   RETURN
  END IF

  filenameOut1 = TRIM(ADJUSTL(prefixOut))//'/mag_field.nc'
  IF (iInducedE == 1 .OR. iConvE == 1) THEN
     filenameOut2 = TRIM(ADJUSTL(prefixOut))//'/e_field.nc'
     filenameOut3 = TRIM(ADJUSTL(prefixOut))//'/grid_points.nc'
  ENDIF


  First_time_call : IF(iCALL == 0) THEN
     CALL check (nf90_create(filenameOut1, nf90_clobber, ncid))

     ! Define dimensions

     CALL check(nf90_def_dim(ncid, 'chi', nthe, chiid))
     CALL check(nf90_def_dim(ncid, 'alpha', npsi, alphaid))
     CALL check(nf90_def_dim(ncid, 'beta', nzeta, betaid))
     CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, timeid))

     ! Define variables
     CALL check(nf90_def_var(ncid, 'chi', nf90_float,(/chiid,timeid/),chivarid))
     CALL check(nf90_put_att(ncid,chivarid,'title','Coordinate along field line'))
     CALL check(nf90_def_var(ncid, 'alpha', nf90_float,(/alphaid,timeid/),alphavarid))
     CALL check(nf90_put_att(ncid,alphavarid,'title','Magnetic flux-like Euler potential'))
     CALL check(nf90_def_var(ncid, 'beta', nf90_float,(/betaid,timeid/),betavarid))
     CALL check(nf90_put_att(ncid,betavarid,'title','Azimuthal angle-like Euler potential'))
     CALL check(nf90_def_var(ncid, 'time', nf90_float,timeid,timevarid))
     CALL check(nf90_put_att(ncid,timevarid,'title','Time'))

     CALL check(nf90_def_var(ncid, 'x', nf90_float, (/chiid,alphaid,betaid,timeid/),xid))
     CALL check(nf90_put_att(ncid,xid,'title','3D array of x locations (GSM)'))
     CALL check(nf90_def_var(ncid, 'y', nf90_float, (/chiid,alphaid,betaid,timeid/),yid))
     CALL check(nf90_put_att(ncid,yid,'title','3D array of y locations (GSM)'))
     CALL check(nf90_def_var(ncid, 'z', nf90_float, (/chiid,alphaid,betaid,timeid/),zid))
     CALL check(nf90_put_att(ncid,zid,'title','3D array of z locations (GSM)'))
     CALL check(nf90_def_var(ncid, 'bx', nf90_float, (/chiid,alphaid,betaid,timeid/),bxid))
     CALL check(nf90_put_att(ncid,bxid,'title','3D array of bx (external) values (GSM)'))
     CALL check(nf90_def_var(ncid, 'by', nf90_float, (/chiid,alphaid,betaid,timeid/),byid))
     CALL check(nf90_put_att(ncid,byid,'title','3D array of by (external) values (GSM)'))
     CALL check(nf90_def_var(ncid, 'bz', nf90_float, (/chiid,alphaid,betaid,timeid/),bzid))
     CALL check(nf90_put_att(ncid,bzid,'title','3D array of bz (external) values (GSM)'))
     CALL check(nf90_def_var(ncid, 'jacobian', nf90_float, (/chiid,alphaid,betaid,timeid/),jacid))
     CALL check(nf90_put_att(ncid,jacid,'title','3D array of jacobian values'))

     CALL check(nf90_def_var(ncid, 'tilt', nf90_float, timeid,tiltid))
     CALL check(nf90_put_att(ncid,tiltid,'title', 'Tilt of dipole axis (degrees)'))

     ! End define mode
     CALL check(nf90_enddef(ncid))


  ELSE ! Open existing NetCDF file
     CALL check(nf90_open(filenameOut1, nf90_write, ncid))

     CALL check( nf90_inq_dimid(ncid, 'chi', chiid))
     CALL check( nf90_inq_dimid(ncid, 'alpha', alphaid))
     CALL check( nf90_inq_dimid(ncid, 'beta', betaid))
     CALL check( nf90_inq_dimid(ncid, 'time', timeid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'chi', chivarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'alpha',alphavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'beta', betavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',timevarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'tilt',tiltid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'x', xid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'y', yid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'z', zid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'bx', bxid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'by', byid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'bz', bzid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jacobian', jacid))

  END IF First_time_call

  ! Write mode - write at all times
  CALL check(nf90_put_var(ncid, chivarid, REAL(chiVal(1:nthe)), START2D))
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START2D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START2D))

  CALL check(nf90_put_var(ncid, tiltid, REAL(tilt_local,sp), START, COUNT))

  CALL check(nf90_put_var(ncid, xid, REAL(x(1:nthe,1:npsi,1:nzeta)*COS(tilt_local(1)*pi_d/180.) + &
       z(1:nthe,1:npsi,1:nzeta)*SIN(tilt_local(1)*pi_d/180.)), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, yid, REAL(y(1:nthe,1:npsi,1:nzeta)), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, zid, REAL(-x(1:nthe,1:npsi,1:nzeta)*SIN(tilt_local(1)*pi_d/180.) + &
       z(1:nthe,1:npsi,1:nzeta)*COS(tilt_local(1)*pi_d/180.)), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, bxid, REAL(bnormal*((bx(1:nthe,1:npsi,1:nzeta)-bxintern(1:nthe,1:npsi,1:nzeta))*&
       COS(tilt_local(1)*pi_d/180.) + &
       (bz(1:nthe,1:npsi,1:nzeta)-bzintern(1:nthe,1:npsi,1:nzeta))*SIN(tilt_local(1)*pi_d/180.))), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, byid, REAL(bnormal*(by(1:nthe,1:npsi,1:nzeta)-byintern(1:nthe,1:npsi,1:nzeta))), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, bzid, REAL(bnormal*((-bx(1:nthe,1:npsi,1:nzeta)+bxintern(1:nthe,1:npsi,1:nzeta))*&
       SIN(tilt_local(1)*pi_d/180.) + &
       (bz(1:nthe,1:npsi,1:nzeta)-bzintern(1:nthe,1:npsi,1:nzeta))*COS(tilt_local(1)*pi_d/180.))), START3D,COUNT3D))
  CALL check(nf90_put_var(ncid, jacid, REAL(jacobian(1:nthe,1:npsi,1:nzeta)), START3D,COUNT3D))

  CALL check(nf90_put_var(ncid, timevarid, time, START, COUNT))

  CALL check(nf90_close(ncid))

  !C  PRINT*, 'iCALL, tilt = ', iCALL, tilt_local(1)

  iCALL = iCALL+1

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

END SUBROUTINE Write_BEandJacob
