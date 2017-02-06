SUBROUTINE Write_pressure_anisotropic

  USE Module1
  USE NETCDF

  IMPLICIT NONE

  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  INTEGER :: ierr, j, k, ncid

  integer :: START(1), COUNT(1)
  integer :: START1D(2), COUNT1D(2)
  INTEGER :: START2D(3), COUNT2D(3) ! For 2-D arrays (+ time)

  INTEGER :: chiid, alphaid, betaid, timeid, chivarid, alphavarid, betavarid, &
       xeqid, yeqid, ppereqid, ppareqid, betaeqid, taueqid, xmidid, zmidid, ppermidid, pparmidid, &
       betamidid, taumidid, xnooid, znooid, ppernooid, pparnooid, betanooid, taunooid, timevarid

  INTEGER, SAVE :: iCALLWP = 0

  REAL :: time(1)

  CHARACTER*500 :: fileName

  time(1) = rHour

  START = (/iCALLWP+1/)
  COUNT = (/1/)

  START1D = (/1,iCALLWP+1/)
  START2D = (/1,1,iCALLWP+1/)


  fileName = trim(adjustl(prefixOut))//'pressure_anisotropic.nc'

  First_time_call : IF(iCALLWP == 0) THEN
     CALL check (nf90_create(filename, nf90_clobber, ncid))

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

     CALL check(nf90_def_var(ncid, 'xEq', nf90_float, (/alphaid,betaid,timeid/),xeqid))
     CALL check(nf90_put_att(ncid,xeqid,'title','2D array of xEq locations '))
     CALL check(nf90_def_var(ncid, 'yEq', nf90_float, (/alphaid,betaid,timeid/),yeqid))
     CALL check(nf90_put_att(ncid,yeqid,'title','2D array of yEq locations '))
     CALL check(nf90_def_var(ncid, 'presPerpEq', nf90_float, (/alphaid,betaid,timeid/),ppereqid))
     CALL check(nf90_put_att(ncid,ppereqid,'title','Perpendicular pressure in the equatorial plane (nPa)'))
     CALL check(nf90_def_var(ncid, 'presParEq', nf90_float, (/alphaid,betaid,timeid/),ppareqid))
     CALL check(nf90_put_att(ncid,ppareqid,'title','Parallel pressure in the equatorial plane (nPa)'))
     CALL check(nf90_def_var(ncid, 'betaPerpEq', nf90_float, (/alphaid,betaid,timeid/),betaeqid))
     CALL check(nf90_put_att(ncid,betaeqid,'title','Perpendicular plasma beta'))
     CALL check(nf90_def_var(ncid, 'tauEq', nf90_float, (/alphaid,betaid,timeid/),taueqid))
     CALL check(nf90_put_att(ncid,taueqid,'title','Tau stability criterion'))

     CALL check(nf90_def_var(ncid, 'xMid', nf90_float, (/chiid,alphaid,timeid/),xmidid))
     CALL check(nf90_put_att(ncid,xmidid,'title','2D array of xMid locations '))
     CALL check(nf90_def_var(ncid, 'zMid', nf90_float, (/chiid,alphaid,timeid/),zmidid))
     CALL check(nf90_put_att(ncid,zmidid,'title','2D array of zmid locations '))
     CALL check(nf90_def_var(ncid, 'presPerpMid', nf90_float, (/chiid,alphaid,timeid/),ppermidid))
     CALL check(nf90_put_att(ncid,ppermidid,'title','Perpendicular pressure in midnight meridian (nPa)'))
     CALL check(nf90_def_var(ncid, 'presParMid', nf90_float, (/chiid,alphaid,timeid/),pparmidid))
     CALL check(nf90_put_att(ncid,pparmidid,'title','Parallel pressure in midnight meridian (nPa)'))
     CALL check(nf90_def_var(ncid, 'betaPerpMid', nf90_float, (/chiid,alphaid,timeid/),betamidid))
     CALL check(nf90_put_att(ncid,betamidid,'title','Perpendicular mid plasma beta'))
     CALL check(nf90_def_var(ncid, 'tauMid', nf90_float, (/chiid,alphaid,timeid/),taumidid))
     CALL check(nf90_put_att(ncid,taumidid,'title','Tau mid stability criterion'))

     CALL check(nf90_def_var(ncid, 'xNoo', nf90_float, (/chiid,alphaid,timeid/),xnooid))
     CALL check(nf90_put_att(ncid,xnooid,'title','2D array of xNoo locations '))
     CALL check(nf90_def_var(ncid, 'zNoo', nf90_float, (/chiid,alphaid,timeid/),znooid))
     CALL check(nf90_put_att(ncid,znooid,'title','2D array of znoo locations '))
     CALL check(nf90_def_var(ncid, 'presPerpNoo', nf90_float, (/chiid,alphaid,timeid/),ppernooid))
     CALL check(nf90_put_att(ncid,ppernooid,'title','Perpendicular pressure in noon meridian (nPa)'))
     CALL check(nf90_def_var(ncid, 'presParNoo', nf90_float, (/chiid,alphaid,timeid/),pparnooid))
     CALL check(nf90_put_att(ncid,pparnooid,'title','Parallel pressure in noon meridian (nPa)'))
     CALL check(nf90_def_var(ncid, 'betaPerpNoo', nf90_float, (/chiid,alphaid,timeid/),betanooid))
     CALL check(nf90_put_att(ncid,betanooid,'title','Perpendicular noon plasma beta'))
     CALL check(nf90_def_var(ncid, 'tauNoo', nf90_float, (/chiid,alphaid,timeid/),taunooid))
     CALL check(nf90_put_att(ncid,taunooid,'title','Tau noon stability criterion'))


     ! End define mode
     CALL check(nf90_enddef(ncid))


  ELSE ! Open existing NetCDF file
     CALL check(nf90_open(filename, nf90_write, ncid))

     CALL check( nf90_inq_dimid(ncid, 'chi', chiid))
     CALL check( nf90_inq_dimid(ncid, 'alpha', alphaid))
     CALL check( nf90_inq_dimid(ncid, 'beta', betaid))
     CALL check( nf90_inq_dimid(ncid, 'time', timeid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'chi', chivarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'alpha',alphavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'beta', betavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',timevarid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xEq', xeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yEq', yeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presPerpEq', ppereqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presParEq', ppareqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'betaPerpEq', betaeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'tauEq', taueqid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xMid', xmidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'zMid', zmidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presPerpMid', ppermidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presParMid', pparmidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'betaPerpMid', betamidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'tauMid', taumidid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xNoo', xnooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'zNoo', znooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presPerpNoo', ppernooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'presParNoo', pparnooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'betaPerpNoo', betanooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'tauNoo', taunooid))

  END IF First_time_call

  ! Write mode - write at all times
  CALL check(nf90_put_var(ncid, chivarid, REAL(chiVal(1:nthe)), START1D))
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START1D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START1D))

  CALL check(nf90_put_var(ncid, xeqid, REAL(x(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yeqid, REAL(y(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, ppereqid, REAL(pper(nThetaEquator,:,1:nzeta)*pnormal),START2D))
  CALL check(nf90_put_var(ncid, ppareqid, REAL(ppar(nThetaEquator,:,1:nzeta)*pnormal),START2D)) 
  CALL check(nf90_put_var(ncid, betaeqid, REAL(2.*pper(nThetaEquator,:,1:nzeta)/bsq(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, taueqid, REAL(tau(nThetaEquator,:,1:nzeta)),START2D))

  CALL check(nf90_put_var(ncid, xmidid, REAL(x(:,:,nZetaMidnight)),START2D))
  CALL check(nf90_put_var(ncid, zmidid, REAL(z(:,:,nZetaMidnight)),START2D))
  CALL check(nf90_put_var(ncid, ppermidid, REAL(pper(:,:,nZetaMidnight)*pnormal),START2D))
  CALL check(nf90_put_var(ncid, pparmidid, REAL(ppar(:,:,nZetaMidnight)*pnormal),START2D)) 
  CALL check(nf90_put_var(ncid, betamidid, REAL(2.*pper(:,:,nZetaMidnight)/bsq(:,:,nZetaMidnight)),START2D))
  CALL check(nf90_put_var(ncid, taumidid, REAL(tau(:,:,nZetaMidnight)),START2D))

  CALL check(nf90_put_var(ncid, xnooid, REAL(x(:,:,2)),START2D))
  CALL check(nf90_put_var(ncid, znooid, REAL(z(:,:,2)),START2D))
  CALL check(nf90_put_var(ncid, ppernooid, REAL(pper(:,:,2)*pnormal),START2D))
  CALL check(nf90_put_var(ncid, pparnooid, REAL(ppar(:,:,2)*pnormal),START2D))
  CALL check(nf90_put_var(ncid, betanooid, REAL(2.*pper(:,:,2)/bsq(:,:,2)),START2D))
  CALL check(nf90_put_var(ncid, taunooid, REAL(tau(:,:,2)),START2D))

  CALL check(nf90_put_var(ncid, timevarid, time, START, COUNT))

  CALL check(nf90_close(ncid))

  iCALLWP = iCALLWP+1

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


END SUBROUTINE Write_pressure_anisotropic
