!******************************************************************************
SUBROUTINE plot_physical     
!    computes and plots physical quantities
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

  USE nrtype
  USE Module1
  USE NETCDF
  use ModRamMain, ONLY: DoWriteCur
  IMPLICIT NONE

  INTERFACE splint_interface
     FUNCTION splint(xa,ya,y2a,x)
       USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
       USE nr, ONLY: locate
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: splint
       INTEGER(I4B) :: khi,klo,n
       REAL(DP) :: a,b,h
     END FUNCTION splint
  END INTERFACE

  INTERFACE spline_interface
     SUBROUTINE spline(x,y,yp1,ypn,y2)
       USE nrtype; USE nrutil, ONLY : assert_eq
       USE nr, ONLY : tridag
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(DP), INTENT(IN) :: yp1,ypn
       REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
       INTEGER(I4B) :: n
       REAL(DP), DIMENSION(SIZE(x)) :: a,b,c,r
     END SUBROUTINE spline
  END INTERFACE

  REAL(DP) :: rrx(npsi,nzeta+1),rry(npsi,nzeta+1),paraj(npsi,nzeta+1), parajDirect(npsi,nzeta+1)
  REAL(DP) :: rrlatx(201),rrlaty(201), wlabels(4)

  INTEGER :: i, iTheta, j, k, iplotx, idx, ixbeg, ixend, kkx, nz1, nz2, nz3, niz, niy, iz, iy, &
       & nzang, kz, k1, k2, nwrite, nearth, jmin, jmax, iplx, kndps, ndps, npsid, &
       & jj, ndpsi, imax, ierr, ncdfid, idealerr, ncid

 INTEGER :: chiid, alphaid, betaid, timeid, chivarid, alphavarid, betavarid, &
       xeqid, yeqid, jcrosseqid, ppareqid, betaeqid, xpoleid, ypoleid, jparid, jpardirectid, xmidid, &
       zmidid, jcrossmidid, pparmidid, &
       betamidid, taumidid, xnooid, znooid, jcrossnooid, pparnooid, betanooid, taunooid, timevarid

 integer :: START(1), COUNT(1)
 integer :: START1D(2), COUNT1D(2)
 INTEGER :: START2D(3), COUNT2D(3) ! For 2-D arrays (+ time)
 
  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  REAL(DP) :: rr1, rr2, zangle, thangle, thangleOnEarth, rr, dza, dya, &
       & za, z1, ya, x1, y1, xmin3, &
       & xmax3, zaxis, dz3, dxz3, xmax2, xmin2, zmax2, zmin2, zmax3, zmin3, curmin, curmax,  &
       & dzang, zang, xtot, ztot, xw, zw, xwrt, ywrt, thet, xxs, zzs, pres, zmax, zmin, dx, dz, &
       & bd, at, yyp
  REAL(DP), PARAMETER :: oneThousandKm = 1000./6371.
  REAL(DP) :: distance(nthe/2 + 1), distance2derivs(nthe/2+1), parCurrent(nthe/2+1), dist, &
       dipoleFactor, dipoleFactor4RE, factorIncrease
  REAL(DP), ALLOCATABLE :: xtemp(:), delb(:), ytemp(:), wtemp(:), ztemp(:)
  REAL :: aMatrix(2,4)
  CHARACTER*100 :: fileName

  INTEGER, SAVE :: iCALLPP = 0

  REAL :: time(1)

  time(1) = rHour
  
  START = (/iCALLPP+1/)
  COUNT = (/1/)

  START1D = (/1,iCALLPP+1/)
  START2D = (/1,1,iCALLPP+1/)

  !ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !ccc  write position and field-aligned current in the polar region
  !ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !cccc  prepare arrays for contour plots on the polar region
  !cccc  paraj is the field-aligned current at northern polar cap, obtained by integration

  !  pio2_d is pi/2, defined in nrtype (DP).

  i = nthe  ! At the north pole

  alphaLoop :  DO k = 2, nzetp
     psiLoop :   DO j = 1, npsi
        rr2 = xx(i,j,k)**2 + z(i,j,k)**2
        rr1 = SQRT(rr2)
        zangle = yy(i,j,k)
        !!  thangle is the polar angle
        thangle = ASIN(z(i,j,k) / rr1)
        IF ((ABS(thangle) > 1000.) .OR. (ABS(zangle) > 1000.)) THEN
           PRINT*, j, k, rr1, z(i,j,k), thangle, zangle
           STOP
        END IF

        thangleOnEarth = ACOS(SQRT((COS(thangle))**2/r0Start))

        !!  radius that corresponds to polar angle is normalized to 1 for every 10 degree in 
        ! latitude zangle is toroidal angle; zangle = 0 is local noon
        IF (INT(r0Start) /= 1) rr = (90._dp - thangleOnEarth * 360._dp / twopi_d) / 10.0_dp
        IF (INT(r0Start) == 1) rr = (90._dp - thangle * 360._dp / twopi_d) / 10.0_dp

        rrx(j,k) = rr * COS(zangle + pio2_d)
        rry(j,k) = rr * SIN(zangle + pio2_d)

        !C if (j==1) print*, 'plot_physical: j, k, zangle, rr, rrx, rry = ', j, k, zangle, rr, rrx(j,k), rry(j,k)

        ! Plot at 1000 km above Earth's surface using interpolation
        !        distance(nThetaEquator) = 0._dp
        !        parCurrent(nThetaEquator) = 0._dp
        !        DO iTheta = nThetaEquator-1, 1, -1
        !           distance(iTheta) = distance(iTheta+1) + SQRT((x(iTheta,j,k)-x(iTheta+1,j,k))**2 + &
        !                (y(iTheta,j,k)-y(iTheta+1,j,k))**2 + (z(iTheta,j,k)-z(iTheta+1,j,k))**2)
        !           parCurrent(iTheta) = bj(nthe - iTheta + 1, j, k)
        !     END DO
        !CALL spline(distance(1:nThetaEquator), parCurrent(1:nThetaEquator), 1.E31_dp, 1.E31_dp, &
        !     distance2derivs(1:nThetaEquator))
        !paraj(j,k) = splint(distance(1:nThetaEquator), parCurrent(1:nThetaEquator), &
        !     distance2derivs(1:nThetaEquator), oneThousandKm)

        !        dist =  SQRT((x(nthe,j,k)-x(nthe-1,j,k))**2 + &
        !     (y(nthe,j,k)-y(nthe-1,j,k))**2 + (z(nthe,j,k)-z(nthe-1,j,k))**2)
        !        paraj(j,k) = bj(nthe,j,k) - oneThousandKm * (bj(nthe,j,k) - bj(nthe-1,j,k)) / dist

        ! Dipole field at the Earth's surface
        dipoleFactor = SQRT(1. + 3. * (SIN(thangleOnEarth))**2) 
        dipoleFactor4RE = SQRT(1. + 3. * (SIN(thangle))**2)
        factorIncrease = dipoleFactor * r0Start**3 / dipoleFactor4RE

        IF (INT(r0Start) /= 1) THEN
           paraj(j,k) = bj(1,j,k) * factorIncrease
           parajDirect(j,k) = jParDirect(1,j,k) * factorIncrease
        ELSE
           paraj(j,k) = bj(1,j,k)
           parajDirect(j,k) = jParDirect(1,j,k)
        END IF
     END DO psiLoop
  END DO alphaLoop

  
  fileName = TRIM(ADJUSTL(prefixOut))//'currents.nc'

  IF (DoWriteCur.eqv..false.) THEN   !VJ  do not write output
   RETURN
  END IF
  
  First_time_call : IF(iCALLPP == 0) THEN
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
     CALL check(nf90_def_var(ncid, 'jCrossEq', nf90_float, (/alphaid,betaid,timeid/),jcrosseqid))
     CALL check(nf90_put_att(ncid,jcrosseqid,'title','Perpendicular current density in the equatorial plane (muA/m^2)'))
    
     CALL check(nf90_def_var(ncid, 'xPole', nf90_float, (/alphaid,betaid,timeid/),xpoleid))
     CALL check(nf90_put_att(ncid,xpoleid,'title','2D array of xEq locations '))
     CALL check(nf90_def_var(ncid, 'yPole', nf90_float, (/alphaid,betaid,timeid/),ypoleid))
     CALL check(nf90_put_att(ncid,ypoleid,'title','2D array of yEq locations '))
     CALL check(nf90_def_var(ncid, 'jParVas', nf90_float, (/alphaid,betaid,timeid/),jparid))
     CALL check(nf90_put_att(ncid,jparid,'title','Parallel current density at ionosphere (Vasyliunas eq.) (muA/m^2)'))
     CALL check(nf90_def_var(ncid, 'jParAmpere', nf90_float, (/alphaid,betaid,timeid/),jpardirectid))
     CALL check(nf90_put_att(ncid,jpardirectid,'title','Parallel current density at ionosphere (Ampere eq.) (muA/m^2)'))

     CALL check(nf90_def_var(ncid, 'xMid', nf90_float, (/chiid,alphaid,timeid/),xmidid))
     CALL check(nf90_put_att(ncid,xmidid,'title','2D array of xMid locations '))
     CALL check(nf90_def_var(ncid, 'zMid', nf90_float, (/chiid,alphaid,timeid/),zmidid))
     CALL check(nf90_put_att(ncid,zmidid,'title','2D array of zmid locations '))
     CALL check(nf90_def_var(ncid, 'jCrossMid', nf90_float, (/chiid,alphaid,timeid/),jcrossmidid))
     CALL check(nf90_put_att(ncid,jcrossmidid,'title','Perpendicular current density in midnight meridian (muA/m^2)'))
   
     CALL check(nf90_def_var(ncid, 'xNoo', nf90_float, (/chiid,alphaid,timeid/),xnooid))
     CALL check(nf90_put_att(ncid,xnooid,'title','2D array of xNoo locations '))
     CALL check(nf90_def_var(ncid, 'zNoo', nf90_float, (/chiid,alphaid,timeid/),znooid))
     CALL check(nf90_put_att(ncid,znooid,'title','2D array of znoo locations '))
     CALL check(nf90_def_var(ncid, 'jCrossNoo', nf90_float, (/chiid,alphaid,timeid/),jcrossnooid))
     CALL check(nf90_put_att(ncid,jcrossnooid,'title','Perpendicular current density in noon meridian (muA/m^2)'))
    
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
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jCrossEq', jcrosseqid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xPole', xpoleid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yPole', ypoleid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jParVas', jparid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jParAmpere', jpardirectid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xMid', xmidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'zMid', zmidid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jCrossMid', jcrossmidid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xNoo', xnooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'zNoo', znooid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'jCrossNoo', jcrossnooid))

  END IF First_time_call


 ! Write mode - write at all times

  CALL check(nf90_put_var(ncid, chivarid, REAL(chiVal(1:nthe)), START1D))
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START1D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START1D))

  CALL check(nf90_put_var(ncid, xeqid, REAL(x(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yeqid, REAL(y(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, jcrosseqid, REAL(phij(nThetaEquator,:,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, xpoleid, REAL(rrx(1:npsi,2:nzetp)),START2D)) 
  CALL check(nf90_put_var(ncid, ypoleid, REAL(rry(1:npsi,2:nzetp)),START2D))
  CALL check(nf90_put_var(ncid, jparid, REAL(paraj(1:npsi,2:nzetp)),START2D))
  CALL check(nf90_put_var(ncid, jpardirectid, REAL(parajDirect(1:npsi,2:nzetp)),START2D))

  CALL check(nf90_put_var(ncid, xmidid, REAL(x(:,:,nZetaMidnight)),START2D))
  CALL check(nf90_put_var(ncid, zmidid, REAL(z(:,:,nZetaMidnight)),START2D))
  CALL check(nf90_put_var(ncid, jcrossmidid, REAL(phij(:,:,nZetaMidnight)),START2D))
  
  CALL check(nf90_put_var(ncid, xnooid, REAL(x(:,:,2)),START2D))
  CALL check(nf90_put_var(ncid, znooid, REAL(z(:,:,2)),START2D))
  CALL check(nf90_put_var(ncid, jcrossnooid, REAL(phij(:,:,2)),START2D))
 
  CALL check(nf90_put_var(ncid, timevarid, time, START, COUNT))

  CALL check(nf90_close(ncid))

  iCALLPP = iCALLPP+1

  ! Writes full B fields and Jacobian to disk, at the end of computation
  IF (iConvGlobal == 1) CALL Write_BEandJacob 
  
  
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



   END SUBROUTINE plot_physical



