!******************************************************************************
!   Dumps plasma density in NetCdf file (time included as unlimited dimension)
!   Chris Jeffery, Feb. 2015
!   Copyright (c) 2016, Los Alamos National Security, LLC
!   All rights reserved.
!******************************************************************************

SUBROUTINE netcdf_density_dump(density_in, NR, NT, time_in)

  USE NETCDF
  USE ModRamMain, ONLY : PathRamOut, Real8_, Real4_
  USE ModPlane, ONLY : L,mlt

  IMPLICIT NONE

  INTEGER :: NCID, STATUS, RADID, MLTID, DENSITYID, DUMMYID, TIMEID, RADVARID, MLTVARID, TIMEVARID
  REAL(kind=Real8_), INTENT(IN) :: density_in(NR, NT)
  REAL(kind=Real8_), INTENT(IN) :: time_in
  REAL(kind=Real8_) :: time(1) 
  REAL(kind=Real4_) :: density(NR,NT)

  INTEGER :: START3(3), COUNT3(3) ! netCDF variable start point
  INTEGER :: START1(1), COUNT1(1)
  INTEGER, INTENT(IN) :: NR, NT
  INTEGER :: ier

  CHARACTER(LEN=200) :: fileNetcdf

  INTEGER, SAVE :: iCALL = 0

  START3 = (/1, 1, iCALL+1/)
  COUNT3 = (/NR, NT, 1/)

  START1 = (/iCALL+1/)
  COUNT1 = (/1/)
  time(1) = time_in

  !C PRINT*, 'netcdf_flux_dump: iCALL, time, species, NR, NT, NE, NPA = ', iCALL(i_species), time_in, i_species, NR, NT, NE, NPA

  !       Raw NetCdf (not EzCdf) to allow dimensionality higher than 3
  !       Create NetCdf file and enter define mode

  fileNetcdf = TRIM(PathRamOut)//'/plasma_density_output.nc'
  
  First_time_call : IF (iCALL == 0) THEN ! First time call

     STATUS = NF90_CREATE(fileNetcdf, NF90_CLOBBER, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)	

     !       Define dimensions
     STATUS = NF90_DEF_DIM(NCID, 'rad', NR, RADID)
     STATUS = NF90_DEF_DIM(NCID, 'mlt', NT, MLTID)
     STATUS = NF90_DEF_DIM(NCID, 'time', NF90_UNLIMITED, TIMEID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
    

     ! Define variables
     CALL CHECK ( NF90_DEF_VAR (NCID, 'plasma_density', NF90_FLOAT, (/RADID,MLTID,TIMEID/), DENSITYID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'time', NF90_DOUBLE, TIMEID, TIMEVARID))
     CALL check (nf90_put_att (NCID, TIMEVARID, 'title', 'Time (s)'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'rad', NF90_DOUBLE, RADID, RADVARID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'mlt', NF90_DOUBLE, MLTID, MLTVARID))
    
     !       End define mode
     STATUS = NF90_ENDDEF(NCID)	
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

     ! Only write these at first call
     CALL check ( NF90_PUT_VAR (NCID, RADVARID, L(1:NR)))
     CALL check ( NF90_PUT_VAR (NCID, MLTVARID, mlt))

  ELSE ! Open existing NetCDF
     STATUS = NF90_OPEN (fileNetcdf, NF90_WRITE, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
     CALL check( nf90_inq_dimid(ncid, "rad", RADID))
     CALL check( nf90_inq_dimid(ncid, "mlt", MLTID))
     CALL check( nf90_inq_dimid(ncid, "time", TIMEID))

     ! Only inquire about density and time; other variables (rad,mlt,energy,pitch_angle_cos) do not change
     CALL CHECK ( NF90_INQ_VARID (NCID, 'plasma_density', DENSITYID))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',TIMEVARID))
     !  CALL CHECK ( NF90_ENDDEF(NCID))	

  END IF First_time_call

  ! Convert double to float
  density = density_in

  ! Write NETCDF variables
  CALL check ( NF90_PUT_VAR (NCID, DENSITYID, density, START3, COUNT3))
  CALL check ( NF90_PUT_VAR (NCID, TIMEVARID, time, START1, COUNT1))


  !       Close file
  CALL CHECK ( NF90_CLOSE(NCID) )

  iCALL = iCALL + 1


CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT*, 'STATUS = ', status
       PRINT *, TRIM(nf90_strerror(status))
       STOP 2
    END IF
  END SUBROUTINE check

END SUBROUTINE netcdf_density_dump
