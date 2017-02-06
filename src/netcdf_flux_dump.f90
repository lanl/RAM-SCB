!****************************************************************************
! Dumps particle flux in NetCdf file (time included as unlimited dimension)
! Author: Sorin Zaharia, Oct. 2009
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!****************************************************************************

SUBROUTINE netcdf_flux_dump(flux_in, NR, rad_in, NT, phi_in, NE, ekev_in, &
     NPA, mu_in, i_species, time_in)

  USE NETCDF
  USE ModRamMain, ONLY : PathRamOut, Real8_, Real4_

  IMPLICIT NONE

  INTEGER :: NCID, STATUS, RADID, PHIID, ENERGID, PAID, FLUXID, DUMMYID, TIMEID, &
       RADVARID, PHIVARID, ENERGVARID, PAVARID, TIMEVARID
  REAL(kind=Real8_), INTENT(IN) :: rad_in(NR+1), phi_in(NT), ekev_in(NE), mu_in(NPA)
  REAL(kind=Real4_), INTENT(IN) :: FLUX_IN(NR, NT, NE, NPA)
  REAL(kind=Real8_), INTENT(IN) :: time_in
  REAL(kind=Real8_) :: time(1) 
  REAL(kind=Real4_) :: flux(NR,NT,NE,NPA)

  INTEGER :: START(5), COUNT(5) ! netCDF variable start point
  INTEGER :: START2(1), COUNT2(1)
  INTEGER, INTENT(IN) :: NR, NT, NE, NPA
  INTEGER, INTENT(IN) :: i_species
  INTEGER :: ier
  CHARACTER(LEN=3) :: sp

  CHARACTER(LEN=200) :: fileNetcdf

  INTEGER, SAVE :: iCALL(4) = 0

  START = (/1, 1, 1, 1, iCALL(i_species)+1/)
  COUNT = (/NR,NT,NE,NPA,1/)

  START2 = (/iCALL(i_species)+1/)
  COUNT2 = (/1/)
  time(1) = time_in

  !C PRINT*, 'netcdf_flux_dump: iCALL, time, species, NR, NT, NE, NPA = ', iCALL(i_species), time_in, i_species, NR, NT, NE, NPA

  !       Raw NetCdf (not EzCdf) to allow dimensionality higher than 3
  !       Create NetCdf file and enter define mode

  SELECT CASE (i_species)
  CASE(1)
     sp = '_e-'
  CASE (2)
     sp = '_H+'
  CASE(3)
     sp = 'He+'
  CASE(4)
     sp = '_O+'
  CASE default
     STOP 'netcdf_flux_dump: select_case problem.'
  END SELECT

  fileNetcdf = TRIM(PathRamOut)//'/flux_output'//sp//'.nc'
  
  First_time_call : IF (iCALL(i_species) == 0) THEN ! First time call

     STATUS = NF90_CREATE(fileNetcdf, NF90_CLOBBER, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)	

     !       Define dimensions
     STATUS = NF90_DEF_DIM(NCID, 'rad', NR, RADID)
     STATUS = NF90_DEF_DIM(NCID, 'phi', NT, PHIID)
     STATUS = NF90_DEF_DIM(NCID, 'energy', NE, ENERGID)
     STATUS = NF90_DEF_DIM(NCID, 'pitch_angle_cos', NPA, PAID)
     STATUS = NF90_DEF_DIM(NCID, 'time', NF90_UNLIMITED, TIMEID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
    

     ! Define variables
     CALL CHECK ( NF90_DEF_VAR (NCID, 'RAM_flux', NF90_FLOAT, (/RADID,PHIID,ENERGID,PAID,TIMEID/), FLUXID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'time', NF90_DOUBLE, TIMEID, TIMEVARID))
     CALL check (nf90_put_att (NCID, TIMEVARID, 'title', 'Time (s)'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'rad', NF90_DOUBLE, RADID, RADVARID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'phi', NF90_DOUBLE, PHIID, PHIVARID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'energy', NF90_DOUBLE, ENERGID, ENERGVARID))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'pitch_angle_cos', NF90_DOUBLE, PAID, PAVARID))

    
     !       End define mode
     STATUS = NF90_ENDDEF(NCID)	
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

     ! Only write these at first call
     CALL check ( NF90_PUT_VAR (NCID, RADVARID, rad_in(1:NR)))
     CALL check ( NF90_PUT_VAR (NCID, PHIVARID, phi_in))
     CALL check ( NF90_PUT_VAR (NCID, ENERGVARID, ekev_in))
     CALL check ( NF90_PUT_VAR (NCID, PAVARID, mu_in))

  ELSE ! Open existing NetCDF
     STATUS = NF90_OPEN (fileNetcdf, NF90_WRITE, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
     CALL check( nf90_inq_dimid(ncid, "rad", RADID))
     CALL check( nf90_inq_dimid(ncid, "phi", PHIID))
     CALL check( nf90_inq_dimid(ncid, "energy", ENERGID))
     CALL check( nf90_inq_dimid(ncid, "pitch_angle_cos", PAID))
     CALL check( nf90_inq_dimid(ncid, "time", TIMEID))

     ! Only inquire about flux and time; other variables (rad,phi,energy,pitch_angle_cos) do not change
     CALL CHECK ( NF90_INQ_VARID (NCID, 'RAM_flux', FLUXID))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',TIMEVARID))
     !  CALL CHECK ( NF90_ENDDEF(NCID))	

  END IF First_time_call

  ! Write NETCDF variables
  CALL check ( NF90_PUT_VAR (NCID, FLUXID, flux_in, START, COUNT))
  CALL check ( NF90_PUT_VAR (NCID, TIMEVARID, time, START2, COUNT2))


  !       Close file
  CALL CHECK ( NF90_CLOSE(NCID) )

  iCALL(i_species) = iCALL(i_species) + 1


CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT*, 'STATUS = ', status
       PRINT *, TRIM(nf90_strerror(status))
       STOP 2
    END IF
  END SUBROUTINE check

END SUBROUTINE netcdf_flux_dump
