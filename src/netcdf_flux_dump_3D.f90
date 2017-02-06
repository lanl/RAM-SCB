!******************************************************************************
! Dumps RAM particle flux in NetCdf file (time included as unlimited dimension)
! Author: Sorin Zaharia, Oct. 2009
! Only at half the pitch angles (size limitations)
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!******************************************************************************

SUBROUTINE netcdf_flux_dump_3D 
  USE NETCDF
  USE MPI
  USE NRTYPE, ONLY : SP, DP
  USE Module1, ONLY : rHour, nthe, npsi, nzeta, chiVal, alphaVal, psiVal, bnormal
  USE Module_RAM, ONLY : flux3DEQ, indexPA, NRAD, NT, NE, NPA
  USE ModRamMain, ONLY : PathScbOut, FLUX, MU, EKEV

  IMPLICIT NONE

  INTEGER :: NCID, STATUS, ENERGID, PAID, FLUXHID, FLUXOID, DUMMYID, TIMEID, &
       RADVARID, PHIVARID, ENERGVARID, PAVARID, TIMEVARID, &
       chiid, alphaid, betaid, chivarid, alphavarid, betavarid
  REAL(SP) :: time(1) 

  REAL(SP) :: start_time, end_time
  INTEGER :: START(6), COUNT(6) ! netCDF variable start point
  INTEGER :: START2(1), COUNT2(1)
  INTEGER :: START2D(2)
  INTEGER :: ier, i, j, k, L, iSpecies

  CHARACTER(LEN=200) :: fileNetcdf

  REAL(SP), ALLOCATABLE :: fluxFullH(:,:,:,:,:) 
  REAL(SP), ALLOCATABLE :: fluxFullO(:,:,:,:,:) 

  INTEGER, SAVE :: iCALL = 0
  
!C  PRINT*, 'Trying to allocate arrays of size ', 8.*REAL(nthe,sp)*npsi*nzeta*NE*NPA/2/1024.**2, ' Mbytes.' ! (has to be less than 2Gb for 32-bit arch).
  IF (.NOT. ALLOCATED(fluxFullH)) ALLOCATE(fluxFullH(nthe,npsi,nzeta,NE,NPA/2), STAT=ier)
  IF (.NOT. ALLOCATED(fluxFullO)) ALLOCATE(fluxFullO(nthe,npsi,nzeta,NE,NPA/2), STAT=ier)

  start_time = MPI_WTIME()
  ! do iSpecies = 2, 4
  DO L = 1, NPA/2
     DO k = 1, nzeta
        DO j = 1, npsi
           DO i = 1, nthe
              IF (indexPA(i,j,k,2*L-1) /= -1) THEN
                 fluxFullH(i,j,k,:,L) = flux3DEQ(2,i,j,k,indexPA(i,j,k,2*L-1))
                 fluxFullO(i,j,k,:,L) = flux3DEQ(4,i,j,k,indexPA(i,j,k,2*L-1))
              ELSE
                 fluxFullH(i,j,k,:,L) = -1.
                 fluxFullO(i,j,k,:,L) = -1.
              END IF
           END DO
        END DO
     END DO
  END DO
  ! end do
  end_time = MPI_WTIME()
 ! PRINT*, 'nfd3D: time  = ', end_time - start_time

  START = (/1, 1, 1, 1, 1, iCALL+1/)
  COUNT = (/nthe,npsi,nzeta,NE,NPA/2,1/)

  START2D = (/1,iCALL+1/)

  START2 = (/iCALL+1/)
  COUNT2 = (/1/)
  time(1) = rHour

  !       Raw NetCdf (not EzCdf) to allow dimensionality higher than 3
  !       Create NetCdf file and enter define mode

  fileNetcdf = TRIM(PathScbOut)//'/flux_full_output.nc'

  First_time_call : IF (iCALL == 0) THEN ! First time call

     STATUS = NF90_CREATE(fileNetcdf, NF90_64BIT_OFFSET, NCID)
     ! STATUS = NF90_CREATE(fileNetcdf, NF90_CLOBBER, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)	

     !       Define dimensions
     STATUS = NF90_DEF_DIM(NCID, 'chi', nthe, CHIID)
     STATUS = NF90_DEF_DIM(NCID, 'alpha', npsi, ALPHAID)
     STATUS = NF90_DEF_DIM(NCID, 'beta', nzeta, BETAID)
     STATUS = NF90_DEF_DIM(NCID, 'energy', NE, ENERGID)
     STATUS = NF90_DEF_DIM(NCID, 'pitch_angle_cos', NPA/2, PAID)
     STATUS = NF90_DEF_DIM(NCID, 'time', NF90_UNLIMITED, TIMEID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)


     ! Define variables
     CALL CHECK ( NF90_DEF_VAR (NCID, 'RAM_flux_H+', NF90_FLOAT, &
          (/CHIID,ALPHAID,BETAID,ENERGID,PAID,TIMEID/), FLUXHID))
   !  CALL check(nf90_put_att(ncid,chivarid,'missing_value', -1))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'RAM_flux_O+', NF90_FLOAT, &
          (/CHIID,ALPHAID,BETAID,ENERGID,PAID,TIMEID/), FLUXOID))
   !  CALL check(nf90_put_att(ncid,chivarid,'missing_value', -1))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'time', NF90_FLOAT, TIMEID, TIMEVARID))
     CALL check (nf90_put_att (NCID, TIMEVARID, 'title', 'Time (s)'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'chi', NF90_FLOAT, CHIID, CHIVARID))
     CALL check(nf90_put_att(ncid,chivarid,'title','Coordinate along field line'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'alpha', NF90_FLOAT, ALPHAID, ALPHAVARID))
     CALL check(nf90_put_att(ncid,alphavarid,'title','Magnetic flux-like Euler potential'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'beta', NF90_FLOAT, BETAID, BETAVARID))
     CALL check(nf90_put_att(ncid,betavarid,'title','Azimuthal angle-like Euler potential'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'energy', NF90_FLOAT, ENERGID, ENERGVARID))
     CALL check(nf90_put_att(ncid,energvarid,'title','RAM energy space (keV)'))
     CALL CHECK ( NF90_DEF_VAR (NCID, 'pitch_angle_cos', NF90_FLOAT, PAID, PAVARID))
     CALL check(nf90_put_att(ncid,pavarid,'title','RAM cosine of pitch angle'))

     !       End define mode
     STATUS = NF90_ENDDEF(NCID)	
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)

     ! Only write these at first call
     CALL check ( NF90_PUT_VAR (NCID, ENERGVARID, REAL(EKEV,SP)))
     CALL check ( NF90_PUT_VAR (NCID, PAVARID, REAL(MU(1:NPA:2),SP)))

  ELSE ! Open existing NetCDF
     STATUS = NF90_OPEN (fileNetcdf, NF90_WRITE, NCID)
     IF (STATUS .NE. NF90_NOERR) CALL HANDLE_ERR(STATUS)
     
   !  CALL check( nf90_inq_dimid(ncid, 'energy', energID))
   !  CALL check( nf90_inq_dimid(ncid, 'pitch_angle_cos', PAID))
   !  CALL check( nf90_inq_dimid(ncid, 'time', TIMEID))

     ! Only inquire about flux, time, alpha, beta, chi; other variables (rad,phi,energy,pitch_angle_cos) do not change
     CALL check( nf90_inq_varid(ncid, 'alpha', alphaVarID))
     CALL check( nf90_inq_varid(ncid, 'beta', betaVarID))
     CALL check( nf90_inq_varid(ncid, 'chi', chiVarID))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'RAM_flux_H+', FLUXHID))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'RAM_flux_O+', FLUXOID))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',TIMEVARID))
     !  CALL CHECK ( NF90_ENDDEF(NCID))	

  END IF First_time_call

  ! Write NETCDF variables
  CALL check(nf90_put_var(ncid, chivarid, REAL(chiVal(1:nthe)), START2D))
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START2D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START2D))

  CALL check ( NF90_PUT_VAR (NCID, FLUXHID, fluxFullH, START, COUNT))
  CALL check ( NF90_PUT_VAR (NCID, FLUXHID, fluxFullO, START, COUNT))
  CALL check ( NF90_PUT_VAR (NCID, TIMEVARID, time, START2, COUNT2))


  !       Close file
  CALL CHECK ( NF90_CLOSE(NCID) )

  iCALL = iCALL + 1

  IF (ALLOCATED(fluxFullH)) DEALLOCATE(fluxFullH, stat=ier)
  IF (ALLOCATED(fluxFullO)) DEALLOCATE(fluxFullO, stat=ier)


CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT*, 'STATUS = ', status
       PRINT *, TRIM(nf90_strerror(status))
       STOP 2
    END IF
  END SUBROUTINE check

END SUBROUTINE netcdf_flux_dump_3D
