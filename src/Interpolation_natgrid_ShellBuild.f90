SUBROUTINE Interpolation_natgrid_Shellbuild(quantitySWMF, lat_local, lon_local, iPoint, nlat, nlon, quantityInterp)

  ! gives values on polar grid, using interpolation

  ! .. Use Statements ..
  USE nrtype, ONLY : SP, DP, pi_d
  USE ModRamCouple, ONLY : nRadSWMF, nLonSWMF, xSWMF, ySWMF, &
       zSWMF, LatSWMF, LonSWMF

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: lat_local(:,:), lon_local(:)
  REAL(DP), INTENT(IN) :: quantitySWMF(:,:)
  REAL(DP), INTENT(OUT) :: quantityInterp(nlat, nlon)
  INTEGER, INTENT(IN) :: iPoint, nlat, nlon

  INTEGER :: ix, j, jy, k, iTotal, ierr, nTotal

  real(dp) :: LatInner(nLonSWMF)
  REAL(DP), ALLOCATABLE :: latScatter(:), lonScatter(:), quantityScatter(:)
  !-----------------------------------------------------------------------

  ! Set total number of points, including extra to provide for continuity 
  ! across 360 to 0 degree boundary and inner boundary ghost cells.
  nTotal = (nRadSWMF+1)*(nLonSWMF+20)

  ! Create inner latitude values that map inside of inner boundary:
  LatInner = 2.*LatSWMF(1,:) - LatSWMF(2,:)

  ! Allocate variables to hold unordered points:
  IF (.NOT. ALLOCATED(latScatter)) ALLOCATE(latScatter(nTotal), stat = ierr)
  IF (.NOT. ALLOCATED(lonScatter)) ALLOCATE(lonScatter(nTotal), stat = ierr)
  IF (.NOT. ALLOCATED(quantityScatter)) ALLOCATE(quantityScatter(nTotal), &
       stat = ierr)  

  ! Begin tracking the current point:
  iTotal = 0

  ! Unpack grid points into 1D vector:
  k_Loop:  DO k = 1, nLonSWMF-1 
     ! Add "safety layer" at inner radius:
     iTotal = iTotal + 1
     latScatter(iTotal) = LatInner(k)
     lonScatter(iTotal) = LonSWMF(1,k)
     quantityScatter(iTotal) = quantitySWMF(1,k)
     ! Add real points outside of inner ring:
     j_loop: DO j = 1, nRadSWMF
        iTotal = iTotal+1
        latScatter(iTotal) = LatSWMF(j,k)
        lonScatter(iTotal) = LonSWMF(j,k)
        quantityScatter(iTotal) = quantitySWMF(j,k)
     END DO j_loop
  END DO k_Loop

  ! Add another set of points for 1 and nLonSWMF.  Small negative longitude 
  ! will ensure that neighbors are found for 0 degrees longitude if field 
  ! lines are curved, and no extrapolation is needed)
  DO k = nLonSWMF-10, nLonSWMF-1 ! 10 longitudes
      ! Add "safety layer" at inner radius:
     iTotal = iTotal + 1
     latScatter(iTotal) = LatInner(k)
     lonScatter(iTotal) = LonSWMF(1,k)
     quantityScatter(iTotal) = quantitySWMF(1,k)
     DO j = 1, nRadSWMF     
        iTotal = iTotal+1
        latScatter(iTotal) = LatSWMF(j,k)
        lonScatter(iTotal) = LonSWMF(j,k) - 360.
        quantityScatter(iTotal) = quantitySWMF(j,k)
     END DO
  END DO
  DO k = 1, 10 ! 10 longitudes
     ! Add "safety layer" at inner radius:
     iTotal = iTotal + 1
     latScatter(iTotal) = LatInner(k)
     lonScatter(iTotal) = LonSWMF(1,k)
     quantityScatter(iTotal) = quantitySWMF(1,k)
     DO j = 1, nRadSWMF     
        iTotal = iTotal+1
        latScatter(iTotal) = LatSWMF(j,k)
        lonScatter(iTotal) = LonSWMF(j,k) + 360.
        quantityScatter(iTotal) = quantitySWMF(j,k)
     END DO
  END DO

  ! Initialize delaunay triangulation:
  CALL NNPNTINITD(iTotal, latScatter, lonScatter, quantityScatter)

  ! Interpolate:
  DO ix = 1, nlat
     DO jy = 1, nlon
        CALL NNPNTD(lat_local(ix,jy),lon_local(jy),quantityInterp(ix,jy))
     END DO
  END DO

  ! Clean up interpolator:
  CALL NNPNTENDD()

  ! Deallocate arrays.
  IF (ALLOCATED(latScatter)) DEALLOCATE(latScatter, stat = ierr)
  IF (ALLOCATED(lonScatter)) DEALLOCATE(lonScatter, stat = ierr)
  IF (ALLOCATED(quantityScatter)) DEALLOCATE(quantityScatter, stat = ierr)

  RETURN


END SUBROUTINE Interpolation_natgrid_Shellbuild


