!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModScbCouple

  use ModScbMain,   ONLY: DP
  use ModRamCouple, ONLY: nRadSWMF, nLonSWMF, nRadSwmfVar, nRadSwmfVarInner, &
                          nPoints, iEnd, MhdLines_IIV, xSWMF, ySWMF, zSWMF, &
                          LatSWMF, LonSWMF, xEqSWMF, yEqSWMF, pEqSWMF, nEqSWMF, &
                          nLinesSWMF

  implicit none

  ! Copied from the old Module_points (now removed)
  INTEGER :: last
  INTEGER :: ilat, ilon
  REAL(DP) :: xx01(100000), yy01(100000), zz01(100000)
  REAL(DP), ALLOCATABLE :: xTsyg(:,:,:), yTsyg(:,:,:), zTsyg(:,:,:), &
                           xTsygNew(:,:,:), yTsygNew(:,:,:), zTsygNew(:,:,:)

  INTEGER  :: nthe, ntheOld, iQuiet, iDay, iTime
  REAL(DP) :: rlatval(500), rlonval(500), psival(500), rlon(500), &
              rlat(500), funcModBeta(500), funcModBetaLarge(5000)
  ! Modifies the alpha Euler potential
  REAL(DP) :: radiusStart

  REAL(DP) :: xpsiin ,xpsiout

  contains
!==============================================================================
!  Boundary Conditions for the 3-D Equilibrium Code from SWMF runs 
!  New version - different choice of Euler potentials (so the last Alpha
!  potential surface is a circle in the equatorial plane 

subroutine build_scb_init

  use ModRamFunctions, ONLY: COSD, SIND, ATAN2D, ASIND

  use ModScbGrids, ONLY: npsi, nzeta

  use ModRamGSL, ONLY: GSL_Interpolation_1D, GSL_Interpolation_2D

  use nrtype, ONLY: DP, SP, pi_d

  IMPLICIT NONE

  ! Choices for the "bunching" of points near a certain local 
  ! time (cZ) and equatorial plane (cTheta)
  INTEGER :: GSLerr
  REAL(DP), PARAMETER   :: constZ = 0.0_dp, constTheta = 0.25_dp  ! For domains for inner mag. modeling
  CHARACTER(LEN = 10) :: statusBound ! "open  " or "closed"
  CHARACTER(LEN = 8) :: nthe_char, npsi_char, nzeta_char
  REAL(DP) :: step, r00, tilt, alat, alon, xstop, ystop, zstop, dd, stepMinim, stepMaxim
  REAL(DP), ALLOCATABLE :: alatMod(:,:), aLonVal(:,:)
  INTEGER :: nlat, nlon, halfLon, i, j, k, nMid
  REAL(DP) :: tailMin, tailMax, tailMpause
  REAL(DP) :: Lat, Lon, minLat, maxLat
  REAL(DP) :: distMax(150),  rlonval2Eq(5000)
  real(dp), dimension(1:5000) :: rlonval2Iono, lontwice, betatwice
  REAL(DP) :: magFlux(5000)
  REAL(DP) :: alpha
  REAL(DP), PARAMETER :: xzero = 6.6_dp ! Regular ellipse r=2; 2.5_dp
  INTEGER :: iDebug, iCounter, nlonLarge
  INTEGER, DIMENSION(3) :: dimlens
  INTEGER :: indexOuter, iFieldLine, iPoint, index

  real(DP), dimension(nLinesSWMF,200) :: x, y, z

  LOGICAL, EXTERNAL :: isNaN


  dimlens = (/8, 1, 1/)
  statusBound = 'closed'
  radiusStart = 1._dp  

  iDebug = 1 ! if not zero, extra quantities are written to the CDF file (SWMF, old points, equatorial points etc.)

  nlat = npsi
  nlon = nzeta-1
  halfLon = INT(nzeta/2)

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !cc  Physics Parameters
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  !  PDYN is the solar wind dynamic pressure (0.5*rho*V^2) in unit of nPa
  !  DST, BYIMF, BZIMF are in units of nT

  DO k = 1, nzeta
     ! rlonval(k) = ((k-1._dp)/real(nzeta-1,DP) * 2.*pi_d + constZ * (SIN(2.*pi_d*(k-1.)/real(nzeta-1,DP) + 0.18_dp*pi_d) - &
     !     SIN(0.18_dp*pi_d))) * 180.0_dp /  pi_d ! Centered away from midnight
     rlonval(k) = ((k-1._dp)/REAL(nzeta-1,DP) * 2.*pi_d + constZ * (SIN(2.*pi_d*(k-1.)/REAL(nzeta-1,DP)))) * 180.0_dp /  pi_d 
     ! Centered at midnight, and more clustered there; rlonval is the actual azimuthal angle on the Earth's surface
  END DO

  ! Place variables into easier-to-use containers.
  x = MhdLines_IIV(:,:,3)
  y = MhdLines_IIV(:,:,4)
  z = MhdLines_IIV(:,:,5)

  nlonLarge = nLonSWMF
  DO k = 1, nLonSWMF
     indexOuter = 2*(nRadSWMF*(k-1) + nRadSWMFVar)-1 ! Index for outer perimeter (ellipse) on equatorial plane
     rlonval2Eq(k) = ATAN2(y(indexOuter,1),x(indexOuter,1)) * 180._dp/pi_d 
     IF (rlonval2Eq(k) < 0._dp) rlonval2Eq(k) = rlonval2Eq(k) + 360.
     ! Debug print out for this section:
     !PRINT*, 'k, x, y, rad, rlonval2: ', k, REAL(x(indexOuter,1),sp), REAL(y(indexOuter,1),sp), &
     !     REAL(SQRT(x(indexOuter,1)**2+y(indexOuter,1)**2),sp), REAL(rlonval2Eq(k),sp)
  END DO

  tilt = 0._dp

  ! Tracing step (between Earth surface and 2.5 RE, where SWMF tracing ends)
  stepMinim = 10._dp
  stepMaxim = 100._dp
  step = 2.5_dp 

  !  R00 is the radius of beginning surface for field line tracing; R00 = 1.0 is the Earth surface
  r00 = radiusStart  
  
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !cc  Input Field Line Parameters
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  xstop = 0._dp
  ystop = 0._dp
  zstop = 0._dp
  distMax = 0._dp
  maxLat = 90._dp
  
  ! Values for inner magnetosphere simulations (the same as used in SWMF for tracing B-field lines at the boundary)
  ! Adjust them in SWMF if found to be going outside magnetopause (or if there's an X-line in the tail)
  tailMin = 3.0_dp
  tailMax = 6.75_dp
  tailMpause = 6.75_dp
  
  iCounter = 0
  
  Longitudes_loop_equator: DO ilon = 1, nLonSWMF
     alon = rlonval2Eq(ilon) ! This is angle on equatorial plane 
     
     ! Find the latitude that maps to tailMax 
     indexOuter = 2*(nRadSWMF*(ilon-1) + nRadSWMFVarInner)-1
     Lat = ASIND(z(indexOuter,nPoints)/ &
          SQRT(x(indexOuter,nPoints)**2+y(indexOuter,nPoints)**2+&
          z(indexOuter,nPoints)**2))
     
     ! Find the magnetic flux (old alpha) for this latitude
     magFlux(ilon) = -xzero**3 * (COSD(Lat))**2/radiusStart  ! Lat is in degrees!!!
     
     ! Find the azimuthal angle at the ionosphere (the same as at sphere, as field is dipole from there to Earth)
     Lon = ATAN2D(y(indexOuter,nPoints), x(indexOuter,nPoints))
     IF (Lon < 0. .AND. ilon>nLonSWMF/4) Lon = Lon + 360. 
     IF (ilon > 3*nLonSWMF/4 .AND. Lon<90.) Lon = Lon + 360.
     
     !write(*,'(i4.4, a, 2F10.5)')iLon, ' Lat, Lon = ', Lat, Lon
     
     rlonval2Iono(ilon) = Lon  
     !C        WRITE(*,'(A, I3, 1X, 9F7.2)') 'ilon, lon, xE, yE, zE,' // &
     !C             'Earth lat, Earth lon, magFlux = ', &
     !C             ilon, alon, x(indexOuter,1), &
     !C             y(indexOuter,1), &
     !C             z(indexOuter,1), Lat, Lon, magFlux(ilon)
     maxLat = MIN(maxLat, Lat)
     !PRINT*, 'ilon, Lat on Earth, maxLat = ', ilon, REAL(Lat,sp), REAL(maxLat,sp)
     IF (ABS(Lat-maxLat) < 1.0E-9_dp) nMid = ilon 
     ! Not necessary the "midnight" distance, it can be noon 
     !if the distance on the dayside is small (e.g. eccentric ellipse)
  END DO Longitudes_loop_equator

  !     rlonval2Iono(nlonLarge+1) = rlonval2Iono(1) + 360._dp

  PRINT*, 'nMid = ', nMid
  DO ilon = 1, nlonLarge
     funcModBetaLarge(ilon) = magFlux(nMid) - magFlux(ilon)
  END DO
  
  ! Not needed any more.
  !funcModBetaLarge(nlonLarge+1) = funcModBetaLarge(1)  ! Periodicity
  
  dd = step / 6371.2_dp       
  minLat = 0._dp ! If looking for the highest latitude

  ! Now do the tracing for minLat
  !C     DO ilon = 1, nlonLarge  ! If want to scan all longitudes
  DO ilon = 1, nLonSWMF-1
     ! Find the latitude (highest/lowest) that maps to tailMin     
     index = 2*( (iLon-1)*nRadSWMF + 1 )-1
     Lat = ASIND(z(index,nPoints)/ &
          SQRT(x(index,nPoints)**2+y(index,nPoints)**2+&
          z(index,nPoints)**2))
     !C     PRINT*, 'iLon, index, x, y, z, radius, Lat = ', iLon,  index, &
     !C          x(index,nPoints), y(index,nPoints), z(index,nPoints), &
     !C          SQRT(x(index,nPoints)**2+y(index,nPoints)**2+&
     !C          z(index,nPoints)**2),  Lat
     
     alpha = -xzero**3 * (COSD(Lat))**2 / radiusStart + funcModBetaLarge(ilon) !alpha is the new Euler potential for this latitude
     ! (the magnetic flux for the higher latitude)
     Lat = 180._dp/pi_d * ACOS(SQRT(-radiusStart*alpha/xzero**3)) ! Corrected by funcModBeta
     !C PRINT*, 'ilon, Lat = ', iLon, Lat
     minLat = MAX(minLat, Lat) ! If looking for the highest latitude
     ! minLat = MIN(minLat, Lat) ! If looking for the lowest latitude
  END DO
  
  PRINT*, 'minLat, maxLat = ', minLat, maxLat
  ! Now use interpolation to get (smaller grid) funcModBeta values
  ! For proper interpolation, need to have the whole range [0, 360] covered by rlonval2Iono(k)
  ! If not, add points using periodicity 
  ! CHANGED: For files where first and last Lon are equal (DTW)
  ! For proper interpolation, first and last values in funcModBetaLarge must
  ! be equal and rlonval2Iono must be periodic in degrees.  Additionally,
  ! the full range of lons (0-360) must be covered.  Easy answer: cover
  ! the range twice to ensure periodicity and coverage.
  if (rlonval2IOno(1)<0) then ! We cover 0 degrees but not 360...
     lontwice(1:nLonLarge)  = rlonval2Iono(1:nLonLarge)
     lontwice(nLonLarge+1:2*nLonLarge) = rlonval2Iono(2:nLonLarge+1) + 360.0
     betatwice(1:nLonLarge) = funcModBetaLarge(1:nLonLarge)
     betatwice(nLonLarge+1:2*nLonLarge)= funcModBetaLarge(2:nLonLarge+1)
  else ! We cover 360 but not 0...
     lontwice(1:nLonLarge-1)           = rlonval2Iono(:nLonLarge-1) - 360.0
     lontwice(nLonLarge:2*nLonLarge-1) = rlonval2Iono(1:nLonLarge)
     betatwice(1:nLonLarge-1)          = funcModBetaLarge(:nLonLarge-1)
     betatwice(nLonLarge:2*nLonLarge-1)= funcModBetaLarge(1:nLonLarge)
  end if

  IF (iDebug /= 0) THEN
     DO k = 1, nlonLarge+5
        PRINT*, 'k, lonLarge(k), funcModBetaLarge(k) = ', k, &
             REAL(rlonval2Iono(k),sp), REAL(funcModBetaLarge(k),sp)
     END DO
  END IF
  
  ! These 2 lines below are really problematic if rlonval2Iono and rlonval do not span the same domain !!! 
  ! (which is very possible with azimuthal field curvature present)
  ! funcModBeta(1) = funcModBetaLarge(1)
  ! funcModBeta(nlon+1) = funcModBetaLarge(nlonLarge+1)
  CALL GSL_Interpolation_1D('Cubic',lontwice(1:2*nLonLarge-1), betatwice(1:2*nLonLarge-1), &
                            rlonval(1:nlon), funcModBeta(1:nlon), GSLerr)
  IF (iDebug /= 0) THEN
     PRINT*, ('i, lon(i), funcModBeta(i) = ', i, REAL(rlonval(i),sp), REAL(funcModBeta(i),sp), i = 1, nlon)
  END IF

  CALL latitudes(minLat, maxLat)

  nthe = 2*nPoints - 1 ! Symmetry for the southern hemisphere
  IF (.NOT. ALLOCATED(xTsyg)) ALLOCATE(xTsyg(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(yTsyg)) ALLOCATE(yTsyg(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(zTsyg)) ALLOCATE(zTsyg(nthe,npsi,nzeta+1))
  
  xTsyg = 0._dp
  yTsyg = 0._dp
  zTsyg = 0._dp
  
  IF (.NOT. ALLOCATED(alatMod)) ALLOCATE(alatMod(nlat, nlon))
  IF (.NOT. ALLOCATED(aLonVal)) ALLOCATE(aLonVal(nlat, nlon))
  
  Latitudes_loop: do  ilat = 1, nlat
     alat = rlatVal(ilat) ! Magnetic latitude of the nMid line
     
     !c  longitude grid do loop
     Longitudes_loop: DO  ilon = 1, nlon
        alon = rlonval(ilon)
        !           WRITE(*,*) 'LON =', alon
        
        ! Actual latitude where the field line starts at longitudes other than nMid is different 
        ! (usually higher); to obtain it:
        alpha = -xzero**3 * (cosd(alat))**2 / radiusStart - funcModBeta(ilon) !alpha is the actual magnetic flux for this beta
        alatMod(ilat,ilon) = 180._dp/pi_d * ACOS(SQRT(-radiusStart*alpha/xzero**3)) ! Will correspond to higher latitude than in the F(beta)=0 case
        
     END DO Longitudes_loop
  END DO Latitudes_loop

  ! Rearrange SWMF (+ extra) points so far so that they are equidistant along a field line
  IF (.NOT. ALLOCATED(xSWMF)) ALLOCATE(xSWMF(nthe,nRadSWMF,nlonSWMF-1))
  IF (.NOT. ALLOCATED(ySWMF)) ALLOCATE(ySWMF(nthe,nRadSWMF,nlonSWMF-1))
  IF (.NOT. ALLOCATED(zSWMF)) ALLOCATE(zSWMF(nthe,nRadSWMF,nlonSWMF-1))
  xSWMF = 99.
  ySWMF = 99.
  zSWMF = 99.
  
  DO k = 1, nLonSWMF-1
     DO j = 1, nRadSWMF
        iFieldLine = 2*((k-1)*nRadSWMF + j)-1 ! index of the field line
        xSWMF(1:nPoints, j, k) = x(iFieldLine, 1:nPoints)
        ySWMF(1:nPoints, j, k) = y(iFieldLine, 1:nPoints)
        zSWMF(1:nPoints, j, k) = z(iFieldLine, 1:nPoints)
     END DO
  END DO

  PRINT*, 'minval xEq, maxval xEq = ', MINVAL(xEqSWMF), MAXVAL(xEqSWMF)
  
  ! Southern hemisphere - symmetry
  DO k = 1, nLonSWMF-1
     DO j = 1, nRadSWMF
        DO iPoint = nPoints, 1, -1
           xSWMF(nPoints+iPoint-1,j,k) = xSWMF(iPoint,j,k)
           ySWMF(nPoints+iPoint-1,j,k) = ySWMF(iPoint,j,k)
           zSWMF(nPoints+iPoint-1,j,k) = zSWMF(iPoint,j,k)
        END DO
        
        DO iPoint = 1, nPoints-1 
           xSWMF(iPoint, j, k) = xSWMF(nthe-iPoint+1,j,k)
           ySWMF(iPoint, j, k) = ySWMF(nthe-iPoint+1,j,k)
           zSWMF(iPoint, j, k) = - zSWMF(nthe-iPoint+1,j,k)
        END DO
     END DO
  END DO
  
  CALL mapthe1d(0._dp) ! equidistance for all points along each field line
  
  ! Southern hemisphere - symmetry again
  DO k = 1, nLonSWMF-1
     DO j = 1, nRadSWMF
        DO iPoint = 1, nPoints-1 
           xSWMF(iPoint, j, k) = xSWMF(nthe-iPoint+1,j,k)
           ySWMF(iPoint, j, k) = ySWMF(nthe-iPoint+1,j,k)
           zSWMF(iPoint, j, k) = - zSWMF(nthe-iPoint+1,j,k)
        END DO
     END DO
  END DO

  IF (j == nRadSWMF) PRINT*, 'k, xEq, yEq = ', k, xSWMF(nPoints,j,k), ySWMF(nPoints,j,k)
  
  ! For iPoint = nPoints, we know point coordinates on the Earth's surface (given by the dipole field)
  ! We do this for the southern hemisphere (for northern - symmetry, see below)
  DO k = 2, nlon+1
     DO j = 1, nlat
        xTsyg(1,j,k) = COSD(alatMod(j,k-1)) * COSD(rlonval(k-1))
        yTsyg(1,j,k) = COSD(alatMod(j,k-1)) * SIND(rlonval(k-1))
        zTsyg(1,j,k) = - SIND(alatMod(j,k-1))
        aLonVal(j,k-1) = rLonVal(k-1)
     END DO
  END DO
  ! Do NN interpolation (scattered) to find x, y, z for each such surface
  DO iPoint = 2, nPoints
     !PRINT*, 'cf: iPoint = ', iPoint
     CALL GSL_Interpolation_2D(LatSWMF(:,:),LonSWMF(:,:),REAL(xSWMF(iPoint,:,:),DP), &
                               aLatMod(:,:), aLonVal(:,:), xTsyg(iPoint,:,2:nLon+1), GSLerr)
     CALL GSL_Interpolation_2D(LatSWMF(:,:),LonSWMF(:,:),REAL(ySWMF(iPoint,:,:),DP), &
                               aLatMod(:,:), aLonVal(:,:), yTsyg(iPoint,:,2:nLon+1), GSLerr)
     CALL GSL_Interpolation_2D(LatSWMF(:,:),LonSWMF(:,:),REAL(zSWMF(iPoint,:,:),DP), &
                               aLatMod(:,:), aLonVal(:,:), zTsyg(iPoint,:,2:nLon+1), GSLerr)
  END DO
  !PRINT*, 'cf: Interpolation finished.'
  
  ! Southern hemisphere
  DO k = 2, nlon+1
     !C        PRINT*, 'cf: k, rlonval(k) = ', k, rlonval(k)
     DO j = 1, nlat
        DO iPoint = 1, nPoints-1
           xTsyg(nthe-iPoint+1,j,k) = xTsyg(iPoint,j,k)
           yTsyg(nthe-iPoint+1,j,k) = yTsyg(iPoint,j,k)
           zTsyg(nthe-iPoint+1,j,k) = -zTsyg(iPoint,j,k)
        END DO
     END DO
  END DO
  
  xTsyg(:,:,1) = xTsyg(:,:,nzeta)
  yTsyg(:,:,1) = yTsyg(:,:,nzeta)
  zTsyg(:,:,1) = zTsyg(:,:,nzeta)
  
  ! Downscale number of grid points to smaller nthe for .cdf file
  ntheOld = nthe
  nthe = 51
  
  WRITE(nthe_char, '(I4)') nthe
  WRITE(npsi_char, '(I4)') npsi
  WRITE(nzeta_char, '(I4)') nzeta
  
  IF (.NOT. ALLOCATED(xTsygNew)) ALLOCATE(xTsygNew(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(yTsygNew)) ALLOCATE(yTsygNew(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(zTsygNew)) ALLOCATE(zTsygNew(nthe,npsi,nzeta+1))
  xTsygNew = 99.
  yTsygNew = 99.
  zTsygNew = 99.
  
  CALL mapthe1d_down(constTheta) ! smaller # of points
  
  ! Periodicity
  xTsygNew(:,:,1) = xTsygNew(:,:,nzeta)
  yTsygNew(:,:,1) = yTsygNew(:,:,nzeta)
  zTsygNew(:,:,1) = zTsygNew(:,:,nzeta)
  xTsygNew(:,:,nzeta+1) = xTsygNew(:,:,2)
  yTsygNew(:,:,nzeta+1) = yTsygNew(:,:,2)
  zTsygNew(:,:,nzeta+1) = zTsygNew(:,:,2)
  
  WRITE(*, '(A, I3, 2X, F8.2)') 'hour, max(R) = ', iTime, SQRT(MAXVAL(xTsyg**2 + yTsyg**2));
  
  IF(ALLOCATED(xTsyg) .AND. ALLOCATED(yTsyg) .AND. ALLOCATED(zTsyg)) &
       DEALLOCATE(xTsyg, yTsyg, zTsyg)
  IF(ALLOCATED(xTsygNew) .AND. ALLOCATED(yTsygNew) .AND. ALLOCATED(zTsygNew))&
       DEALLOCATE(xTsygNew, yTsygNew, zTsygNew)
  IF(ALLOCATED(xSWMF) .AND. ALLOCATED(ySWMF) .AND. ALLOCATED(zSWMF))  &
       DEALLOCATE(xSWMF, ySWMF, zSWMF)
  
END subroutine build_scb_init

!==============================================================================
!LOGICAL FUNCTION isNaN(a) 
!  USE nrtype, ONLY : SP, DP
!  REAL(SP) :: a 
!  IF (a.NE.a .OR. a*0._dp.NE.0._dp) THEN 
!     isnan = .TRUE. 
!  ELSE 
!     isnan = .FALSE. 
!  END IF
!  RETURN 
!END FUNCTION isnan

!==============================================================================
SUBROUTINE latitudes(min_lat, max_lat)
  ! const_Z is doing AMR in zeta

  use ModScbGrids, ONLY: npsi

  USE nrtype, ONLY : DP, PI_D

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: min_lat, max_lat
  INTEGER :: nc, j
  REAL(DP), PARAMETER :: pow = 1.0_dp
  REAL(DP), EXTERNAL :: ftor
  REAL(DP) :: xzero3, re, psibc
  REAL(DP) :: xzero, bnormal, psiin, psiout, aa, bb, psitot, xpsiin1, xpsiout1, &
              xpsitot, psis, xpl, cost, cost2, thet

  !NAMELIST/eqdat/nthe,npsi,nzeta,xpsiin,xpsiout

  re = radiusStart ! If the latitudes are NOT from Earth's surface, but from the surface of the sphere of radius radiusStart
  xzero = 6.6

  xzero3 = xzero**3

  xpsiin = re / ((COS(pi_d*min_lat/180._dp))**2)
  xpsiout = re / ((COS(pi_d*max_lat/180._dp))**2)

  !cc.. bnormal is Earth's dipole magnetic field at xzero in equator (in nT)
  !cc.. for xzero = 6.6 R_E, bnormal = 107.83 nT
  bnormal = 0.31_dp / xzero3 * 1.E5_dp

  psiin = -xzero3/xpsiin
  psiout = -xzero3/xpsiout
  psitot = psiout-psiin
  aa = -xzero3
  bb = 0.

  xpsiout1=xpsiout
  xpsiin1=xpsiin

  xpsitot = xpsiout1 - xpsiin1

  !c  Need to make sure that psival is a monotonically increasing function of j
  !c  define psival grids that correspond to dipole psivals for j=1 and j=nps
  !c  define psival grids that correspond to equal equatorial distance grids in the midnight sector

  !c...  the (i,j,k) coordinate system is related to (theta,psis,zeta) flux coordinate
  !c...   by theta = (i-1) * pi / (nthe-1), 0 < theta < pi
  !c...   by zeta = phi + delta(i,j,k), 0 < zeta < 2*pi
  !c...   by psis = (j-1) / (npsi-1), 0 < psis < 1.0
  !c...      psiin < psival < psiout

  DO  j = 1, npsi
     psis = REAL(j-1, DP) / REAL(npsi-1, DP)
     ! xpl = xpsiin1 + xpsitot * 0.1 * (psis + 9.*SIN(0.5*pi_d*psis))
     xpl = xpsiin1 + xpsitot * psis**pow

     !     xpl = xpsiin + 0.5_dp*xpsitot * (psis + sin(0.5_dp*pi_d*psis))
     psival(j) = aa / xpl + bb * (xpl-re)  ! Actual Euler potential value for this surface
  END DO
  !cc corresponding polar angle on ! the tracing start surface; Earth's surface


  re = radiusStart ! added, to compute latitudes for tracing from radiusStart correctly


  nc = 0
  DO j = 1, npsi
     psibc = psival(j)
     cost2 = (re * psibc) / (-xzero3)
     cost = SQRT(cost2)
     thet = ACOS(cost)
     nc = nc + 1
     rlatval(nc) = thet * 90. / (0.5 * pi_d)
  END DO

  ! With the new Euler potential choice, the idea is that we can choose any
  ! latitudes on the midnight 
  ! (or equivalent line), but then for other local times they are found using
  ! the F(beta) function

  RETURN

END SUBROUTINE latitudes

!==============================================================================
SUBROUTINE mapthe1d(const_theta)

  use ModRamFunctions, ONLY: ASIND, ATAN2D, COSD
  USE ModRamGSL,       ONLY: GSL_Interpolation_1D

  use nrtype, ONLY: DP, pi_d

  IMPLICIT NONE

  INTEGER :: GSLerr
  REAL(DP), INTENT(IN)   :: const_theta
  REAL(DP), ALLOCATABLE :: distance(:), xOld(:), yOld(:), zOld(:), distance2derivsX(:), &
                           distance2derivsY(:), distance2derivsZ(:)
  REAL(DP), ALLOCATABLE :: theval(:), thetaVal(:), chiVal(:)
  REAL(DP), ALLOCATABLE :: lonSWMFEarth(:,:)
  INTEGER :: ierr, i, ii, j, k
  LOGICAL, EXTERNAL :: isNaN

  !   move theta coordinates along one line

  ALLOCATE(distance(nthe))
  ALLOCATE(xOld(nthe))
  ALLOCATE(yOld(nthe))
  ALLOCATE(zOld(nthe))
  ALLOCATE(distance2derivsX(nthe))
  ALLOCATE(distance2derivsY(nthe))
  ALLOCATE(distance2derivsZ(nthe))
  ALLOCATE(theval(nthe))
  ALLOCATE(thetaVal(nthe))
  ALLOCATE(chiVal(nthe))
  distance = 0.0; xOld = 0.0; yOld = 0.0; zOld = 0.0
  distance2derivsX = 0.0; distance2derivsY = 0.0; distance2derivsZ = 0.0
  theVal = 0.0; thetaVal = 0.0; chiVal = 0.0

  ierr = 0

  ! print*, 'mapthe1d: nthe = ', nthe

  IF (.NOT. ALLOCATED(LatSWMF)) ALLOCATE(LatSWMF(nRadSWMF,nLonSWMF))
  IF (.NOT. ALLOCATED(LonSWMF)) ALLOCATE(LonSWMF(nRadSWMF,nLonSWMF))
  IF (.NOT. ALLOCATED(LonSWMFEarth)) ALLOCATE(LonSWMFEarth(nRadSWMF,nLonSWMF))
  LatSWMF = 0.0; LonSWMF = 0.0; LonSWMFEarth = 0.0

  DO k = 1, nLonSWMF-1
     DO j = 1, nRadSWMF
        LatSWMF(j,k) = ASIND(zSWMF(nthe,j,k) / &
             SQRT(xSWMF(nthe,j,k)**2+ySWMF(nthe,j,k)**2+&
             zSWMF(nthe,j,k)**2))
        LonSWMF(j,k) = ATAN2D(ySWMF(nthe,j,k), xSWMF(nthe,j,k))
        LonSWMFEarth(j,k) = ATAN2D(ySWMF(nPoints,j,k), xSWMF(nPoints,j,k))
        IF (LonSWMF(j,k) < 0. .AND. k > nLonSWMF/4) LonSWMF(j,k) = LonSWMF(j,k) + 360.
        IF (LonSWMF(j,k) < 90. .AND. k > 3*nLonSWMF/4) LonSWMF(j,k) = LonSWMF(j,k) + 360.
        if (LonSWMFEarth(j,k) < 0.) LonSWMFEarth(j,k) = LonSWMFEarth(j,k) + 360.
        !if (iDebug) WRITE(*,'(A, 2I4, 3F12.5)') 'mapthe1d: k, j, LatI, LonI,
        !LonE = ', &
        !     k, j, LatSWMF(j,k), LonSWMF(j,k), LonSWMFEarth(j,k)

        !       if (j == 14) print*, 'mapthe1d: j, k, x(1,j,k) = ', j, k,
        !       xSWMF(1,j,k)
        distance(1) = 0._dp
        !        PRINT*, 'mapthe1d: j, k, 2, x, y, z: ', j, k,
        !        REAL(xSWMF(2,j,k),sp), &
        !             REAL(ySWMF(2,j,k),sp), REAL(zSWMF(2,j,k),sp)
        xOld(:) = xSWMF(:,j,k)
        yOld(:) = ySWMF(:,j,k)
        zOld(:) = zSWMF(:,j,k)

        DO  i = 2, nthe
           distance(i) = distance(i-1) + SQRT((xSWMF(i,j,k)-xSWMF(i-1,j,k))**2 + &
                (ySWMF(i,j,k)-ySWMF(i-1,j,k))**2+(zSWMF(i,j,k)-zSWMF(i-1,j,k))**2)
           !C         if (j == 14 .and. k == 4) PRINT*, 'mapthe1d: j, k, i, x,
           !y, z, distance, diff: ', j, k, i, REAL(xOld(i),sp), &
           !C              REAL(yOld(i),sp), REAL(zOld(i),sp),
           !REAL(distance(i),sp), real(distance(i)-distance(i-1),sp)
        END DO

        DO  ii = 1, nthe
           theval(ii) = distance(nthe) * REAL(ii - 1,dp) / REAL(nthe - 1,dp)
        END DO

        thetaVal = theval * pi_d / distance(nthe)
        chiVal = (thetaVal + const_theta * SIN(2.*thetaVal)) * distance(nthe)/pi_d

        CALL GSL_Interpolation_1D('Cubic',distance(:),xOld(:),chiVal(2:nthe-1),xSWMF(2:nthe-1,j,k),GSLerr)
        CALL GSL_Interpolation_1D('Cubic',distance(:),yOld(:),chiVal(2:nthe-1),ySWMF(2:nthe-1,j,k),GSLerr)
        CALL GSL_Interpolation_1D('Cubic',distance(:),zOld(:),chiVal(2:nthe-1),zSWMF(2:nthe-1,j,k),GSLerr)

     END DO
     !C if (k == 4) print*, 'mapthe1d: xEq(:,k) = ', xSWMF(nPoints,:,k)
  END DO

  DEALLOCATE(distance, xOld, yOld, zOld, distance2derivsX, distance2derivsY, distance2derivsZ, &
       theval, thetaVal, chiVal)

  RETURN

END SUBROUTINE mapthe1d

!==============================================================================
SUBROUTINE mapthe1d_down(const_theta)

  use ModRamGSL,   ONLY: GSL_Interpolation_1D
  use ModScbGrids, ONLY: nzeta, npsi, nthe

  use nrtype, ONLY: DP, pi_d

  IMPLICIT NONE

  INTEGER :: GSLerr
  REAL(DP), INTENT(IN)   :: const_theta
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: distance, xOld, yOld, zOld, distance2derivsX, &
                                         distance2derivsY, distance2derivsZ
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: theval, thetaVal, chiVal
  INTEGER :: ierr, i, ii, j, k
  LOGICAL, EXTERNAL :: isNaN

  ALLOCATE(distance(ntheOld))
  ALLOCATE(xOld(ntheOld))
  ALLOCATE(yOld(ntheOld))
  ALLOCATE(zOld(ntheOld))
  ALLOCATE(distance2derivsX(ntheOld))
  ALLOCATE(distance2derivsY(ntheOld))
  ALLOCATE(distance2derivsZ(ntheOld))

  ALLOCATE(theval(nthe))
  ALLOCATE(thetaVal(nthe))
  ALLOCATE(chiVal(nthe))

  !   move theta coordinates along one line

  ierr = 0

  DO k = 2, nzeta
     DO j = 1, npsi
        distance(1) = 0._dp
        xOld(:) = xTsyg(:,j,k)
        yOld(:) = yTsyg(:,j,k)
        zOld(:) = zTsyg(:,j,k)
        DO  i = 2, ntheOld
           distance(i) = distance(i-1) + SQRT((xTsyg(i,j,k)-xTsyg(i-1,j,k))**2 + &
                (yTsyg(i,j,k)-yTsyg(i-1,j,k))**2+(zTsyg(i,j,k)-zTsyg(i-1,j,k))**2)
        END DO

        ! New coordinates
        DO  ii = 1, nthe
           theval(ii) = distance(ntheOld) * REAL(ii - 1,dp) / REAL(nthe - 1,dp)
        END DO
        thetaVal = theval * pi_d / distance(ntheOld)
        chiVal = (thetaVal + const_theta * SIN(2.*thetaVal)) * distance(ntheOld)/pi_d ! dimensions of distance

        CALL GSL_Interpolation_1D('Cubic',distance(1:ntheOld),xOld(1:ntheOld), &
                                  chiVal(2:nthe-1),xTsygNew(2:nthe-1,j,k),GSLerr)
        CALL GSL_Interpolation_1D('Cubic',distance(1:ntheOld),yOld(1:ntheOld), &
                                  chiVal(2:nthe-1),yTsygNew(2:nthe-1,j,k),GSLerr)
        CALL GSL_Interpolation_1D('Cubic',distance(1:ntheOld),zOld(1:ntheOld), &
                                  chiVal(2:nthe-1),zTsygNew(2:nthe-1,j,k),GSLerr)

        xTsygNew(1,j,k) = xTsyg(1,j,k)
        xTsygNew(nthe,j,k) = xTsyg(ntheOld,j,k)
        yTsygNew(1,j,k) = yTsyg(1,j,k)
        yTsygNew(nthe,j,k) = yTsyg(ntheOld,j,k)
        zTsygNew(1,j,k) = zTsyg(1,j,k)
        zTsygNew(nthe,j,k) = zTsyg(ntheOld,j,k)
     END DO
  END DO

  DEALLOCATE(distance)
  DEALLOCATE(xOld)
  DEALLOCATE(yOld)
  DEALLOCATE(zOld)
  DEALLOCATE(distance2derivsX)
  DEALLOCATE(distance2derivsY)
  DEALLOCATE(distance2derivsZ)
  DEALLOCATE(theval)
  DEALLOCATE(thetaVal)
  DEALLOCATE(chiVal)

  RETURN
END SUBROUTINE mapthe1d_down

End Module ModScbCouple
