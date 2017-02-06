!  Boundary Conditions for the 3-D Equilibrium Code from SWMF runs 
!  New version - different choice of Euler potentials (so the last Alpha
!  potential surface is a circle in the equatorial plane 

subroutine build_scb_init

  USE nrtype, ONLY : dp, sp, pi_d
  USE Module_points
  USE ezcdf
  use ModRamFunctions
  use ModRamCouple, ONLY: nRadSWMF, nLonSWMF, nRadSwmfVar, nRadSwmfVarInner, &
       nPoints, iEnd, MhdLines_IIV, &
       xSWMF, ySWMF, zSWMF, LatSWMF, LonSWMF, &
       xEqSWMF, yEqSWMF, pEqSWMF, nEqSWMF, nLinesSWMF

  IMPLICIT NONE

  INTERFACE spline
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
  INTERFACE splint
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

  INTERFACE Interpolation_natgrid_ShellBuild
     SUBROUTINE Interpolation_natgrid_ShellBuild(quantitySWMF, lat_local, lon_local, iPoint, nlat, nlon, quantityInterp)
       USE nrtype, ONLY : DP, pi_d
       USE ModRamCouple, ONLY : nRadSWMF, nLonSWMF, xSWMF, ySWMF, zSWMF
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: lat_local(:,:), lon_local(:)
       REAL(DP), INTENT(IN) :: quantitySWMF(:,:)
       REAL(DP), INTENT(OUT) :: quantityInterp(nlat, nlon)
       INTEGER, INTENT(IN) :: iPoint, nlat, nlon
     END SUBROUTINE Interpolation_natgrid_ShellBuild
  END INTERFACE

  INTERFACE Spline_1D
     SUBROUTINE Spline_1D_periodic(x1, f, z1, func, iErrDomain)
       USE EZspline_obj ! import the modules
       USE EZspline  
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind; same as DP in nrtype
       REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8
       REAL(r8), DIMENSION(:), INTENT(IN) :: x1 ! independent variable
       REAL(r8), DIMENSION(:), INTENT(IN) :: f
       REAL(r8), DIMENSION(:), INTENT(IN) :: z1
       REAL(r8), DIMENSION(:), INTENT(OUT) :: func ! interpolated values
       integer, intent(OUT) :: iErrDomain
       INTEGER n1, ier, BCS1(2)
       TYPE(EZspline1_r8) :: spline_o ! 1-D EZspline object
       INTEGER k1
     END SUBROUTINE Spline_1D_periodic
  END INTERFACE


  ! Choices for the "bunching" of points near a certain local 
  ! time (cZ) and equatorial plane (cTheta)
  REAL(DP), PARAMETER   :: constZ = 0.0_dp, constTheta = 0.25_dp  ! For domains for inner mag. modeling
  INTEGER, PARAMETER :: noHrs = 24*30-1 ! Number of hours (snapshots) ! For Nov. 2002 run
  CHARACTER(LEN = 10) :: statusBound ! "open  " or "closed"
  CHARACTER(LEN = 8) :: pdyn_char, byimf_char, bzimf_char, dst_char, r00_char, xpsiin_char, &
       xpsiout_char, constz_char, consttheta_char, nthe_char, npsi_char, nzeta_char
  CHARACTER(LEN = 4) :: time_char
  REAL(DP)           :: step, r00, tilt, alat, alon, rtilt, x0, y0, z0, &
       xstop, ystop, zstop, dd, stepMinim, stepMaxim
  REAL(DP), ALLOCATABLE :: alatMod(:,:)
  INTEGER :: nlat, nlon, halfLon, i, j, k, netcdfId, nMid
  REAL(DP) :: start_time, end_time
  REAL(DP) :: tailMin, tailMax, tailMpause, aMajor, eccent, radius
  REAL(DP) :: Lat, LatSphere, Lon, minLat, maxLat, rSphere
  REAL(DP) :: distMax(150),  rlonval2Eq(5000), distance2derivs(5000)
  real(dp), dimension(1:5000) :: rlonval2Iono, lontwice, betatwice
  REAL(DP) :: magFlux(5000)
  REAL(DP) :: alpha
  REAL(DP), PARAMETER :: xzero = 6.6_dp, rSuper = 2.0_dp ! Regular ellipse r=2; 2.5_dp
  REAL(DP) :: rXY, rXYSWMF, rad
  INTEGER :: ierr, iDebug, iCounter, nlonLarge
  INTEGER, DIMENSION(3) :: dimlens
  INTEGER, DIMENSION(3) :: dimlens2
  INTEGER :: indexOuter, iFieldLine, iPoint, iRad, index

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
     rlonval2Eq(k) = DATAN2(y(indexOuter,1),x(indexOuter,1)) * 180._dp/pi_d 
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
25 CONTINUE      

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
  ! Now use spline interpolation to get (smaller grid) funcModBeta values
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
  
  !C CALL spline(rlonval2Iono(1:nlonLarge+5),funcModBetaLarge(1:nlonLarge+5), &
  !C      1.E31_DP, 1.E31_DP, distance2derivs(1:nlonLarge+5))
  !C DO i = 1, nlon
  !C    funcModBeta(i) = splint(rlonval2Iono(1:nlonLarge+5),funcModBetaLarge(1:nlonLarge+5), &
  !C         distance2derivs(1:nlonLarge+5), rlonval(i))
  !C IF (iDebug /= 0) PRINT*, 'i, lon(i), funcModBeta(i) = ', i, REAL(rlonval(i),sp), REAL(funcModBeta(i),sp)
  !C END DO 
  
  ! These 2 lines below are really problematic if rlonval2Iono and rlonval do not span the same domain !!! 
  ! (which is very possible with azimuthal field curvature present)
  ! funcModBeta(1) = funcModBetaLarge(1)
  ! funcModBeta(nlon+1) = funcModBetaLarge(nlonLarge+1)
  CALL Spline_1D(lontwice(1:2*nLonLarge-1), betatwice(1:2*nLonLarge-1), &
       rlonval(1:nlon), funcModBeta(1:nlon), ierr)
  IF (iDebug /= 0) THEN
     PRINT*, ('i, lon(i), funcModBeta(i) = ', i, REAL(rlonval(i),sp), REAL(funcModBeta(i),sp), i = 1, nlon)
  END IF

  CALL latitudes(minLat, maxLat)

  nthe = 2*nPoints - 1 ! Symmetry for the southern hemisphere
  IF (.NOT. ALLOCATED(xTsyg)) ALLOCATE(xTsyg(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(yTsyg)) ALLOCATE(yTsyg(nthe,npsi,nzeta+1))
  IF (.NOT. ALLOCATED(zTsyg)) ALLOCATE(zTsyg(nthe,npsi,nzeta+1))
  xTsyg = 99.
  yTsyg = 99.
  zTsyg = 99.
  
  WRITE(xpsiin_char, '(SP, F6.2)') xpsiin
  WRITE(xpsiout_char, '(SP, F6.2)') xpsiout
  WRITE(r00_char, '(SP, F8.4)') r00
  WRITE(constz_char, '(SP, F6.2)') constZ
  WRITE(consttheta_char, '(SP, F6.2)') constTheta
  
  CALL cdf_open(netcdfId, 'SWMF_config_test.cdf', 'w')
  
  dimlens(1) = LEN(statusBound)
  CALL  cdf_define (netcdfId, 'statusBC', dimlens, 'char')
  dimlens(1) = LEN(pdyn_char)
  CALL  cdf_define (netcdfId, 'pdyn', dimlens, 'char')
  dimlens(1) = LEN(byimf_char)
  CALL  cdf_define (netcdfId, 'byimf', dimlens, 'char')
  dimlens(1) = LEN(bzimf_char)
  CALL  cdf_define (netcdfId, 'bzimf', dimlens, 'char')
  dimlens(1) = LEN(dst_char)
  CALL  cdf_define (netcdfId, 'DST', dimlens, 'char')
  dimlens(1) = LEN(xpsiin_char)
  CALL  cdf_define (netcdfId, 'xpsi_in', dimlens, 'char')
  dimlens(1) = LEN(xpsiout_char)
  CALL  cdf_define (netcdfId, 'xpsi_out', dimlens, 'char')
  dimlens(1) = LEN(r00_char)
  CALL  cdf_define (netcdfId, 'rStart', dimlens, 'char')
  dimlens(1) = LEN(constz_char)
  CALL  cdf_define (netcdfId, 'Constz', dimlens, 'char')
  dimlens(1) = LEN(consttheta_char)
  CALL  cdf_define (netcdfId, 'Consttheta', dimlens, 'char')
  dimlens(1) = LEN(nthe_char)
  CALL  cdf_define (netcdfId, 'ntheta', dimlens, 'char')
  dimlens(1) = LEN(npsi_char)
  CALL  cdf_define (netcdfId, 'npsi', dimlens, 'char')
  dimlens(1) = LEN(nzeta_char)
  CALL  cdf_define (netcdfId, 'nzeta', dimlens, 'char')
  CALL cdfSetatt(netcdfId, 'n_zeta', 'Number of zeta points')
  
  xTsyg = 0._dp
  yTsyg = 0._dp
  zTsyg = 0._dp
  
  dimlens2(1) = 1
  dimlens2(2) = 1
  dimlens2(3) = 1
  CALL cdf_define(netcdfId, 'tilt', dimlens2, 'R4')
  CALL cdfSetatt(netcdfId, 'tilt', 'Tilt of the magnetic axis (degrees)')
  
  dimlens2(1) = 6
  dimlens2(2) = 1
  dimlens2(3) = 1
  CALL cdf_define(netcdfId, 'wTsyg', dimlens2, 'R4')
  
  CALL CPU_TIME(start_time)    ! start timer, measured in seconds
  
  IF (.NOT. ALLOCATED(alatMod)) ALLOCATE(alatMod(nlat, nlon), STAT = ierr)
  
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
  IF (.NOT. ALLOCATED(xSWMF)) ALLOCATE(xSWMF(nthe,nRadSWMF,nlonSWMF-1), stat = ierr)
  IF (.NOT. ALLOCATED(ySWMF)) ALLOCATE(ySWMF(nthe,nRadSWMF,nlonSWMF-1), stat = ierr)
  IF (.NOT. ALLOCATED(zSWMF)) ALLOCATE(zSWMF(nthe,nRadSWMF,nlonSWMF-1), stat = ierr)
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
     END DO
  END DO
  ! Do NATGRID interpolation (scattered) to find x, y, z for each such surface
  DO iPoint = 2, nPoints
     !PRINT*, 'cf: iPoint = ', iPoint
     CALL Interpolation_natgrid_ShellBuild(REAL(xSWMF(iPoint,:,:), DP), alatMod(:,:), rlonval(1:nlon), iPoint, nlat, nlon, &
          xTsyg(iPoint,1:nlat,2:nlon+1))
     CALL Interpolation_natgrid_ShellBuild(REAL(ySWMF(iPoint,:,:), DP), alatMod(:,:), rlonval(1:nlon), iPoint, nlat, nlon, &
          yTsyg(iPoint,1:nlat,2:nlon+1))
     CALL Interpolation_natgrid_ShellBuild(REAL(zSWMF(iPoint,:,:), DP), alatMod(:,:), rlonval(1:nlon), iPoint, nlat, nlon, &
          zTsyg(iPoint,1:nlat,2:nlon+1))
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
  
  IF (iDebug /= 0) THEN
     dimlens2(1) = nRadSWMF
     dimlens2(2) = nLonSWMF-1
     dimlens2(3) = 1
     CALL  cdf_define (netcdfId, 'xEqSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'yEqSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'nEqSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'pEqSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'LonIonSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'LatIonSWMF', dimlens2, 'R4')
     dimlens2(1) = nlat
     dimlens2(2) = nlon
     dimlens2(3) = 1
     CALL  cdf_define (netcdfId, 'LatGrid', dimlens2, 'R4')
     dimlens2(1) = nlon
     dimlens2(2) = 1
     dimlens2(3) = 1
     CALL  cdf_define (netcdfId, 'LonGrid', dimlens2, 'R4')
     dimlens2(1) = ntheOld
     dimlens2(2) = npsi
     dimlens2(3) = nzeta + 1
     CALL  cdf_define (netcdfId, 'xTsygOld', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'yTsygOld', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'zTsygOld', dimlens2, 'R4')    
     dimlens2(1) = ntheOld
     dimlens2(2) = nRadSWMF
     dimlens2(3) = nLonSWMF-1
     CALL  cdf_define (netcdfId, 'xSWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'ySWMF', dimlens2, 'R4')
     CALL  cdf_define (netcdfId, 'zSWMF', dimlens2, 'R4')
     write(*,*)'x, y to file=', xSWMF(1:nPoints, 5,nLonSWMF-1), ySWMF(1:nPoints, 5,nLonSWMF-1)
  END IF
  
  dimlens2(1) = nthe
  dimlens2(2) = npsi
  dimlens2(3) = nzeta + 1
  CALL  cdf_define (netcdfId, 'xTsyg', dimlens2, 'R4')
  CALL  cdf_define (netcdfId, 'yTsyg', dimlens2, 'R4')
  CALL  cdf_define (netcdfId, 'zTsyg', dimlens2, 'R4')  
  
  ! Writing block - do not define further netcdf quantities
  IF (iDebug /= 0) THEN
     CALL  cdf_write (netcdfId, 'xEqSWMF', REAL(xEqSWMF,sp))
     CALL  cdf_write (netcdfId, 'yEqSWMF', REAL(yEqSWMF,sp))
     CALL  cdf_write (netcdfId, 'pEqSWMF', REAL(pEqSWMF,sp))
     CALL  cdf_write (netcdfId, 'nEqSWMF', REAL(nEqSWMF,sp))
     CALL  cdf_write (netcdfId, 'LonIonSWMF', REAL(LonSWMF-1,sp))
     CALL  cdf_write (netcdfId, 'LatIonSWMF', REAL(LatSWMF,sp))
     CALL  cdf_write (netcdfId, 'LatGrid', REAL(alatMod,sp))
     CALL  cdf_write (netcdfId, 'LonGrid', REAL(rlonval(1:nlon),sp))
  END IF
  
  !    CALL  cdf_write (netcdfId, 'statusBC', statusBound) 
  !    CALL  cdf_write (netcdfId, 'pdyn', pdyn_char)
  !    CALL  cdf_write (netcdfId, 'byimf', byimf_char)
  !    CALL  cdf_write (netcdfId, 'bzimf', bzimf_char)
  !    CALL  cdf_write (netcdfId, 'DST', dst_char)
  CALL  cdf_write (netcdfId, 'xpsi_in', xpsiin_char)
  CALL  cdf_write (netcdfId, 'xpsi_out', xpsiout_char)
  CALL  cdf_write (netcdfId, 'rStart', r00_char)
  CALL  cdf_write (netcdfId, 'Constz', constz_char)
  CALL  cdf_write (netcdfId, 'Consttheta', consttheta_char)
  CALL  cdf_write (netcdfId, 'ntheta', nthe_char)
  CALL  cdf_write (netcdfId, 'npsi', npsi_char)
  CALL  cdf_write (netcdfId, 'nzeta', nzeta_char)
  CALL  cdf_write (netcdfId, 'tilt', REAL(0.0,dp))
  
  CALL  cdf_write (netcdfId, 'xTsyg', REAL(xTsygNew(:,1:npsi,:),sp))
  CALL  cdf_write (netcdfId, 'yTsyg', REAL(yTsygNew(:,1:npsi,:),sp))
  CALL  cdf_write (netcdfId, 'zTsyg', REAL(zTsygNew(:,1:npsi,:),sp))
  
  IF (iDebug /= 0) THEN
     CALL  cdf_write (netcdfId, 'xTsygOld', REAL(xTsyg,sp))
     CALL  cdf_write (netcdfId, 'yTsygOld', REAL(yTsyg,sp))
     CALL  cdf_write (netcdfId, 'zTsygOld', REAL(zTsyg,sp))
     
     CALL  cdf_write (netcdfId, 'xSWMF', REAL(xSWMF,sp))
     CALL  cdf_write (netcdfId, 'ySWMF', REAL(ySWMF,sp))
     CALL  cdf_write (netcdfId, 'zSWMF', REAL(zSWMF,sp))
  END IF
  
  ! close the cdf file 
  CALL cdf_close(netcdfId)
  
  WRITE(*, '(A, I3, 2X, F8.2)') 'hour, max(R) = ', iTime, SQRT(MAXVAL(xTsyg**2 + yTsyg**2));
  
  IF(ALLOCATED(xTsyg) .AND. ALLOCATED(yTsyg) .AND. ALLOCATED(zTsyg)) &
       DEALLOCATE(xTsyg, yTsyg, zTsyg)
  IF(ALLOCATED(xTsygNew) .AND. ALLOCATED(yTsygNew) .AND. ALLOCATED(zTsygNew))&
       DEALLOCATE(xTsygNew, yTsygNew, zTsygNew)
  IF(ALLOCATED(xSWMF) .AND. ALLOCATED(ySWMF) .AND. ALLOCATED(zSWMF))  &
       DEALLOCATE(xSWMF, ySWMF, zSWMF)
  
  !CALL CPU_TIME(end_time)          ! stop timer

END subroutine build_scb_init

LOGICAL FUNCTION isNaN(a) 
  USE nrtype, ONLY : SP, DP
  REAL(SP) :: a 
  IF (a.NE.a .OR. a*0._dp.NE.0._dp) THEN 
     isnan = .TRUE. 
  ELSE 
     isnan = .FALSE. 
  END IF
  RETURN 
END FUNCTION isnan
