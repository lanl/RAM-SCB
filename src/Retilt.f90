PROGRAM Retilt
  USE ezcdf
  USE nrtype, ONLY : SP, DP, pi

  IMPLICIT NONE

  !  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  INTEGER :: ierr, nid1, i, iDay, iHour, noTimes
  INTEGER, PARAMETER :: nthe=51, npsi=35, nzeta=61, NPA=72
  INTEGER, DIMENSION(3)   :: dimlens 
  INTEGER, DIMENSION(3)   :: dimlens2
  CHARACTER*500 :: filenameIn, filenameOut
  REAL(SP) :: xOld(nthe,npsi,nzeta), yOld(nthe,npsi,nzeta), zOld(nthe,npsi,nzeta), &
       bXOld(nthe,npsi,nzeta), bYOld(nthe,npsi,nzeta), bZOld(nthe,npsi,nzeta), b(nthe,npsi,nzeta), &
       xNew(nthe,npsi,nzeta), yNew(nthe,npsi,nzeta), zNew(nthe,npsi,nzeta), &
       bXNew(nthe,npsi,nzeta), bYNew(nthe,npsi,nzeta), bZNew(nthe,npsi,nzeta), &
       xSM(nthe,npsi,nzeta), ySM(nthe,npsi,nzeta), zSM(nthe,npsi,nzeta), &
       bXSM(nthe,npsi,nzeta), bYSM(nthe,npsi,nzeta), bZSM(nthe,npsi,nzeta)
  REAL(SP) :: I_value_Yue(npsi,nzeta,NPA)
  REAL(SP) :: PA(NPA)
  REAL(SP) :: pdyn, byimf, bzimf, dst, tiltOld, tiltNew
  REAL(SP) :: rHour
  REAL(SP) :: wTsyg(6)
  REAL(SP), DIMENSION(npsi,nzeta) :: xEq, yEq, BEq, JacobianEq
  REAL(SP), DIMENSION(nthe,npsi) :: xMid, zMid, BMid, JacobianMid, xNoo, zNoo, BNoo, JacobianNoo
  REAL(SP) :: chi(nthe)
  REAL(SP), DIMENSION(npsi) :: alpha, latitude
  REAL(SP) :: beta(nzeta)
  REAL(DP) :: AAA(10), SPS, CPS, BBB(23)
  REAL(DP) :: G(66), H(66), REC(66)
  CHARACTER*5 :: HourDec
  COMMON /GEOPACK1/ AAA,SPS,CPS,BBB
  COMMON /GEOPACK2/ G,H,REC

  iDay = 7
  noTimes = 12*24

  DO i = 0, noTimes-1
     rHour = REAL(i)/12._sp
     iHour = FLOOR(rHour)
     WRITE(HourDec, '(F5.2)') rHour
     if (iHour < 10) HourDec(1:1) = '0'
     PRINT*, 'iHour, rHour, HourDec = ', iHour, rHour, HourDec

     filenameIn = 'Day07/OldB/mag_field_'//HourDec//'.cdf'

     CALL cdf_open(nid1, TRIM(filenameIn), 'r') 
     CALL cdf_read(nid1, 'pdyn', pdyn)
     CALL cdf_read (nid1, 'byimf', byimf)
     CALL cdf_read (nid1, 'bzimf', bzimf)
     CALL cdf_read (nid1, 'dst', dst)
     CALL cdf_read (nid1, 'tilt', tiltOld)
     CALL cdf_read(nid1, 'wTsyg', wTsyg)
     CALL cdf_read(nid1, 'xEq', xEq)
     CALL cdf_read(nid1, 'yEq', yEq)
     CALL cdf_read(nid1, 'BEq', BEq)
     CALL cdf_read(nid1, 'JacobianEq', JacobianEq)
     CALL cdf_read(nid1, 'xMid', xMid)
     CALL cdf_read(nid1, 'zMid', zMid)
     CALL cdf_read(nid1, 'BMid', BMid)
     CALL cdf_read(nid1, 'xNoo', xNoo)
     CALL cdf_read(nid1, 'zNoo', zNoo)
     CALL cdf_read(nid1, 'BNoo', BNoo)
     CALL cdf_read(nid1, 'JacobianNoo', JacobianNoo)
     CALL cdf_read(nid1, 'chi', chi)
     CALL cdf_read(nid1, 'alpha', alpha)
     CALL cdf_read(nid1, 'latitude', latitude)
     CALL cdf_read(nid1, 'beta', beta)
     CALL  cdf_read (nid1, 'xFull', xOld)
     CALL  cdf_read (nid1, 'yFull', yOld)
     CALL  cdf_read (nid1, 'zFull', zOld)
     CALL  cdf_read (nid1, 'bX', bxOld)
     CALL  cdf_read (nid1, 'bY', byOld)
     CALL  cdf_read (nid1, 'bZ', bzOld)
     CALL cdf_read(nid1, 'bFull', b)
     CALL cdf_read(nid1, 'PA', PA)
     CALL cdf_read(nid1, 'IEq_Yue', I_value_Yue)
     CALL cdf_close(nid1)

     ! Untilt
     xSM = xOld*cosd(tiltOld) - zOld*sind(tiltOld)
     ySM = yOld
     zSM = xOld*sind(tiltOld) + zOld*cosd(tiltOld)
     bxSM = bxOld*cosd(tiltOld) - bzOld*sind(tiltOld)
     bySM = byOld
     bzSM = bxOld*sind(tiltOld) + bzOld*cosd(tiltOld)

     ! Get the right tilt

     CALL RECALC(2002,305+iDay-1,iHour,0,0) ! Only every hour change in tilt, as in the BCs every hour 
     tiltNew = 180./pi * BBB(4)

     print*, 'tiltOld, tiltNew = ', tiltOld, tiltNew

     ! Apply the proper tilt
     xNew = xSM*cosd(tiltNew) + zSM*sind(tiltNew)
     yNew = ySM
     zNew = -xSM*sind(tiltNew) + zSM*cosd(tiltNew)
     bxNew = bxSM*cosd(tiltNew) + bzSM*sind(tiltNew)
     byNew = bySM
     bzNew = -bxSM*sind(tiltNew) + bzSM*cosd(tiltNew)

     filenameOut = 'Day07/mag_field_'//HourDec//'.cdf'
     CALL cdf_open(nid1, TRIM(filenameOut), 'w') 

     dimlens(1) = 1
     dimlens(2) = 1
     dimlens(3) = 1
     CALL cdf_define(nid1, 'pdyn', dimlens, 'R4')
     CALL cdf_define(nid1, 'byimf', dimlens, 'R4')
     CALL cdf_define(nid1, 'bzimf', dimlens, 'R4')
     CALL cdf_define(nid1, 'dst', dimlens, 'R4')
     CALL cdf_define(nid1, 'tilt', dimlens, 'R4')
     dimlens(1) = 6
     CALL cdf_define(nid1, 'wTsyg', dimlens, 'R4')

     dimlens(1) = npsi
     dimlens(2) = nzeta
     dimlens(3) = 1
     CALL cdf_define(nid1, 'xEq', dimlens, 'R4')
     CALL cdf_define(nid1, 'yEq', dimlens, 'R4')
     CALL cdf_define(nid1, 'BEq', dimlens, 'R4')
     CALL cdf_define(nid1, 'JacobianEq', dimlens, 'R4')

     dimlens(1) = nthe
     dimlens(2) = npsi
     dimlens(3) = 1
     CALL cdf_define(nid1, 'xMid', dimlens, 'R4')
     CALL cdf_define(nid1, 'zMid', dimlens, 'R4')
     CALL cdf_define(nid1, 'BMid', dimlens, 'R4')
     CALL cdf_define(nid1, 'JacobianMid', dimlens, 'R4')
     CALL cdf_define(nid1, 'xNoo', dimlens, 'R4')
     CALL cdf_define(nid1, 'zNoo', dimlens, 'R4')
     CALL cdf_define(nid1, 'BNoo', dimlens, 'R4')
     CALL cdf_define(nid1, 'JacobianNoo', dimlens, 'R4')

     dimlens(1) = nthe
     dimlens(2) = 1
     dimlens(3) = 1
     CALL cdf_define(nid1, 'chi', dimlens, 'R4')
     CALL cdf_setatt(nid1, 'chi', 'Coordinate along the field - equidistant; dimension nchi')

     dimlens(1) = npsi
     dimlens(2) = 1
     dimlens(3) = 1
     CALL cdf_define(nid1, 'alpha', dimlens, 'R4')
     CALL cdf_setatt(nid1, 'alpha', 'One Euler potential - magnetic flux; dimension nalpha')
     CALL cdf_define(nid1, 'latitude', dimlens, 'R4')
     CALL cdf_setatt(nid1, 'latitude', 'Latitudes corresponding to Euler potential alpha; dimension nalpha')

     dimlens(1) = nzeta
     dimlens(2) = 1
     dimlens(3) = 1
     CALL cdf_define(nid1, 'beta', dimlens, 'R4')
     CALL cdf_setatt(nid1, 'beta', 'Other Euler potential - equiv. to azimuthal angle; periodic; dimension nbeta')

     dimlens(1) = nthe
     dimlens(2) = npsi
     dimlens(3) = nzeta

     CALL cdf_define(nid1, 'xFull', dimlens, 'R4')
     CALL cdf_define(nid1, 'yFull', dimlens, 'R4')
     CALL cdf_define(nid1, 'zFull', dimlens, 'R4')

     CALL cdf_define(nid1, 'bX', dimlens, 'R4')
     CALL cdf_define(nid1, 'bY', dimlens, 'R4')
     CALL cdf_define(nid1, 'bZ', dimlens, 'R4')
     CALL cdf_define(nid1, 'bFull', dimlens, 'R4')

     dimlens(1) = NPA
     dimlens(2) = 1
     dimlens(3) = 1
     CALL cdf_define(nid1, 'PA', dimlens, 'R4')

     dimlens(1) = npsi
     dimlens(2) = nzeta
     dimlens(3) = NPA
     CALL  cdf_define(nid1, 'IEq_Yue', dimlens, 'R4') 

     CALL cdf_write(nid1, 'chi', chi(1:nthe))
     CALL cdf_write(nid1, 'alpha', alpha(1:npsi))
     CALL cdf_write(nid1, 'latitude', latitude(1:npsi))
     CALL cdf_write(nid1, 'beta', beta(1:nzeta))

     CALL cdf_write(nid1, 'pdyn', pdyn)
     CALL cdf_write(nid1, 'byimf', byimf)
     CALL cdf_write(nid1, 'bzimf', bzimf)
     CALL cdf_write(nid1, 'dst', dst)
     CALL cdf_write(nid1, 'tilt', tiltNew)
     CALL cdf_write(nid1, 'wTsyg', wTsyg)

     CALL cdf_write(nid1, 'xEq', xEq)  
     CALL cdf_write(nid1, 'yEq', yEq)  
     CALL cdf_write(nid1, 'BEq', bEq)  
     CALL cdf_write(nid1, 'JacobianEq', JacobianEq)  

     CALL cdf_write(nid1, 'xMid', xMid)
     CALL cdf_write(nid1, 'zMid', zMid)  
     CALL cdf_write(nid1, 'BMid', BMid)  
     CALL cdf_write(nid1, 'JacobianMid', JacobianMid)  
     CALL cdf_write(nid1, 'xNoo', xNoo)  
     CALL cdf_write(nid1, 'zNoo', zNoo)  
     CALL cdf_write(nid1, 'BNoo', BNoo)  
     CALL cdf_write(nid1, 'JacobianNoo', JacobianNoo)  

     CALL cdf_write(nid1, 'xFull', xNew)
     CALL cdf_write(nid1, 'yFull', yNew)
     CALL cdf_write(nid1, 'zFull', zNew)

     CALL cdf_write(nid1, 'bFull', b)  
     CALL cdf_write(nid1, 'bX', bXNew)  
     CALL cdf_write(nid1, 'bY', bYNew)  
     CALL cdf_write(nid1, 'bZ', bZNew)  

     CALL cdf_write(nid1, 'PA', PA(1:NPA))
     CALL  cdf_write(nid1, 'IEq_Yue', I_value_Yue(1:npsi,1:nzeta,1:NPA)) 

     CALL cdf_close(nid1)

  END DO

END PROGRAM Retilt
