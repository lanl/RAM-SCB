SUBROUTINE Euler_derivatives

  USE nrtype, ONLY : DP, pi_d, twopi_d
  USE nr, ONLY : qtrap
  USE Module1, ONLY : nthe, nThetaEquator, bf, psival, alphaVal, chiVal, &
       npsi, nzeta, nzetp, nZetaMidnight, x, y, z, bnormal, iHour, &
       iHourBegin, Hour, chiVal, rank, psi, alfa, alpha_Cyl, beta_Cyl, &
       dAlphadT, dBetadT
  USE Module_RAM
  USE ezcdf
  USE ModIoUnit, ONLY: UNITTMP_

  IMPLICIT NONE

  INTERFACE Interpolation_Delaunay_Euler_2D
     SUBROUTINE Interpolation_Delaunay_Euler_2D(r_local, azim_local, alpha_l, beta_l, alphaCyl_l, betaCyl_l)
       USE nrtype, ONLY : DP, pi_d
       USE Module1, ONLY : x, y, z, nThetaEquator
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: r_local(:), azim_local(:)
       REAL(DP), INTENT(IN) :: alpha_l(:), beta_l(:)
       REAL(DP), INTENT(OUT) :: alphaCyl_l(:,:), betaCyl_l(:,:)
     END SUBROUTINE Interpolation_Delaunay_Euler_2D
  END INTERFACE

  INTERFACE locate
     INTEGER FUNCTION locate(xx,x)
       USE nrtype
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xx
       REAL(DP), INTENT(IN) :: x
     END FUNCTION locate
  END INTERFACE

  INTEGER :: I, IC, iCount, ierr, IP1, j1, k1, NWK, ncdfId
  INTEGER, DIMENSION(3) :: dimlens 
  INTEGER, PARAMETER :: LIMIT = 10000, INF = 1
  REAL(DP) :: switch, end_time, yLowerLimit, chiHalf, bfHalf, &
       alphaT, alphaTM1, betaT, betaTM1, radius, angle
  REAL(DP) :: r0(npsi, nzeta)
  REAL(DP) :: MUBOUN, RWU, RE, HMIN, ME, DL1, DR, IR1, DPHI, CLC, dyDummy
  INTEGER, PARAMETER :: nXRaw_l = 61, nYRaw_l = 101  !Fine fixed grid needed for accurate inductive E computation
  REAL(DP) :: radRaw(nXRaw_l), azimRaw(nYRaw_l)
  CHARACTER*100 :: fileNameEuler

  dimlens = (/1, 1, 1/)

  DO j1 = 1, nXRaw_l
     radRaw(j1) = 2.05_dp + (6.45_dp-2.05_dp) * REAL(j1-1,DP)/REAL(nXRaw-1,DP)
  END DO
  DO k1 = 1, nYRaw_l
     azimRaw(k1) = 0.25_dp/REAL(nYRaw_l-1,DP) + twopi_d * REAL(k1-1,DP)/REAL(nYRaw_l-1,DP) ! To straddle the Euler potential grid
  END DO

  ! Interpolate data to POLAR coordinates 
  CALL Interpolation_Delaunay_Euler_2D(radRaw, azimRaw(1:nYRaw_l-1), psival(1:npsi),  &
       alphaVal(2:nzeta), alpha_Cyl(1:nXRaw_l,1:nYRaw_l-1, iHour), beta_Cyl(1:nXRaw_l,1:nYRaw_l-1, iHour))

  !C beta_Cyl(1:nXRaw, (nYRaw+1)/2, iHour) = beta_Cyl(1:nXRaw, 1, iHour) + twopi_d 

  ! Calculate dAlpha/dt and dBeta/dt (Euler potential time derivatives)
  DO k = 2, nzeta
     DO j = 1, npsi
        ! Find in which cylindrical sector the point (x, y) is, and find alpha, beta at T-1 there
        radius = SQRT(x(nThetaEquator,j,k)**2 + y(nThetaEquator,j,k)**2)
        angle = ASIN(y(nThetaEquator,j,k) / radius) !C + pi_d
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
             angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
        IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
             angle = pi_d - ASIN(y(nThetaEquator,j,k) / radius)
        IF ((x(nThetaEquator,j,k) > 0) .AND. (y(nThetaEquator,j,k) <= 0)) &
             angle = twopi_d + ASIN(y(nThetaEquator,j,k) / radius)
        j1 = locate(radRaw, radius)
        k1 = locate(azimRaw(1:nYRaw_l-2), angle)
        !C if (j == 1) print*, 'k, j1, k1, alphaVal, betaCyl = ', k, j1, k1, alphaVal(k), beta_Cyl(j1,k1, iHour); CALL FLUSH(6)
        IF (iHour > iHourBegin) THEN
           alphaTM1 = alpha_Cyl(j1, k1, iHour-1)
           betaTM1 = beta_Cyl(j1,k1, iHour-1)
        END IF
        ! Find alpha, beta there now at T
        alphaT = alpha_Cyl(j1, k1, iHour)
        betaT = beta_Cyl(j1,k1, iHour)

        IF (iHour > iHourBegin) THEN
           IF ((betaTM1 < 0.5_dp*pi_d .AND. betaT > pi_d) .OR. (betaTM1 > pi_d .AND. betaT < 0.5_dp*pi_d)) THEN  ! Fix periodicity issue
              PRINT*, 'j, k, betaTM1, betaT = ', j, k, betaTM1, betaT
              ! stop
           END IF
        END IF

        dAlphadT(nThetaEquator,j,k) = (alphaT - alphaTM1)/3600._dp
        dBetadT(nThetaEquator,j,k) = (betaT - betaTM1)/3600._dp
     END DO
  END DO

  dAlphadT(:,:,1) = dAlphadT(:,:,nzeta)
  dBetadT(:,:,1) = dBetadT(:,:,nzeta)
  dAlphadT(:,:,nzetp) = dAlphadT(:,:,2)
  dBetadT(:,:,nzetp) = dBetadT(:,:,2)

  ! READ(Hour, '(I3)') iHour
  PRINT*, 'Hour = ', Hour
  PRINT*, 'iHour = ', iHour
  filenameEuler = 'Euler_output_'//TRIM(Hour)//'.dat'
  PRINT*, 'filenameEuler = ', filenameEuler
  OPEN(UNITTMP_, file = TRIM(filenameEuler), action = 'write', status = 'replace')
  WRITE(UNITTMP_, 15) 'T= ', iHour
  WRITE(UNITTMP_, 25) '   Lsh        Azim     alpha(Lsh,MLT)  beta(Lsh,MLT)'
  DO i = 1, nXRaw_l
     DO j1 = 1, nYRaw_l
        WRITE(UNITTMP_, 21) radRaw(i), azimRaw(j1), alpha_Cyl(i,j1, iHour), beta_Cyl(i,j1, iHour)
     END DO
  END DO
  CLOSE(UNITTMP_)
  CALL FLUSH(UNITTMP_) 

15 FORMAT(A, I3)
2 FORMAT(F8.2, F10.1, 3X, I2, 6X, E12.4, 1X, E12.4, 1X, E12.4)
21 FORMAT(F9.3, F9.3, 6X, E12.4, 1X, E12.4)
25 FORMAT(A)


  RETURN

END SUBROUTINE Euler_derivatives
   
