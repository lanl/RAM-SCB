!******************************************************************************
SUBROUTINE pressure(entropy_local, vol_local, icount_local)
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

  USE nrtype 
  USE Module1
  USE ModRamMain, ONLY: PathSwmfOut, DoWritePres
  USE ModIoUnit,  ONLY: UNITTMP_

  IMPLICIT NONE

  INTERFACE
     FUNCTION SavGol(pres)
       USE nrtype, ONLY : dp
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: pres(:,:)
       REAL(DP), DIMENSION(SIZE(pres,1),SIZE(pres,2))  :: SavGol
     END FUNCTION SavGol
     FUNCTION SavGol7(pres)
       USE nrtype, ONLY : dp
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: pres(:,:)
       REAL(DP), DIMENSION(SIZE(pres,1),SIZE(pres,2))  :: SavGol7
     END FUNCTION SavGol7
  END INTERFACE

  INTERFACE Spline_2D_derivatives
     SUBROUTINE Spline_2D_derivs(x1, x2, f, derivsX1, derivsX2)
       USE EZspline_obj ! import the modules
       USE EZspline  
       use nrtype, ONLY : DP
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP 
       REAL(r8), DIMENSION(:), INTENT(IN) :: x1, x2 ! independent variables
       REAL(r8), DIMENSION(:,:), INTENT(IN) :: f, derivsX1, derivsX2  ! function values
     END SUBROUTINE Spline_2D_derivs
  END INTERFACE

  INTERFACE Spline_2D_point
     SUBROUTINE Spline_2D_point(x_1, x_2, f, x, y, func, ierrDomain)
       USE EZspline_obj ! import the modules
       USE EZspline  
       use nrtype, ONLY : DP
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP 
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1, x_2 ! independent variable
       REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
       REAL(r8), DIMENSION(:,:), INTENT(OUT) :: func   ! interpolated values 
       REAL(r8), DIMENSION(:,:), INTENT(IN)    :: x, y  ! grid of points for output
       INTEGER, INTENT(out) :: ierrDomain
     END SUBROUTINE Spline_2D_point
  END INTERFACE

  INTERFACE EZ_Spline_3D_derivatives
     SUBROUTINE Spline_coord_derivs(x_1, x_2, x_3, f_input, derivsX1, derivsX2, derivsX3)
       USE Module1
       USE EZspline_obj ! import the modules
       USE EZspline  
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_2
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_3 ! independent variables
       REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f_input
       REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: derivsX1, derivsX2, derivsX3 
     END SUBROUTINE Spline_coord_derivs
  END INTERFACE

  INTERFACE smooth
     SUBROUTINE smooft(Y, N, PTS)
       USE nrtype
       REAL(DP) :: Y(1024)
       INTEGER :: N, PTS
     END SUBROUTINE smooft
  END INTERFACE

  INTERFACE locate
     FUNCTION locate(xx,x)
       USE nrtype
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xx
       REAL(DP), INTENT(IN) :: x
       INTEGER  :: locate
     END FUNCTION locate
  END INTERFACE

  INTERFACE polinomial_interpolation  ! (and extrapolation)
     SUBROUTINE polint(xa,ya,x,y,dy)
       USE nrtype; USE nrutil, ONLY : assert_eq, iminloc, nrerror
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa, ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y, dy
     END SUBROUTINE polint
  END INTERFACE

  REAL(DP), INTENT(IN) :: entropy_local(:,:), vol_local(:,:)
  INTEGER, INTENT(IN) :: icount_local
  INTEGER                       :: i, iloopOut, ierflg, j, j1, k1, jSKBoundary, k, ierr, ierrDom, idealerr, m1, mstate, n1
  REAL(DP), EXTERNAL            :: pressureRad, anisEq

  REAL(DP) :: press(npsi, nzeta+1), dPresdRho(npsi, nzeta+1), dPresdZeta(npsi, nzeta+1), &
       xEq(npsi, nzeta+1), yEq(npsi, nzeta+1),  &
       aratio(npsi, nzeta+1), aratioOld(npsi, nzeta+1), &
       aLiemohn(npsi, nzeta+1), dSqPresdRhoSq(npsi,nzeta+1), dSqPresdZetaSq(npsi,nzeta+1), &
       dSqPresdRhodZeta(npsi,nzeta+1), pperEq(npsi,nzeta+1), pparEq(npsi,nzeta+1), &
       pperEqOld(npsi,nzeta+1), pparEqOld(npsi,nzeta+1), &
       radGridEq(npsi, nzeta), angleGridEq(npsi,nzeta)
  REAL(DP) :: radius, angle, bEqSq, aN, pperN, pparN
  REAL(DP) :: distance(npsi), distance2derivs(npsi)
  REAL(DP) :: yyp, factorChange, &
       gParam, pEq, ratioB, rBI, bd, colatitudeMid, dipoleFactorMid(nthe,npsi), &
       colatitudeNoo, dipoleFactorNoo(nthe,npsi), pressureNonL
  REAL(DP), ALLOCATABLE :: coeffLsq(:), coeffLsqGeotail(:)
  INTEGER :: j_local, k_local, iplx, iChange, numberCoeffLsqGeo, numberCoeffLsqDMSP
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: dBsqdRho, dBsqdZeta
  REAL(DP), ALLOCATABLE:: BigBracketPsi(:,:,:), &
       BigBracketAlpha(:,:,:), dBBdRho(:,:,:), dBBdZeta(:,:,:), dummy1(:,:,:), dummy2(:,:,:)
  REAL(DP) :: rCenter, rr1, rr2, thangle, zangle, pMin, pMax, deltaCS, deltaPhi, deltaPhi2, pressSK, &
       delta1, delta2, x1, x2, pUp, pDown, pUp2, pDown2, coeffUp, coeffDown
  REAL(DP) :: press1D(npsi), pressMid(npsi)
  REAL(DP), PARAMETER :: tiny = 1e-6_dp
  REAL(DP) :: dydummy
  REAL(SP) :: dummyLine(10)
  INTEGER, PARAMETER :: nXRoe = 17, nYRoe = 14, nEnergRoe = 12, nPARoe = 18
  INTEGER, PARAMETER :: nXRoeGeo = 8 ! Index of first Roeder radius > 6.6 RE (or less, if overlapping is chosen) !
  !C (more if GEO data to be more efficient in determining the fit)
  !C INTEGER, PARAMETER :: nXRaw = 121, nYRaw = 49 ! For DMSP runs
  REAL(DP) :: xRaw(nXRaw,nYRaw), YRaw(nXRaw,nYRaw),pressProtonPerRaw(nXRaw,nYRaw), pressProtonParRaw(nXRaw,nYRaw), &
       pressOxygenPerRaw(nXRaw,nYRaw), pressOxygenParRaw(nXRaw,nYRaw), pressHeliumPerRaw(nXRaw,nYRaw), &
       pressHeliumParRaw(nXRaw,nYRaw), pressPerRaw(nXRaw,nYRaw), pressParRaw(nXRaw,nYRaw), &
       pressEleParRaw(nXRaw,nYRaw), pressElePerRaw(nXRaw,nYRaw), &     !Vania
       radRaw(nXRaw), azimRaw(nYRaw), ratioRaw(nXRaw,nYRaw), &
       radRoe(nXRoe), azimRoe(nYRoe), energRoe(0:nEnergRoe), PARoe(nPARoe), fluxRoe(nXRoe, nYRoe, nEnergRoe, 18), &
       pressProtonPerRoe(nXRoe, nYRoe), pressProtonParRoe(nXRoe, nYRoe), &
       pressPerRoe(nXRoe, nYRoe), pressParRoe(nXRoe, nYRoe), ratioRoe(nXRoe, nYRoe)
  !C INTEGER, PARAMETER :: nXRawExt = 171, nYRawExt = 49 ! For DMSP runs
  INTEGER, PARAMETER :: nXRawExt = nXRaw+2, nYRawExt = nAzimRAM
  ! To extend the domain to 10 RE, w/ the standard RAM resolution we need extra 14 radial cells, thus nXRawExt = nXRaw+14 = 33
  ! To extend the domain to 11 RE, w/ the standard RAM resolution we need extra 18 radial cells, thus nXRawExt = nXRaw+18 = 37
  ! (The domain over which pressure is defined has to include the domain delimited by magnetic boundaries)

  REAL(DP) :: xRawExt(nXRawExt,nYRawExt), YRawExt(nXRawExt,nYRawExt),pressPerRawExt(nXRawExt,nYRawExt), & 
       pressParRawExt(nXRawExt,nYRawExt), radRawExt(nXRawExt), azimRawExt(nYRawExt), ratioRawExt(nXRawExt,nYRawExt)
  REAL(DP), PARAMETER :: l0 = 50._dp
  CHARACTER*93 :: firstLine, secondLine
  CHARACTER*200 :: header
  INTEGER, PARAMETER :: ISLIM = nXRaw*nYRaw, NUMXOUT = npsi, NUMYOUT = nzeta-1
  !  INTEGER, PARAMETER :: IDIM = 2*NUMXOUT*NUMYOUT
  !  REAL(dp) :: X_neighbor(ISLIM), Y_neighbor(ISLIM), Z_neighbor(ISLIM), indexPsi(npsi,nzeta), &
  !       indexAlpha(npsi,nzeta)

  ! REAL(DP) :: XI(NUMXOUT), YI(NUMYOUT)
  ! REAL     ::        XP(NUMXOUT), YP(NUMYOUT), ZP(NUMXOUT,NUMYOUT)
  ! INTEGER :: IWORK(IDIM)
  INTEGER :: ier, iCount_neighbor, iDomain
  REAL(DP) :: w1, w2, w3, w4, w5, w6, w7, w8, w9
  REAL(DP), PARAMETER :: Rweight = 0.1_dp, gammaEnt = 5./3.

  INTEGER, PARAMETER :: lwrk = 50000, lwrk1 = 500000, lwrk2 = 500000
  INTEGER :: iopt(3), iopt1, ider(2), nu, nv
  INTEGER, SAVE :: nxout, nyout, nxoutPer, nxoutPar, nyoutPer, nyoutPar
  REAL :: pressPerRawRowExt(nXRawExt*nYRawExt), pressParRawRowExt(nXRawExt*nYRawExt)

  REAL :: wrk(lwrk), wrk1(lwrk1), wrk2(lwrk2)
  REAL :: fpResids, fpResidsPer, fpResidsPar
  REAL :: smoothFactor, smoothFactorPer, smoothFactorPar  
  INTEGER, PARAMETER :: kwrk = 50000, kwrk1 = 5000, nuest = nXRaw+7, nvest = nYRaw+7 !C nuest = nXRawExt + 7, nvest = nYRawExt+7
  INTEGER, PARAMETER :: kx = 3, ky = 3  ! Must be 3 for polar, can vary for surfit
  !C  INTEGER, PARAMETER :: nxest = 24, nyest = 12, nmax = MAX(nxest, nyest)
  INTEGER, PARAMETER :: nxest = 15, nyest = 15, nmax = MAX(nxest, nyest) 
  ! For Roeder expansion, not wise to go for larger nx, ny as it might force an unnatural spline
  REAL :: coeff((nuest-4)*(nvest-4))
  REAL, SAVE :: coeff1((nxest-kx-1)*(nyest-ky-1)), coeff2((nxest-kx-1)*(nyest-ky-1))
  INTEGER :: iwrk(kwrk), iwrk1(kwrk1)
  REAL :: tu(nuest), tv(nvest)
  REAL, SAVE :: tx(nxest), ty(nyest), txPer(nxest), txPar(nxest), tyPer(nxest), tyPar(nxest)
  REAL :: radCenter, radDisk, radMin, radMax, phiBeg, phiEnd, z0, val 
  REAL :: t, tout, ydriv, epsFit, epsdriv, deltaDev

  INTEGER, PARAMETER :: number = 5929, mlat_range = 121, mlon_range = 49
  !C INTEGER, PARAMETER :: mlat_range_Y = 48, mlon_range_Y = 25, number_Y = mlat_range_Y*mlon_range_Y ! Prev. case
  !C INTEGER, PARAMETER :: mlat_range_Y = 95, mlon_range_Y = 49, number_Y = mlat_range_Y*mlon_range_Y ! Higher res
  INTEGER, PARAMETER :: mlat_range_Y = 47, mlon_range_Y = 49, number_Y = mlat_range_Y*mlon_range_Y
  INTEGER :: indexLatMax
  INTEGER, PARAMETER :: nuestY = mlat_range_Y+7, nvestY = mlon_range_Y+7
  REAL(DP) :: radRawY(mlat_range_Y), azimRawY(mlon_range_Y), &
       pressPerRawY(mlat_range_Y, mlon_range_Y), pressParRawY(mlat_range_Y, mlon_range_Y), &
       ro(mlat_range_Y, mlon_range_Y), mlto(mlat_range_Y, mlon_range_Y)
  REAL :: tuY(nuestY), tvY(nvestY), coeffY((nuestY-4)*(nvestY-4))
  REAL,  ALLOCATABLE :: r(:), phi(:), xSp(:), ySp(:), pValue(:), pValuePer(:), pValuePar(:), weight(:), &
       weightPer(:), weightPar(:), u(:), v(:)
  REAL(DP), ALLOCATABLE :: xGeo(:), yGeo(:), radGeo(:), angleGeo(:), factorPerGeo(:), factorParGeo(:), pPerGeo(:), pParGeo(:)
  REAL(DP) :: factorPer, factorPar
  INTEGER  :: iloop, m, n, ierralloc, mlon2, mlat2, nc, nGeo
  REAL(DP) :: f_sum_sq
  REAL(DP), ALLOCATABLE :: f_vec(:), lat(:), latDummy(:), lon(:), mlt(:), pres(:), pressureIono(:,:)
  REAL(DP)  :: p1_main, p2_main, rad, lon2(100), lat2(100) , yTemp, wTemp, pTemp, presMax, dataTemp(4), &
       pressAt10, coeffIncrease
  EXTERNAL :: fdriv
  INTEGER, EXTERNAL :: is_nan ! C function
  REAL, EXTERNAL :: radFunc, evapol
  REAL(DP), EXTERNAL :: pRoeRad
  REAL(DP) :: xAr(48), yAr(48)
  REAL(DP) :: xSWMF(48,48), ySWMF(48,48), pressSWMF(48,48), rhoSWMF(48,48) 
  CHARACTER*4 :: ST3
  ! LOGICAL :: isnand ! intrinsic for PGF

  ! PRINT*, 'Beginning of pressure call; icount_local = ', icount_local

  iCountPressureCall = iCountPressureCall + 1 ! global variable, counts how many times pressure is called

  DO  j = 1,npsi
     DO  i = 1,nthe
        DO  k = 2,nzeta
           xx(i,j,k) = SQRT(x(i,j,k)**2 + y(i,j,k)**2)
           yy(i,j,k) = x(i,j,k) / xx(i,j,k)
           yyp = yy(i,j,k)
           yy(i,j,k) = ACOS(yyp)
        END DO
        !  avoid k = 2 for y near zero so that phi can be negative, but small
        IF(y(i,j,2) < 0.0) yy(i,j,2)=-yy(i,j,2)
        DO  k=3,nzeta
           IF(y(i,j,k) < 0.0) yy(i,j,k)=twopi_d-yy(i,j,k)
        END DO
        DO  k=2,nzeta
           dela(k) = alfa(i,j,k) - yy(i,j,k)
        END DO
        dela(1) = dela(nzeta)
        dela(nzetp) = dela(2)
        xx(i,j,1) = xx(i,j,nzeta)
        xx(i,j,nzetp) = xx(i,j,2)
        yy(i,j,1) = alfa(i,j,1) - dela(1)
        yy(i,j,nzetp) = alfa(i,j,nzetp) - dela(nzetp)
     END DO
  END DO


1001 Isotropy_choice:  IF (isotropy == 1) THEN    ! isotropic case
     CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), press(:,1:nzeta), &
          dPresdRho(:,1:nzeta), dPresdZeta(:,1:nzeta))
     press(:,nzetp) = press(:,2)
     press(:,1) = press(:,nzeta)
     dPresdRho(:,nzetp) = dPresdRho(:,2)
     dPresdRho(:,1) = dPresdRho(:,nzeta)
     dPresdZeta(:,nzetp) = dPresdZeta(:,2)
     dPresdZeta(:,1) = dPresdZeta(:,nzeta)

     CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), dPresdRho(:,1:nzeta), &
          dSqPresdRhoSq(:,1:nzeta), dSqPresdRhodZeta(:,1:nzeta))
     CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), dPresdZeta(:,1:nzeta), &
          dSqPresdRhodZeta(:,1:nzeta), dSqPresdZetaSq(:,1:nzeta))
     dSqPresdRhoSq(:,nzetp) = dSqPresdRhoSq(:,2)
     dSqPresdRhoSq(:,1) = dSqPresdRhoSq(:,nzeta)
     dSqPresdZetaSq(:,nzetp) = dSqPresdZetaSq(:,2)
     dSqPresdZetaSq(:,1) = dSqPresdZetaSq(:,nzeta)

     DO i = 1, nthe
        DO j = 1, npsi
           dpdPsi(i,j,1:nzetp) = 1._dp / f(j) * dPresdRho(j,1:nzetp)
           IF (iOuterMethod == 2) dSqPdPsiSq(i,j,1:nzetp) = 1._dp / f(j)**2 * dSqPresdRhoSq(j,1:nzetp)
        END DO
        DO k = 1, nzetp
           dpdAlpha(i,1:npsi, k) = dPresdZeta(1:npsi,k) / fzet(k)
           IF (iOuterMethod == 2) dSqPdAlphaSq(i,1:npsi,k) = dSqPresdZetaSq(1:npsi,k) / fzet(k)**2
        END DO
        pressure3D(i,1:npsi,1:nzetp) = press(1:npsi, 1:nzetp)
     END DO


  ELSE    ! Anisotropic pressure case
     IF (iPressureChoice == 5) THEN ! Calculation using RAM pressures
        OPEN(UNITTMP_, file=TRIM(fileNamePressure), status = 'OLD')
        READ(UNITTMP_, '(A)') firstLine
        READ(UNITTMP_, '(A)') secondLine
        ! print*, 'rank, 1st line: ', rank, firstLine; call flush(6)
        ! print*, 'rank, 2nd line:', rank, secondLine; call flush(6)
        ! print*, 'HERE'; call flush(6)
        DO j1 = 1, nXRaw
           DO k1 = 1, nYRaw
              READ(UNITTMP_, *) radRaw(j1), azimRaw(k1), pressProtonPerRaw(j1,k1), pressProtonParRaw(j1,k1), &
                   pressOxygenPerRaw(j1,k1), pressOxygenParRaw(j1,k1), pressHeliumPerRaw(j1,k1), &
                   pressHeliumParRaw(j1,k1)
           END DO
        END DO

        azimRaw = azimRaw * 360./24 * pi_d / 180._dp ! In radians

        pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw) ! from keV/cm^3 to nPa
        pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw) ! from keV/cm^3 to nPa

        ratioRaw = pressOxygenPerRaw/pressProtonPerRaw

        CLOSE(UNITTMP_)

        radRawExt(1:nXRaw) = radRaw(1:nXRaw)
        DO j1 = nXRaw+1, nXRawExt
           radRawExt(j1) = radRaw(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw(nXRaw)-radRaw(1))/(REAL(nXRaw-1, DP))
        END DO

        azimRawExt(1:nYRawExt) = azimRaw(1:nYRaw) ! nYRaw = nYRawExt

        pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
        pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
        ratioRawExt(1:nXRaw,:) = ratioRaw(1:nXRaw,:)
        
        DO k1 = 1, nYRawExt
           DO j1 = nXRaw+1, nXRawExt
              ! Alternatively, extrapolate the pressure assuming SK dependence in regions where we don't know it instead of f(R) extrapolation; 
              ! or, decrease the order of the polynomial extrapolation
              pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
                   (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
              pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
                   (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53)) 
              ratioRawExt(j1,k1) = pressOxygenPerRaw(nXRaw,k1)/pressProtonPerRaw(nXRaw,k1)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), pressPerRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), pressPerRawExt(j1,k1),dydummy)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), pressParRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), pressParRawExt(j1,k1),dydummy)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), ratioRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), ratioRawExt(j1,k1),dydummy)
           END DO
        END DO

        !  anisotropic pressure functions
        DO k = 2, nzeta
           DO j = 1, npsi

              radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
              angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d

              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
                   angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)

              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
                   angle = - ASIN(y(nThetaEquator,j,k) / radius)

              radGrid(j,k) = radius
              angleGrid(j,k) = angle

           END DO
        END DO

        IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing (possibly multiple) for the pressure
           pressPerRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressPerRawExt(1:nXRawExt,1:nYRawExt))
           pressParRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressParRawExt(1:nXRawExt,1:nYRawExt))
        END IF

        !Cubic spline interpolation
        CALL Spline_2D_point(radRawExt**2, azimRawExt, pressPerRawExt, &
             radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain) 
        CALL Spline_2D_point(radRawExt**2, azimRawExt, pressParRawExt, &
             radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain) 
        CALL Spline_2D_point(radRawExt**2, azimRawExt, ratioRawExt, &
             radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), ratioEq(1:npsi,2:nzeta), iDomain) 
        ratioEq(1:npsi,1) = ratioEq(1:npsi,nzeta)
        IF (iDomain > 0) THEN
           PRINT*, 'Stop; problem with pressure domain; iDomain = ', iDomain
           STOP
        END IF

        ! Sometimes the interpolation can give very small negative values very 
        ! near the Earth; inconsequential
        WHERE(pperEq < 0.0) pperEq = 1e-2_dp
        WHERE(pparEq < 0.0) pparEq = 1e-2_dp
        DO k = 2, nzeta
           DO j = 1, npsi
              IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2 RE from Earth
                 pperEq(j,k) = pressPerRaw(1,1)
                 pparEq(j,k) = pressParRaw(1,1)
              END IF
           END DO
        END DO

       pperEq(:,nzeta+1) = pperEq(:,2)
       pparEq(:,nzeta+1) = pparEq(:,2)
       pperEq(:,1) = pperEq(:,nzeta)
       pparEq(:,1) = pparEq(:,nzeta)

       pperEq = pperEq/pnormal
       pparEq = pparEq/pnormal


        DO k = 1, nzeta
           DO j = 1, npsi

              pEq = (2.*pperEq(j,k) + pparEq(j,k)) / 3._dp 
              !C if (rank == 0) print*, 'j, k, pEq = ', j, k, pEq; call flush(6)
              aratio(j,k) = pperEq(j,k) / pparEq(j,k) - 1._dp  
              aLiemohn(j,k) = - aratio(j,k) / (aratio(j,k)+1_dp)

              DO i = 1, nthe
                 ratioB = bf(nThetaEquator,j,k) / bf(i,j,k)
                 !      if (rank == 0 .and. i == 2) print*, 'j, k, ratioB = ', j, k, ratioB; call flush(6)

                 IF (iLossCone == 2) THEN
                    ! New reference values (Liemohn)
                    rBI = MAX(bf(1,j,k)/bf(i,j,k), 1._dp+1.E-9_dp)  ! Must be larger than 1, i.e. the field at "Earth" higher than last field value 
                    !C print*, 'i, j, k, rBI = ', i, j, k, rBI
                    pparN = pparEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    pperN = pperEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    aN = pparN/pperN - 1._dp
                    ! if (rank == 0 .and. i == nthe) then
                    ! write(*, '(A, 1X, 2I3, 1X, 4F9.3)') 'j, k, bf, ratioB, ppar, pper = ', j, k, bf(i,j,k), ratioB, pparN, pperN; CALL FLUSH(6)
                    ! stop
                    ! end if
                    ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * SQRT((rBI-1._dp)/(rBI-ratioB)) * &
                         (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
                    pper(i,j,k) = ppar(i,j,k) / (1._dp + aN*ratioB)
                    !C if (rank == 0 .and. i == 1) print*, 'j, k, rBI, ppar, pper = ', j, k, rBI, ppar(i,j,k), pper(i,j,k); call flush(6)
                 END IF
                 IF (iLossCone == 1) THEN
                    gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                    ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                    pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
                 END IF
                 sigma(i,j,k) = 1._dp + (pper(i,j,k)-ppar(i,j,k)) / bsq(i,j,k)
                 tau(i,j,k) = 1._dp - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * &
                      pper(i,j,k)/ppar(i,j,k)
              END DO
              press(j,k) = pEq
           END DO
        END DO

     ELSEIF (iPressureChoice == 6) THEN ! RAM pressures & Roeder model extension; only protons for now in the 
        ! extension formula, but apply the formula to the total (H+, He++, O+) pressures

        !  radius, angle in flux coordinates
        DO k = 2, nzeta
           DO j = 1, npsi
              radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
              angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
                   angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
                   angle = - ASIN(y(nThetaEquator,j,k) / radius)
              radGrid(j,k) = radius
              angleGrid(j,k) = angle
           END DO
        END DO

        OPEN(UNITTMP_, file=TRIM(fileNamePressure), status = 'OLD')
        READ(UNITTMP_, '(A)') firstLine
        READ(UNITTMP_, '(A)') secondLine
        DO j1 = 1, nXRaw
           DO k1 = 1, nYRaw
              READ(UNITTMP_, *) radRaw(j1), azimRaw(k1), pressProtonPerRaw(j1,k1), pressProtonParRaw(j1,k1), &
                   pressOxygenPerRaw(j1,k1), pressOxygenParRaw(j1,k1), pressHeliumPerRaw(j1,k1), &
                   pressHeliumParRaw(j1,k1), pressElePerRaw(j1,k1), pressEleParRaw(j1,k1)
           END DO
        END DO
        azimRaw = azimRaw * 360./24 * pi_d / 180._dp
        pressPerRaw = 0.16_dp * (pressProtonPerRaw + 1.*pressOxygenPerRaw + 1.*pressHeliumPerRaw + 1.*pressElePerRaw) 
! from keV/cm^3 to nPa
        pressParRaw = 0.16_dp * (pressProtonParRaw + 1.*pressOxygenParRaw + 1.*pressHeliumParRaw + 1.*pressEleParRaw) 
! from keV/cm^3 to nPa
        ! ratioRaw = pressOxygenPerRaw/pressProtonPerRaw
        CLOSE(UNITTMP_)

        radRawExt(1:nXRaw) = radRaw(1:nXRaw)
        DO j1 = nXRaw+1, nXRawExt
           radRawExt(j1) = radRaw(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw(nXRaw)-radRaw(1))/(REAL(nXRaw-1, DP))
           !   print*, 'pressure: j1, radRawExt = ', j1, radRawExt(j1)
        END DO

        !        DO j1 = 1, nXRawExt
        !           print*, 'pressure: j1, radRawExt(j1) = ', j1, radRawExt(j1)
        !        END DO

        !        print*, 'pressure: nYRaw, nYRawExt = ', nYRaw, nYRawExt
        azimRawExt(1:nYRawExt) = azimRaw(1:nYRaw) ! nYRaw = nYRawExt
        !        print*, 'pressure: azimRawExt(1) = ', azimRawExt(1)
        pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
        pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
        ! ratioRawExt(1:nXRaw,:) = ratioRaw(1:nXRaw,:)

        ! Try at 6.75 RE to have the same pressure as @ 6.5 RE (until RAM boundary issues are solved)
        !C        pressPerRawExt(nXRaw+1,:) = pressPerRawExt(nXRaw,:)
        !C        pressParRawExt(nXRaw+1,:) = pressParRawExt(nXRaw,:)

        DO k1 = 1, nYRawExt
           DO j1 = nXRaw+1, nXRawExt
              ! Alternatively, extrapolate the pressure assuming SK dependence in regions where we don't know it instead of f(R) extrapolation; 
              ! or, decrease the order of the polynomial extrapolation
              !C      pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
              !C           (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
              !C      pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
              !C           (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53)) 

              pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
              pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) *  pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
              !            PRINT*, 'pressure: j1, k1, presspar = ', j1, k1, pressParRawExt(j1,k1)
              ! ratioRawExt(j1,k1) = pressOxygenPerRaw(nXRaw,k1)/pressProtonPerRaw(nXRaw,k1)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), pressPerRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), pressPerRawExt(j1,k1),dydummy)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), pressParRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), pressParRawExt(j1,k1),dydummy)
              !C CALL polint(radRawExt(nXRaw-2:nXRaw), ratioRaw(nXRaw-2:nXRaw,k1), radRawExt(j1), ratioRawExt(j1,k1),dydummy)
           END DO
        END DO

        IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing (possibly multiple) for the pressure
           pressPerRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressPerRawExt(1:nXRawExt,1:nYRawExt))
           pressParRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressParRawExt(1:nXRawExt,1:nYRawExt))
        END IF


        ! Piecewise cubic spline interpolation; alternative - put all points scattered and do surfit
204     CALL Spline_2D_point(radRawExt, azimRawExt, pressPerRawExt, &
             radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain) 
        CALL Spline_2D_point(radRawExt, azimRawExt, pressParRawExt, &
             radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain) 
        ! CALL Spline_2D_point(radRawExt**2, azimRawExt, ratioRawExt, &
        !      radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), ratioEq(1:npsi,2:nzeta), iDomain) 
        ! ratioEq(1:npsi,1) = ratioEq(1:npsi,nzeta)
        IF (iDomain > 0) THEN
           PRINT*, 'Stop; problem with pressure domain; iDomain = ', iDomain
           STOP
        END IF

        ! Now outside geosynchronous orbit - if the DIERCKX routines are used
        !C        DO k = 2, nzeta
        !C          DO j = 1, npsi 
        !C              if (radGrid(j,k) > 11._dp) then ! Outside geosynchronous orbit
        !C                 ! evapol (with DIERCKX polar) not implemented yet
        !C                 val = evapol(txPer,nxoutPer,tyPer,nyoutPer,coeff1,radFunc, REAL(-(radGrid(j,k)-radCenter)*cos(angleGrid(j,k))), &
        !C                      REAL(-(radGrid(j,k)-radCenter)*sin(angleGrid(j,k))))
        !C               !C  CALL bispev(txPer(1:nxoutPer),nxoutPer,tyPer(1:nyoutPer),nyoutPer,coeff1(1:(nxoutPer-kx-1)*(nyoutPer-ky-1)),&
        !C               !C       kx,ky,REAL(radGrid(j,k)),1,REAL(angleGrid(j,k)),1, val,wrk2,lwrk2,iwrk1,kwrk1,ier)
        !C                 pperEq(j,k) = (REAL(val, DP))

        !C                 val = evapol(txPar,nxoutPar,tyPar,nyoutPar,coeff2,radFunc,REAL(-(radGrid(j,k)-radCenter)*cos(angleGrid(j,k))), &
        !C                      REAL(-(radGrid(j,k)-radCenter)*sin(angleGrid(j,k))))
        !C                 !C CALL bispev(txPar(1:nxoutPar),nxoutPar,tyPar(1:nyoutPar),nyoutPar,coeff2(1:(nxoutPar-kx-1)*(nyoutPar-ky-1)),&
        !C                 !C     kx,ky,REAL(radGrid(j,k)),1,REAL(angleGrid(j,k)),1, val,wrk2,lwrk2,iwrk1,kwrk1,ier)
        !C                 pparEq(j,k) = (REAL(val, DP))
        !C                 ! print*, 'j, k, pperEq, pparEq = ', j, k, pperEq(j,k), pparEq(j,k)
        !C              end if
        !C              
        !C           END DO
        !C        END DO


        !C if (rank == 0) print*, 'bispev - ier = ', ier
        IF (ALLOCATED(pValue)) DEALLOCATE(pValue, STAT=ierr)

        pperEq = pperEq / pnormal
        pparEq = pparEq / pnormal

        ! Sometimes the interpolation can give very small negative values very 
        ! near the Earth; inconsequential
        WHERE(pperEq < 0.0) pperEq = MINVAL(pressPerRaw) ! 1e-1_dp/pnormal
        WHERE(pparEq < 0.0) pparEq = MINVAL(pressParRaw) ! 1e-1_dp/pnormal
        DO k = 1, nzeta
           DO j = 1, npsi
              IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2RE from Earth
                 pperEq(j,k) = pperEq(4,k) ! MINVAL(pressPerRaw)/pnormal ! 1e-1_dp/pnormal
                 pparEq(j,k) = pparEq(4,k) ! MINVAL(pressParRaw)/pnormal ! 1e-1_dp/pnormal
              END IF
           END DO
        END DO
        pperEq(:,nzeta+1) = pperEq(:,2)
        pparEq(:,nzeta+1) = pparEq(:,2)

        ! Do smoothing on the raw pressure, not here
        ! IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing
        !   pperEq(1:npsi,2:nzeta+1) = SavGol(pperEq(1:npsi,2:nzeta+1))
        !   pparEq(1:npsi,2:nzeta+1) = SavGol(pparEq(1:npsi,2:nzeta+1))
        ! END IF

        pperEq(:,1) = pperEq(:,nzeta)
        pparEq(:,1) = pparEq(:,nzeta)

        DO k = 1, nzeta
           DO j = 1, npsi

       !       IF (pperEq(j,k) < 0. .OR. pparEq(j,k) < 0.) THEN
       !          PRINT*, 'pressure: j, k, pper, ppar = ', j, k, pperEq(j,k), pparEq(j,k)
       !       END IF

              pEq = (2.*pperEq(j,k) + pparEq(j,k)) / 3._dp 
              !    if (rank == 0) print*, 'j, k, pEq = ', j, k, pEq; call flush(6)
              aratio(j,k) = pperEq(j,k) / pparEq(j,k) - 1._dp  
              aLiemohn(j,k) = - aratio(j,k) / (aratio(j,k)+1_dp)

              DO i = 1, nthe
                 ratioB = bf(nThetaEquator,j,k) / bf(i,j,k)
                 !      if (rank == 0 .and. i == 2) print*, 'j, k, ratioB = ', j, k, ratioB; call flush(6)
     !C            IF (is_nan(DBLE(ratioB)) /= 0) THEN
     !C               PRINT*, 'pressure: ratioB j, k = ', j, k
     !C               STOP
     !C            END IF
                 IF (iLossCone == 2) THEN
                    ! New reference values (Liemohn)
                    rBI = MAX(bf(1,j,k)/bf(i,j,k), 1._dp+1.E-9_dp)  ! Must be larger than 1, i.e. the field at "Earth" higher than last field value 
                    !C print*, 'i, j, k, rBI = ', i, j, k, rBI
                    pparN = pparEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    pperN = pperEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    aN = pparN/pperN - 1._dp
                    ! if (rank == 0 .and. i == nthe) then
                    ! write(*, '(A, 1X, 2I3, 1X, 4F9.3)') 'j, k, bf, ratioB, ppar, pper = ', j, k, bf(i,j,k), ratioB, pparN, pperN; CALL FLUSH(6)
                    ! stop
                    ! end if
                    ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * SQRT((rBI-1._dp)/(rBI-ratioB)) * &
                         (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
                    pper(i,j,k) = ppar(i,j,k) / (1._dp + aN*ratioB)
                    !C if (rank == 0 .and. i == 1) print*, 'j, k, rBI, ppar, pper = ', j, k, rBI, ppar(i,j,k), pper(i,j,k); call flush(6)
                 END IF
                 IF (iLossCone == 1) THEN
                    !     print*, 'pressure: j, k, aratio(j,k), ratioB = ', j, k, aratio(j,k), ratioB
                    gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                    ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                    pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
   !C                 IF (is_nan(SQRT(DBLE(gParam))) /= 0) THEN
   !C                    PRINT*, 'pressure: j, k = ', j, k
   !C                    STOP
   !C                 END IF
                 END IF
                 sigma(i,j,k) = 1._dp + (pper(i,j,k)-ppar(i,j,k)) / bsq(i,j,k)
                 tau(i,j,k) = 1._dp - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * &
                      pper(i,j,k)/ppar(i,j,k)

        !C         IF (ppar(i,j,k) < 0. .OR. pper(i,j,k) < 0.) THEN
        !C         PRINT*, 'pressure: i, j, k, pper, ppar = ', i, j, k, pper(i,j,k), ppar(i,j,k)
        !C      END IF

              END DO
              press(j,k) = pEq
              !C           if (is_nan(press(j,k))) THEN
              !C              print*, 'pressure: j, k = ', j, k
              !C              STOP
              !C           end if
           END DO
        END DO


     ELSEIF (iPressureChoice == 11) THEN  ! SWMF pressures outside 6.5 RE; RAM pressures inside

        !  radius, angle in flux coordinates
        DO k = 2, nzeta
           DO j = 1, npsi
              radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
              angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
                   angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
                   angle = - ASIN(y(nThetaEquator,j,k) / radius)
              radGrid(j,k) = radius
              angleGrid(j,k) = angle
           END DO
        END DO

        OPEN(UNITTMP_, file=TRIM(fileNamePressure), status = 'OLD')
        READ(UNITTMP_, '(A)') firstLine
        READ(UNITTMP_, '(A)') secondLine
        DO j1 = 1, nXRaw
           DO k1 = 1, nYRaw
              READ(UNITTMP_, *) radRaw(j1), azimRaw(k1), pressProtonPerRaw(j1,k1), pressProtonParRaw(j1,k1), &
                   pressOxygenPerRaw(j1,k1), pressOxygenParRaw(j1,k1), pressHeliumPerRaw(j1,k1), &
                   pressHeliumParRaw(j1,k1)
           END DO
        END DO
        azimRaw = azimRaw * 360./24 * pi_d / 180._dp
        pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw) ! from keV/cm^3 to nPa
        pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw)  ! from keV/cm^3 to nPa
        ! ratioRaw = pressOxygenPerRaw/pressProtonPerRaw
        CLOSE(UNITTMP_)

        radRawExt(1:nXRaw) = radRaw(1:nXRaw)
        DO j1 = nXRaw+1, nXRawExt
           radRawExt(j1) = radRaw(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw(nXRaw)-radRaw(1))/(REAL(nXRaw-1, DP))
        END DO

        azimRawExt(1:nYRawExt) = azimRaw(1:nYRaw) ! nYRaw = nYRawExt
        pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
        pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)

        DO j = 1, 48
           xAr(j) = -11.75 + 23.5*(j-1)/47.
           yAr(j) = -11.75 + 23.5*(j-1)/47.
        END DO

        WRITE(ST3,'(I4)') iSWMF
        IF(iSWMF.LE.999) ST3(1:1)='0'
        IF(iSWMF.LE.99) ST3(2:2)='0'
        IF(iSWMF.LE.9) ST3(3:3)='0'
        !C        print*, 'iSWMF = ', iSWMF

        OPEN(UNITTMP_, file = trim(PathSwmfOut)//'/pressure_geo_'// &
             ST3//'.dat', form = 'unformatted', status = 'old', &
             action = 'read')
        ! PRINT*, 'PRESSURE: calling ', 'pressure_geo_'//ST3//'.dat'
        READ(UNITTMP_) xSWMF
        READ(UNITTMP_) ySWMF
        READ(UNITTMP_) pressSWMF
        READ(UNITTMP_) rhoSWMF
        CLOSE(UNITTMP_)

        ! Verify that the grid is correct
        DO k = 1, 48
           DO j = 1, 48
              !C            print*, 'PRESSURE:',  j, k, xSWMF(j,k), ySWMF(j,k), pressSWMF(j,k), rhoSWMF(j,k)
              !C            if (xSWMF(j,k)**2 + ySWMF(j,k)**2 < 4.) pressSWMF(j,k) = minval(pressSWMF)
              IF (ABS(xSWMF(j,k)-xAr(j))>1e-3_dp .OR. ABS(ySWMF(j,k)-yAr(k))>1e-3_dp) STOP 'Problem.'
           END DO
        END DO

        DO k = 1, nYRawExt
           DO j = 1, nXRawExt
              xRawExt(j,k) = - radRawExt(j) * COS(azimRawExt(k))
              yRawExt(j,k) = - radRawExt(j) * SIN(azimRawExt(k))
           END DO
        END DO

        !C    IF (iSWMF == 0) THEN ! Initially, SWMF pressure in whole domain
        !C       CALL Spline_2D_point(xAr, yAr, pressSWMF, xRawExt(:,:), &
        !C            yRawExt(:,:), pressPerRawExt(:,:), iDomain) 
        !C       pressParRawExt = pressPerRawExt ! Isotropy
        !C           ELSE 
        CALL Spline_2D_point(xAr, yAr, pressSWMF, xRawExt(nXRaw+1:nXRawExt,:), &
             yRawExt(nXRaw+1:nXRawExt,:), pressPerRawExt(nXRaw+1:nXRawExt,:), iDomain) 
        CALL Spline_2D_point(xAr, yAr, pressSWMF, xRawExt(nXRaw+1:nXRawExt,:), &
             yRawExt(nXRaw+1:nXRawExt,:), pressParRawExt(nXRaw+1:nXRawExt,:), iDomain) 
        !C    END IF

        ! For better smoothing at 6.5 - 6.75 RE (if boundary is there), extend the RAM pressure (same or extrapolation)
        pressPerRawExt(nXRaw+1,:) = 0.5*(pressPerRawExt(nXRaw,:)+pressPerRawExt(nXRaw+1,:))
        pressParRawExt(nXRaw+1,:) = 0.5*(pressParRawExt(nXRaw,:)+pressParRawExt(nXRaw+1,:))

        IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing for the pressure
           pressPerRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressPerRawExt(1:nXRawExt,1:nYRawExt))
           pressParRawExt(1:nXRawExt,1:nYRawExt) = SavGol7(pressParRawExt(1:nXRawExt,1:nYRawExt))
        END IF

        CALL Spline_2D_point(radRawExt, azimRawExt, pressPerRawExt, &
             radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain) 

        !C    IF (iSWMF == 0) THEN
        !C       pparEq = pperEq
        !C    ELSE
        CALL Spline_2D_point(radRawExt, azimRawExt, pressParRawExt, &
             radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain) 
        !C    END IF

        pperEq = pperEq / pnormal
        pparEq = pparEq / pnormal

        ! Sometimes the interpolation can give very small negative values very 
        ! near the Earth; inconsequential
        WHERE(pperEq < 0.0) pperEq = 1e-1_dp/pnormal
        WHERE(pparEq < 0.0) pparEq = 1e-1_dp/pnormal
        DO k = 1, nzeta
           DO j = 1, npsi
              IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2 RE from Earth
                 pperEq(j,k) = pressPerRaw(1,1)/pnormal ! 1e-1_dp/pnormal
                 pparEq(j,k) = pressParRaw(1,1)/pnormal ! 1e-1_dp/pnormal
              END IF
           END DO
        END DO

        pperEq(:,nzeta+1) = pperEq(:,2)
        pparEq(:,nzeta+1) = pparEq(:,2)
        pperEq(:,1) = pperEq(:,nzeta)
        pparEq(:,1) = pparEq(:,nzeta)

        DO k = 1, nzeta
           DO j = 1, npsi
              pEq = (2.*pperEq(j,k) + pparEq(j,k)) / 3._dp 
              !C if (rank == 0) print*, 'j, k, pEq = ', j, k, pEq; call flush(6)
              aratio(j,k) = pperEq(j,k) / pparEq(j,k) - 1._dp  
              aLiemohn(j,k) = - aratio(j,k) / (aratio(j,k)+1_dp)

              DO i = 1, nthe
                 ratioB = bf(nThetaEquator,j,k) / bf(i,j,k)
                 IF (iLossCone == 2) THEN
                    ! New reference values (Liemohn)
                    rBI = MAX(bf(1,j,k)/bf(i,j,k), 1._dp+1.E-9_dp)  ! Must be larger than 1, i.e. the field at "Earth" higher than last field value 
                    pparN = pparEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    pperN = pperEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    aN = pparN/pperN - 1._dp
                    ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * SQRT((rBI-1._dp)/(rBI-ratioB)) * &
                         (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
                    pper(i,j,k) = ppar(i,j,k) / (1._dp + aN*ratioB)
                 END IF
                 IF (iLossCone == 1) THEN
                    gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                    ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                    pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
                 END IF
                 sigma(i,j,k) = 1._dp + (pper(i,j,k)-ppar(i,j,k)) / bsq(i,j,k)
                 tau(i,j,k) = 1._dp - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * &
                      pper(i,j,k)/ppar(i,j,k)
                 !             print*, 'pressure: k, j, i, ppar, pper: ', k, j, i, ppar(i,j,k), pper(i,j,k)
                 !       IF (isnand(pper(i,j,k))) STOP 'pressure: NaN.'
              END DO
              press(j,k) = pEq
           END DO
        END DO



     ELSEIF (iPressureChoice == 7) THEN   ! Roeder model (not coupling with RAM)

        DO n1 = 1, nPARoe
           PARoe(n1) = 180._dp * REAL(n1-1, dp) / (nPARoe-1)
        END DO

        pressProtonPerRoe = 0._dp
        pressProtonParRoe = 0._dp

        OPEN(UNITTMP_, file='polar_cammice_comp_protons.txt', status = 'OLD')
        DO j = 1, 7
           READ(UNITTMP_, '(a)') header
        END DO
        !C print*, 'rank, header: ', rank, header; call flush(6)
        energRoe(0) = 0._dp
        DO j1 = 2, nXRoe
           DO k1 = 2, nYRoe-1
              DO m1 = 1, nEnergRoe
                 READ(UNITTMP_, '(f5.2, f6.1, f7.1, 18(f15.1))') radRoe(j1), azimRoe(k1), energRoe(m1), fluxRoe(j1, k1, m1, 1:18)
                 ! Problem with Roeder's model
                 !C   if (radRoe(j1) == 8.25 .and. azimRoe(k1) == 11.0) WRITE(*, '(f5.2, f6.1, f7.1, 18(f15.1))') radRoe(j1), azimRoe(k1), energRoe(m1), fluxRoe(j1, k1, m1, 1:18)
              END DO
              ! Obtain pressures
              DO m1 = 1, nEnergRoe
                 DO n1 = 1, nPARoe
                    ! In nPa:
                    pressProtonPerRoe(j1,k1) = pressProtonPerRoe(j1,k1) + 0.16_dp * pi_d * 4.56_dp * 1E-8_dp * &
                         SQRT(energRoe(m1)) * (energRoe(m1)-energRoe(m1-1)) * &
                         (SIN(pi_d/180._dp*PARoe(n1)))**3 * 10._dp*pi_d/180._dp * fluxRoe(j1, k1, m1, n1)
                    pressProtonParRoe(j1,k1) = pressProtonParRoe(j1,k1) + 0.16_dp * twopi_d * 4.56_dp *&
                         1E-8_dp * SQRT(energRoe(m1)) * (energRoe(m1)-energRoe(m1-1)) * &
                         (COS(pi_d/180._dp*PARoe(n1)))**2 * SIN(pi_d/180._dp*PARoe(n1)) * 10._dp*pi_d/180._dp &
                         * fluxRoe(j1, k1, m1, n1)
                 END DO
              END DO
              !C              WRITE(*, '(A, 2F7.2, 3F11.3)') 'rad, MLT, pper, ppar, ratioSK = ', radRoe(j1), azimRoe(k1), &
              !C                   pressProtonPerRoe(j1, k1), pressProtonParRoe(j1,k1), &
              !C                   pressProtonPerRoe(j1,k1) / (89.*EXP(-0.59*radRoe(j1)) + 8.9*radRoe(j1)**(-1.53))
           END DO
        END DO
        CLOSE(UNITTMP_)

        radRoe(1) = 2._dp
        azimRoe(1) = 0._dp
        azimRoe(nYRoe) = 24._dp

        pressProtonPerRoe(1,:) = pressProtonPerRoe(2,:) ! Values at R=2 RE
        pressProtonParRoe(1,:) = pressProtonParRoe(2,:)
        pressProtonPerRoe(:,1) = 0.5_dp*(pressProtonPerRoe(:,2)+pressProtonPerRoe(:,nYRoe-1))
        pressProtonParRoe(:,1) = 0.5_dp*(pressProtonParRoe(:,2)+pressProtonParRoe(:,nYRoe-1))
        pressProtonPerRoe(:,nYRoe) = 0.5_dp*(pressProtonPerRoe(:,2)+pressProtonPerRoe(:,nYRoe-1))
        pressProtonParRoe(:,nYRoe) = 0.5_dp*(pressProtonParRoe(:,2)+pressProtonParRoe(:,nYRoe-1)) 

        azimRoe = azimRoe * pi_d / 12._dp ! From MLT to radians
        !        ratioRoe = pressProtonPeRaw/pressProtonPerRaw

        DO k = 2, nzeta
           DO j = 1, npsi
              radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
              angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
                   angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)
              IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
                   angle = - ASIN(y(nThetaEquator,j,k) / radius)
              radGrid(j,k) = radius
              angleGrid(j,k) = angle
              ! if (j == 1) print*, 'pressure: j, k, radGrid, angleGrid = ', j, k, radGrid(j,k), angleGrid(j,k)
           END DO
        END DO

        CALL Spline_2D_point(radRoe**2, azimRoe, pressProtonPerRoe, &
             radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain) 
        CALL Spline_2D_point(radRoe**2, azimRoe, pressProtonParRoe, &
             radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain) 
        !    CALL Spline_2D_point(radRoe**2, azimRoe, ratioRoe, &
        !         radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), ratioEq(1:npsi,2:nzeta), iDomain) 
        !    ratioEq(1:npsi,1) = ratioEq(1:npsi,nzeta)
        IF (iDomain > 0) THEN
           PRINT*, 'Stop; problem with pressure domain; iDomain = ', iDomain
           STOP
        END IF


        pperEq(:,nzeta+1) = pperEq(:,2)
        pparEq(:,nzeta+1) = pparEq(:,2)
        pperEq(:,1) = pperEq(:,nzeta)
        pparEq(:,1) = pparEq(:,nzeta)

        ! From nPa to normalized code pressure units
        pperEq = pperEq / pnormal
        pparEq = pparEq / pnormal

        ! Sometimes the interpolation can give very small negative values very 
        ! near the Earth; inconsequential (field there is dipolar)
        WHERE(pperEq < 0.0) pperEq = 1e-2_dp/pnormal
        WHERE(pparEq < 0.0) pparEq = 1e-2_dp/pnormal

        DO k = 1, nzeta
           DO j = 1, npsi
              pEq = (2.*pperEq(j,k) + pparEq(j,k)) / 3._dp 
              ! if (k == nZetaMidnight) PRINT*, 'j, k, radGrid(j,k), pEq = ', j, k, radGrid(j,k), pEq
              aratio(j,k) = pperEq(j,k) / pparEq(j,k) - 1._dp  
              aLiemohn(j,k) = - aratio(j,k) / (aratio(j,k)+1_dp)
              DO i = 1, nthe
                 ratioB = bf(nThetaEquator,j,k) / bf(i,j,k)
                 IF (iLossCone == 2) THEN
                    rBI = MAX(bf(1,j,k)/bf(i,j,k), 1._dp+1.E-9_dp)  ! Must be larger than 1, i.e. the field 
                    ! at "Earth" higher than last field value 
                    pparN = pparEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    pperN = pperEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                    aN = pparN/pperN - 1._dp
                    ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * SQRT((rBI-1._dp)/(rBI-ratioB)) * &
                         (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
                    pper(i,j,k) = ppar(i,j,k) / (1._dp + aN*ratioB)
                 END IF
                 IF (iLossCone == 1) THEN
                    gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                    ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                    pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
                 END IF
                 sigma(i,j,k) = 1._dp + (pper(i,j,k)-ppar(i,j,k)) / bsq(i,j,k)
                 tau(i,j,k) = 1._dp - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * &
                      pper(i,j,k)/ppar(i,j,k)
              END DO
           END DO
        END DO

     ELSE     
        STOP 'PROBLEM in pressure.f90'
     END IF


     ! Block for reducing anisotropy
     IF (iReduceAnisotropy == 1) THEN 
        DO k = 1, nzeta
           DO j = 1, npsi
              Mirror_unstable:  IF (tau(nThetaEquator,j,k) < 0._dp) THEN
                 pEq = press(j,k)
                 bEqSq = bsq(nThetaEquator,j,k)
                 pperEq(j,k) = 1./6. * (3.*pEq - bEqSq + SQRT(bEqSq**2 + 12.*bEqSq*pEq + 9.*pEq**2))
                 pparEq(j,k) = 3.*pEq - 2.*pperEq(j,k)
                 aratio(j,k) = pperEq(j,k)/pparEq(j,k) - 1.
                 aLiemohn(j,k) = - aratio(j,k) / (aratio(j,k)+1._dp)
                 ! print*, 'pressure: j, k, aratio = ', j, k, aratio(j,k)
                 DO i = 1, nthe
                    ratioB = bf(nThetaEquator,j,k) / bf(i,j,k)
                    IF (iLossCone == 2) THEN
                       rBI = MAX(bf(1,j,k)/bf(i,j,k), 1._dp+1.E-9_dp)  
                       pparN = pparEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                       pperN = pperEq(j,k) * (1._dp - (ratioB+aLiemohn(j,k)*ratioB)/(rBI+aLiemohn(j,k)*ratioB))
                       aN = pparN/pperN - 1._dp
                       ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * &
                            SQRT((rBI-1._dp)/(rBI-ratioB)) * (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
                       pper(i,j,k) = ppar(i,j,k) / (1._dp + aN*ratioB)
                    END IF
                    IF (iLossCone == 1) THEN
                       gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                       ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                       pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
                    END IF
                    ! New sigmas and taus
                    sigma(i,j,k) = 1.0 + (pper(i,j,k)-ppar(i,j,k))/bsq(i,j,k)
                    tau(i,j,k) = 1. - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * pper(i,j,k)/ppar(i,j,k)
                 END DO
              END IF Mirror_unstable
           END DO
        END DO
     END IF

     IF (iOuterMethod == 2) THEN
        ALLOCATE(BigBracketPsi(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(BigBracketAlpha(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(dBBdRho(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(dBBdZeta(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(dummy1(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(dummy2(nthe,npsi,nzeta), stat = ierr)
     END IF

     CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, pper(1:nthe,1:npsi,1:nzeta), &
          dPperdTheta, dPperdRho, dPperdZeta)
     CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, bsq(1:nthe,1:npsi,1:nzeta), &
          dBsqdTheta, dBsqdRho, dBsqdZeta)

     IF (iOuterMethod == 2)  THEN ! Newton method
        DO k = 1, nzeta
           DO j = 1, npsi
              DO i = 1, nthe
                 BigBracketAlpha(i,j,k) = (-1./sigma(i,j,k) * dPperdAlpha(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j)**2 * fzet(k) * (gradRhoSq(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradRhoGradZeta(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k))*0.5*dBsqdTheta(i,j,k)) - &
                      (1. - sigma(i,j,k)) / sigma(i,j,k) * 0.5 * dBsqdAlpha(i,j,k))
                 BigBracketPsi(i,j,k) = (1./sigma(i,j,k) * dPperdPsi(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j) * fzet(k)**2 * (gradRhoGradZeta(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradZetaSq(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k)) * 0.5_dp * dBsqdTheta(i,j,k)) + &
                      (1.-sigma(i,j,k)) / sigma(i,j,k) * 0.5_dp * dBsqdPsi(i,j,k))
              END DO
           END DO
        END DO
        CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, BigBracketAlpha(1:nthe,1:npsi,1:nzeta), &
             dummy1, dummy2, dBBdZeta)    ! for Newton method
        CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, BigBracketPsi(1:nthe,1:npsi,1:nzeta), &
             dummy1, dBBdRho, dummy2)    ! for Newton method
     END IF

     DO j = 1, npsi
        dPperdPsi(:,j,:) = 1./f(j) * dPperdRho(:,j,:)
        IF (iOuterMethod == 2) dBBdPsi(:,j,:) = dBBdRho(:,j,:) / f(j)
        dBsqdPsi(:,j,:) = 1./f(j) * dBsqdRho(:,j,:)
     END DO

     DO k = 1, nzeta
        dPperdAlpha(:,:,k) = 1. / fzet(k) * dPperdZeta(:,:,k)
        DO j = 1, npsi
           
        END DO
        IF (iOuterMethod == 2) dBBdAlpha(:,:,k) = dBBdZeta(:,:,k) / fzet(k)
        dBsqdAlpha(:,:,k) = 1. / fzet(k) * dBsqdZeta(:,:,k)
     END DO

     IF (iOuterMethod == 2) THEN
        IF(ALLOCATED(BigBracketPsi)) DEALLOCATE(BigBracketPsi, stat = idealerr)
        IF(ALLOCATED(BigBracketAlpha)) DEALLOCATE(BigBracketAlpha, stat = idealerr)
        IF(ALLOCATED(dBBdRho)) DEALLOCATE(dBBdRho, stat = idealerr)
        IF(ALLOCATED(dBBdZeta)) DEALLOCATE(dBBdZeta, stat = idealerr)
        IF(ALLOCATED(dummy1)) DEALLOCATE(dummy1, stat = idealerr)
        IF(ALLOCATED(dummy2)) DEALLOCATE(dummy2, stat = idealerr)
     END IF

  END IF Isotropy_choice

  DO j = 1, npsi
     DO i = 1, nthe
        colatitudeMid = ATAN2(xx(i,j,nZetaMidnight), z(i,j,nZetaMidnight))
        colatitudeNoo = ATAN2(xx(i,j,2), z(i,j,2))
        dipoleFactorMid(i,j) = SQRT(1. + 3. * (COS(colatitudeMid))**2) / (SIN(colatitudeMid))**6
        dipoleFactorNoo(i,j) = SQRT(1. + 3. * (COS(colatitudeNoo))**2) / (SIN(colatitudeMid))**6 
     END DO
  END DO

  IF (rank == 0) THEN
     Anisotropy: IF (isotropy /= 1) THEN ! VJ Write pressure data to disk, only once
        if (iConvGlobal==1.and.DoWritePres) CALL Write_pressure_anisotropic
     ELSE IF(isotropy == 1) THEN
        if (iConvGlobal==1.and.DoWritePres) CALL Write_pressure_isotropic
     END IF Anisotropy
  END IF

  !C STOP ! Temporary

  RETURN

END SUBROUTINE pressure

REAL FUNCTION radFunc(x)
  IMPLICIT NONE
  REAL :: x
  radFunc = 10. - 5. ! This will ensure the DIERCKX spline is defined up to R=10 RE

END FUNCTION radFunc

