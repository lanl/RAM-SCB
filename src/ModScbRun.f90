  MODULE ModScbRun
  ! Contains subroutines responsible for making the SCB calculations
  
    use ModScbVariables, ONLY: radEqMidNew, GradRhoSq, GradZetaSq, GradRhoGradZeta, &
                               GradRhoGradTheta, GradThetaGradZeta, GradPsiGradAlpha, &
                               dPPerdAlpha, dBBdAlpha, dBBdPsi, dPPerdPsi, dBsqdAlpha, &
                               dBsqdPsi, dBsqdTheta, bfInitial, pressure3D, bj, phij, &
                               ppar, pper, tau, sigma, dPdAlpha, dPdPsi, dSqPdAlphaSq, &
                               dSqPdPsiSq, DstDps, DstDpsInsideGeo, DstBiot, DstBiotInsideGeo, &
                               kmax, nisave, nitry, iteration, iConvGlobal, lconv, zetaVal
  
    implicit none
    save

    contains
  
!==============================================================================
  SUBROUTINE scb_run
    !!!! Module Variables
    use ModRamParams,    ONLY: boundary, electric, NameBoundMag
    use ModRamVariables, ONLY: KP  
    USE ModScbMain,      ONLY: damp, iSm, nrelax, numit, relax, thresh
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: alfa, alfaSav1, alfaSav2, psi, psiSav1, psiSav2, psiVal, &
                               alphaVal, blendAlpha, blendPsi, iAlphaMove, iPsiMove, &
                               decreaseConvAlpha, decreaseConvPsi, errorAlpha, errorPsi, &
                               diffmx, errorAlphaPrev, errorPsiPrev, x, y, z, sumb, &
                               sumdb, jacobian, xzero3, psiin, psiout, psitot, &
                               xpsiin, xpsiout, f, fp, fluxVolume, psiPrev, alfaPrev, &
                               constZ, fzet, fzetp, chiVal, chi, thetaVal, constTheta
    use ModScbParams,    ONLY: decreaseConvAlphaMin, decreaseConvPsiMin, blendMin, &
                               decreaseConvAlphaMax, decreaseConvPsiMax, blendMax, &
                               blendAlphaInit, blendPsiInit, MinSCBIterations, &
                               iAMR, isEnergDetailNeeded, isFBDetailNeeded, &
                               method, isotropy
    !!!! Module Subroutine/Functions
    USE ModScbInit,     ONLY: computeBandJacob_Initial
    USE ModScbEuler,    ONLY: alfges, psiges, mapalpha, mappsi, directAlpha, &
                              iterateAlpha, directPsi, iteratePsi, psiFunctions, &
                              InterpolatePsiR
    USE ModScbEquation, ONLY: newk, newj, metric, metrica ! LHS and RHS equations
    USE ModScbIO,       ONLY: computational_domain, Write_Convergence_Anisotropic, &
                              create_BField
    !!!! Share Modules
    USE ModIOUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, twopi_d
  
    IMPLICIT NONE
  
    INTEGER  :: iconv, nisave1, ierr, iCountEntropy
    INTEGER  :: i, j, k
    REAL(DP) :: sumdbconv, errorfirstalpha, diffmxfirstalpha, &
                errorfirstpsi, diffmxfirstpsi, dphi, phi
    REAL(DP) :: sumb1, sumdb1, diffmx1, xpsitot, xpl, psis
    REAL(DP) :: entropyFixed(npsi,nzeta)
    REAL(DP), ALLOCATABLE, SAVE :: xPrev(:,:,:), yPrev(:,:,:), zPrev(:,:,:)
 
    REAL(DP), PARAMETER :: pow = 1.0_dp, TINY = 1.E-15_dp

    decreaseConvAlpha = decreaseConvAlphaMin + (decreaseConvAlphaMax - decreaseConvAlphaMin) &
                        *(MIN(Kp,6._dp))**2/36.
    decreaseConvPsi   = decreaseConvPsiMin + (decreaseConvPsiMax - decreaseConvPsiMin) &
                        *(MIN(Kp,6._dp))**2/36.

!!!!! Compute the SCB computational domain
    call create_BField
    psiin   = -xzero3/xpsiin
    psiout  = -xzero3/xpsiout
    psitot  = psiout-psiin
    xpsitot = xpsiout - xpsiin
    DO j = 1, npsi
       psis = REAL(j-1, DP) / REAL(npsi-1, DP)
       xpl = xpsiin + xpsitot * psis**pow
       psival(j) = -xzero3 / xpl
       f(j) = (xzero3 / xpl**2) * xpsitot * pow * psis**(pow-1.)
       fp(j) = 0._dp ! If pow = 1
    END DO
    call psiges

    dphi  = twopi_d/REAL(nzeta-1, DP)
    DO k = 1, nzeta+1
       phi         = REAL(k-2, DP) * dphi
       alphaVal(k) = phi + constz*sin(phi) ! Concentration at midnight
       fzet(k)     = 1._dp + constz*COS(phi)
       fzetp(k)    = - constz * SIN(phi)
    END DO
    call alfges

    chiVal = (thetaVal + constTheta * SIN(2._dp*thetaVal))
    DO i = 1, nthe
       DO j = 1, npsi
          DO k = 2, nzeta
             chi(i,j,k) = chiVal(i)
          END DO
       END DO
    END DO
    chi(:,:,1) = chi(:,:,nzeta)
    chi(:,:,nzeta+1) = chi(:,:,2)


! Shouldn't derivative's be taken with respect to the actual distribution????
! If so, we need to change all thetaVal -> chiVal and zetaVal -> alphaVal
! As a first attempt let's try it here -ME
    thetaVal = chiVal
    zetaVal = alphaVal(1:nzeta)
!!!!!

    sumb1 = 0._dp
    sumdb1 = 0._dp
    diffmx1 = 0._dp
    entropyFixed = 0._dp
    fluxVolume = 0._dp
    blendAlpha = blendAlphaInit
    blendPsi   = blendPsiInit
  
    iteration = 1
    iCountEntropy = 1
    iconv = 0
    sumdbconv = 0.0_dp
    iConvGlobal = 0

    IF (iAMR == 1) THEN
       CALL findR
       CALL InterpolatePsiR
       CALL mappsi(0) ! Full mapping needed, psis changed
       CALL psifunctions
       CALL maptheta
    ENDIF

    Outeriters: DO

       call computeBandJacob_initial
       CALL metrica

       ! Define the right-hand side of the betaEuler equation
       CALL pressure

       ! Moved the convergence test to after the SCB computation completes
       IF (iteration == 1 .AND. isFBDetailNeeded == 1) THEN
          CALL Write_Convergence_Anisotropic('00')
          ! to see how far from equilibrium we are before computing
       END IF

       SCB_CALCULATION: IF (method /= 3) then
          psiPrev(:,:,:) = psi(:,:,:)
          alfaPrev(:,:,:) = alfa(:,:,:)

          CALL newk
  
          errorAlphaPrev = errorAlpha
  
          IF (iteration==1) THEN
             alfaSav1(:,:,:) = alfa(:,:,:)
             alfaSav2(:,:,:) = alfa(:,:,:) ! Before the calculation
          END IF

          IF (MOD(iteration,nrelax) == 0) blendAlpha = relax * blendAlpha

          blendAlpha = MAX(blendAlpha,blendMin)
          blendAlpha = MIN(blendAlpha,blendMax)
          SELECT CASE (method)
            CASE(1)
               CALL directAlpha
            CASE(2)
               CALL iterateAlpha
            END SELECT
  
          sumb1 = sumb
          sumdb1 = sumdb
          diffmx1 = diffmx
          nisave1 = nisave  
          errorAlpha = diffmx
  
          IF (iteration == 1) THEN
             errorfirstalpha = sumdb1
             diffmxfirstalpha = diffmx1
             errorAlphaPrev = errorAlpha
          END IF
          IF (sumdb1 < sumdbconv) sumdbconv = sumdb1
  
          iAlphaMove = 1
  
          IF (.NOT. ALLOCATED(xPrev)) ALLOCATE(xPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
          IF (.NOT. ALLOCATED(yPrev)) ALLOCATE(yPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
          IF (.NOT. ALLOCATED(zPrev)) ALLOCATE(zPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
          xPrev = x
          yPrev = y
          zPrev = z

          Move_points_in_alpha_theta: DO
             ! move zeta grid points along constant alphaEuler and theta lines
             CALL mapalpha(iSm)
             ! move theta grid points along constant alphaEuler and zeta lines
             CALL maptheta
             CALL metrica

             IF (MINVAL(jacobian) < 0._dp) THEN
                blendAlpha = damp * blendAlpha
                if (blendAlpha.lt.blendMin) then
                   call CON_stop('Failed to converge Alpha potential with minimum blend size')
                endif
                PRINT*, 'CE: Cycling alpha_theta pts, blendAlpha = ', blendAlpha
                alfaPrev(:,:,:) = alfa(:,:,:)
                IF (MOD(iteration,2)==1) THEN
                   alfa(:,:,:) = alfa(:,:,:)*blendAlpha + (1.-blendAlpha)*alfaSav1(:,:,:)
                ELSE
                   alfa(:,:,:) = alfa(:,:,:)*blendAlpha + (1.-blendAlpha)*alfaSav2(:,:,:)
                END IF
                ! Revert to previous point configuration
                x = xPrev
                y = yPrev
                z = zPrev
                CYCLE Move_points_in_alpha_theta
             END IF
             EXIT Move_points_in_alpha_theta
          END DO Move_points_in_alpha_theta
  
          IF (MOD(iteration,2)==1) THEN
             alfaSav1(:,:,:) = alfa(:,:,:)
          ELSE
             alfaSav2(:,:,:) = alfa(:,:,:)
          END IF

          IF (iAMR == 1) THEN
             CALL findR
             CALL InterpolatePsiR
             CALL mappsi(0)  ! Full mapping needed, changed psis
             CALL psifunctions
             CALL maptheta
          ENDIF

          !IF (isotropy == 1) CALL entropy(entropyFixed, fluxVolume, iCountEntropy)
  
          CALL metric
          IF (MINVAL(jacobian) < 0._dp) STOP 'CE: metric problem.'


          call computeBandJacob_initial
          CALL pressure
  
          !c  define the right-hand side of the alphaEuler equation
          CALL newj
  
          errorPsiPrev = errorPsi
          IF (iteration==1) THEN
             psiSav1(:,:,:) = psi(:,:,:)
             psiSav2(:,:,:) = psi(:,:,:) ! Before the calculation
          END IF

          IF (MOD(iteration,nrelax) == 0) blendPsi = relax * blendPsi

          blendPsi = MAX(blendPsi,blendMin)
          blendPsi = MIN(blendPsi,blendMax)
          SELECT CASE (method)
            CASE(1)
              CALL directPsi
            CASE(2)
              CALL iteratePsi
          END SELECT

          errorPsi = diffmx
          IF (iteration==1) errorPsiPrev = errorPsi
  
          iPsiMove = 1
  
          xPrev = x
          yPrev = y
          zPrev = z
          Move_points_in_psi_theta: DO
             ! Call metrica to find out if the jacobian is well behaved
             CALL mappsi(iSm)
             CALL maptheta
             CALL metric
             IF (MINVAL(jacobian) < 0._dp) THEN
                blendPsi = damp * blendPsi
                if (blendPsi.lt.blendMin) then
                   call CON_stop('Failed to converge Psi potential with minimum blend size')
                endif
                PRINT*, 'CE: Cycling psi_theta pts, blendPsi = ', blendPsi
                psiPrev(:,:,:) = psi(:,:,:)
                IF (MOD(iteration,2)==1) THEN
                   psi(:,:,:) = psi(:,:,:)*blendPsi + (1.-blendPsi)*psiSav1(:,:,:)
                ELSE
                   psi(:,:,:) = psi(:,:,:)*blendPsi + (1.-blendPsi)*psiSav2(:,:,:)
                END IF
                ! Revert to previous point configuration
                x = xPrev
                y = yPrev
                z = zPrev
                CYCLE Move_points_in_psi_theta
             END IF
             EXIT Move_points_in_psi_theta
          END DO Move_points_in_psi_theta
  
          IF (MOD(iteration,2)==1) THEN
             psiSav1(:,:,:) = psi(:,:,:)
          ELSE
             psiSav2(:,:,:) = psi(:,:,:)
          END IF

          ! PRINT*, 'CE: blendPsi = ', blendPsi
  
          IF (iteration == 1) THEN
             errorfirstpsi = sumdb
             diffmxfirstpsi = diffmx
          END IF

          WRITE(*,*) ' itout ',' blendAlpha ',' blendPsi ',' itAlpha ',' diffAlpha ',' errorAlpha ',&
                     ' itPsi ',' diffPsi ',' errorPsi '
          WRITE(*,*) iteration, blendAlpha, blendPsi, nisave1,sumdb1,errorAlpha/twopi_d,nisave,sumdb,errorPsi/MAXVAL(ABS(psival))  !C relative diffmx errors now
  
          IF (iAMR == 1) THEN
             CALL findR
             CALL InterpolatePsiR
             CALL mappsi(0)  ! Full mapping needed, changed psis
             CALL psiFunctions
             CALL maptheta
          ENDIF

          ! Need to set an actual convergence criteria using the JxB and GradP
          ! values calculated in the Write_Convergence_Anisotropic subroutine.
          ! For now we will just use the residuals of the SOR calculations
          If ((errorAlpha.lt.decreaseConvAlpha).and.(errorPsi.lt.decreaseConvPsi)) then
             iConvGlobal = 1
          endif

          IF (((iteration.lt.numit).AND.(iConvGlobal.eq.0)).OR.(iteration.lt.MinSCBIterations)) THEN
             iteration = iteration + 1
             CYCLE Outeriters
          END IF

          IF (iConvGlobal == 1) THEN
             PRINT*, 'Approaching convergence.'
             sumdbconv = sumdb1
             EXIT Outeriters
          END IF
 
       END IF SCB_CALCULATION
  
       EXIT Outeriters
    END DO Outeriters
  
    if (method /= 3) then
       iConv = 1
       iConvGlobal = 1
  
       !   The end of the iterative calculation
       PRINT*, iteration, " outer iterations performed."
       PRINT*, ' '
    end if
  
    IF (boundary /= 'SWMF') PRINT*, "End of calculation."

    IF (iteration > numit) lconv = 1
    nitry = nisave
  
    ! The following block should be uncommented for applications where equal-arc-length is needed 
    ! or desirable
    !C  constTheta = 0.0_dp ! For equal-arc length
    !C  chiVal = (thetaVal + constTheta * SIN(2.*thetaVal)) 
    !C  CALL maptheta
    !C  CALL metrica(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9) 
    ! This will compute the new Bfield on the new grid, to get the right pressure mapping
  
    CALL pressure 
    CALL entropy(entropyFixed, fluxVolume, iCountEntropy)

    ! Compute physical quantities: currents, field components etc..
    CALL metrics
    ! extrapolate to the fixed boundary
    CALL bounextp

    IF (isFBDetailNeeded == 1) CALL Write_Convergence_Anisotropic('01')
    ! Computes energies and Dst from DPS relation, write to disk (+ Biot-Savart values) 
    ! Remove for speed
    IF (isotropy == 0 .AND. isEnergDetailNeeded == 1) CALL dps_general

    RETURN
  
  END SUBROUTINE scb_run
  
!==============================================================================
  SUBROUTINE findR
    ! Obtains equidistant x-s on the equatorial plane, at a local time where the
    ! plasma gradients are strong (e.g. between midnight and dusk
    ! during a storm main phase)
    !!!! Module Variables
    USE ModScbParams,    ONLY: iAzimOffset
    USE ModScbGrids,     ONLY: npsi, nzeta
    use ModScbVariables, ONLY: x, y, nThetaEquator, nZetaMidnight
    !!!! NR Modules
    use nrtype, ONLY: DP

    IMPLICIT NONE
  
    REAL(DP) :: deltaR, distConsecFluxSqOld, distConsecFluxSq
    REAL(DP) :: radius(npsi)
    INTEGER :: j, k
  
    distConsecFluxSq = 0._dp
    distConsecFluxSqOld = 0._dp
  
    IF (iAzimOffset == 2) THEN
       DO k = 2, nzeta
          do j = 2, npsi
             distConsecFluxSq = (x(nThetaEquator,j,k) - x(nThetaEquator,j-1,k))**2 &
                               +(y(nThetaEquator,j,k) - y(nThetaEquator,j-1,k))**2
             IF (distConsecFluxSq > distConsecFluxSqOld) THEN
                distConsecFluxSqOld = distConsecFluxSq
                kMax = k ! (3*nzeta)/8 < what???? not even an integer with nzeta = 61
             END IF
          END DO
       end DO
       DO j = 1, npsi
          radius(j) = SQRT(x(nThetaEquator,j,kMax)**2 + &
               y(nThetaEquator,j,kMax)**2)
       END DO
    ELSE IF (iAzimOffset == 1) THEN
       DO j = 1, npsi
          radius(j) = SQRT(x(nThetaEquator,j,nZetaMidnight)**2 + &
               y(nThetaEquator,j,nZetaMidnight)**2)
       END DO
    END IF
  
    deltaR = (radius(npsi) - radius(1)) / REAL(npsi-1,dp)
  
    DO j = 2, npsi
       radEqMidNew(j) = radius(1) + REAL(j-1,dp)*deltaR
    END DO
  
    radEqMidNew(1) = radius(1)
  
    RETURN
  
  END SUBROUTINE findR
  
!==============================================================================
  SUBROUTINE energy
    !!!! Module Variables
    USE ModScbParams,    ONLY: isotropy
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, dr, dt, dpPrime
    use ModScbVariables, ONLY: pressure3D, x, jacobian, bsq
    !!!! NR Modules
    use nrtype, ONLY: DP

    IMPLICIT NONE
  
    INTEGER  :: i, j, k, iplx
    REAL(DP) :: magneticEnergy, thermalEnergy, totalEnergy, volumeTotal
  
    magneticEnergy = 0.0_dp
    thermalEnergy = 0.0_dp
    totalEnergy = 0.0_dp
    volumeTotal = 0.0_dp
  
    IF (isotropy == 1) THEN
       DO i = 2, nthe-1
          DO j = 2, npsi-1
             DO k = 2, nzeta
                IF ((ABS(x(i,j,k)) > 4._dp)) THEN
                   ! Only for |X| > 5 R_E, closer distances not considered 
                   magneticEnergy = magneticEnergy + jacobian(i,j,k) * dr * dpPrime * dt * &
                        & 0.5_dp * bsq(i,j,k)
                   thermalEnergy = thermalEnergy + jacobian(i,j,k) * dr * dpPrime * dt * &
                        & 1.5_dp * pressure3D(i,j,k)
                   volumeTotal = volumeTotal + jacobian(i,j,k) * dr * dpPrime * dt
                END IF
             END DO
          END DO
       END DO
    END IF
  
    totalEnergy = magneticEnergy + thermalEnergy
  
    WRITE(*, '(A, 1X, E12.3)') 'Magnetic energy = ', magneticEnergy
    WRITE(*, '(A, 1X, E12.3)') 'Thermal energy = ', thermalEnergy
    WRITE(*, '(A, 1X, E12.3)') 'Total energy = ', totalEnergy
    WRITE(*, '(A, 1X, E12.3)') 'Total volume = ', volumeTotal
  
    RETURN
  
  END SUBROUTINE energy
  
!==============================================================================
  SUBROUTINE entropy(ent_local, vol_local, iteration_local)
    !!!! Module Variables 
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, dt
    use ModScbVariables, ONLY: x, y, z, xx, yy, jacobian, bf, nThetaEquator, f, &
                               fzet, rhoVal, thetaVal, zetaVal, psiVal, pjconst, &
                               r0Start
    !!!! Module Subroutines/Functions
    use ModScbSpline, ONLY: Spline_2D_derivs, Spline_coord_derivs
    !!!! NR Modules
    use nrtype, ONLY: DP
    use ezcdf

    IMPLICIT NONE
  
    real(DP) :: ent_local(:,:), vol_local(:,:)
    integer, intent(IN) :: iteration_local
  
    INTEGER :: i, j, k, ierr, idealerr, ncdfId
    REAL(DP) :: yyp, phi, deltaPhi
    REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: derivXTheta, derivXRho, derivXZeta, &
         & derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
         & gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
         gradThetaY, gradThetaZ, gradThetaSq, derivBsqRho, derivBsqZeta
    ! gradRhoSq, gradRhoGradZeta are global
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dVoldXEq, dVoldYEq, dVoldZeta, dVoldAlpha, dVoldRho, &
         dVoldPsi, dEntdXEq, dEntdYEq, dEntdZeta, dEntdAlpha, dEntdRho, &
         dEntdPsi, facVasGlobal, secondTermB
    INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
    REAL(DP) :: delS
    REAL(DP) :: rr1, rr2, zangle, thangle, thangleOnEarth, rr, dza, dya
    REAL(DP) ::  dipoleFactor, dipoleFactor4RE, factorIncrease
  
    ALLOCATE(dVoldZeta(npsi,nzeta), dVoldAlpha(npsi,nzeta), &
         dVoldRho(npsi,nzeta), dVoldPsi(npsi,nzeta), dEntdRho(npsi,nzeta), dEntdPsi(npsi,nzeta), &
         dEntdZeta(npsi,nzeta), dEntdAlpha(npsi,nzeta))
  
    vol_local = 0.0_dp
    if (iteration_local == 1) ent_local = 0.0_dp
  
    DO j = 1, npsi
       DO k = 1, nzeta
          DO i = 1, nThetaEquator-1  ! only (nThetaEquator-1) integration intervals
             vol_local(j,k) = vol_local(j,k) + jacobian(i,j,k) * 1./(f(j)*fzet(k)) * dt
          END DO
          if (iteration_local == 1) ent_local(j,k) = pressure3D(nThetaEquator,j,k) * (vol_local(j,k))**(5./3.)
       END DO
    END DO
  
    ! Compute grad(fluxVolume) components
  
    CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), vol_local(:,1:nzeta), &
         dVoldRho(:,1:nzeta), dVoldZeta(:,1:nzeta))
    CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), ent_local(:,1:nzeta), &
         dEntdRho(:,1:nzeta), dEntdZeta(:,1:nzeta))
  
    DO j = 1, npsi
       dVoldPsi(j,:) = 1._dp / f(j) * dVoldRho(j,:)
       dEntdPsi(j,:) = 1._dp / f(j) * dEntdRho(j,:)
    END DO
  
    DO k = 1, nzeta
       dVoldAlpha(:,k) = dVoldZeta(:,k) / fzet(k)
       dEntdAlpha(:,k) = dEntdZeta(:,k) / fzet(k)
    END DO
  
    ! Allocate derivXTheta etc.
  
    ALLOCATE(derivXTheta(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivXRho(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivXZeta(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivYTheta(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivYRho(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivYZeta(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivZTheta(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivZRho(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(derivZZeta(nthe, npsi, nzeta), STAT = ierr)
  
    CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), derivXTheta, &
         derivXRho, derivXZeta)
    CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), derivYTheta, &
         derivYRho, derivYZeta)
    CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), derivZTheta, &
         derivZRho, derivZZeta)
    ! Now I have all the point derivatives
  
    ! Time to build the Jacobian
  
    jacobian = derivXRho * (derivYZeta * derivZTheta - derivYTheta * derivZZeta) + derivXZeta * &
         & (derivYTheta * derivZRho - derivYRho * derivZTheta) + derivXTheta * &
         & (derivYRho * derivZZeta - derivYZeta * derivZRho)
  
    ! allocate gradRhoX, etc.
  
    ALLOCATE(gradRhoX(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradRhoY(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradRhoZ(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradZetaX(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradZetaY(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradZetaZ(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradThetaX(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradThetaY(nthe, npsi, nzeta), STAT = ierr)
    ALLOCATE(gradThetaZ(nthe, npsi, nzeta), STAT = ierr)
  
    gradRhoX = (derivYZeta * derivZTheta - derivYTheta * derivZZeta) / jacobian
    gradRhoY = (derivZZeta * derivXTheta - derivZTheta * derivXZeta) / jacobian
    gradRhoZ = (derivXZeta * derivYTheta - derivXTheta * derivYZeta) / jacobian
  
    gradZetaX = (derivYTheta * derivZRho - derivYRho * derivZTheta) / jacobian
    gradZetaY = (derivZTheta * derivXRho - derivZRho * derivXTheta) / jacobian
    gradZetaZ = (derivXTheta * derivYRho - derivXRho * derivYTheta) / jacobian
  
    gradThetaX = (derivYRho * derivZZeta - derivYZeta * derivZRho) / jacobian
    gradThetaY = (derivZRho * derivXZeta - derivZZeta * derivXRho) / jacobian
    gradThetaZ = (derivXRho * derivYZeta - derivXZeta * derivYRho) / jacobian
  
    gradRhoSq = gradRhoX**2 + gradRhoY**2 + gradRhoZ**2
    gradRhoGradZeta = gradRhoX * gradZetaX + gradRhoY * gradZetaY + gradRhoZ * gradZetaZ
    gradRhoGradTheta = gradRhoX * gradThetaX + gradRhoY * gradThetaY + gradRhoZ * gradThetaZ
  
    gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ
  
    gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2
  
    DEALLOCATE(derivXTheta, STAT = idealerr)
    DEALLOCATE(derivXRho, STAT = idealerr)
    DEALLOCATE(derivXZeta, STAT = idealerr)
    DEALLOCATE(derivYTheta, STAT = idealerr)
    DEALLOCATE(derivYRho, STAT = idealerr)
    DEALLOCATE(derivYZeta, STAT = idealerr)
    DEALLOCATE(derivZTheta, STAT = idealerr)
    DEALLOCATE(derivZRho, STAT = idealerr)
    DEALLOCATE(derivZZeta, STAT = idealerr)
  
    ALLOCATE(dVoldXEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dVoldYEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dEntdXEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dEntdYEq(npsi,nzeta), stat = ierr)
  
    DO j = 1, npsi
       DO k = 1, nzeta
          dVoldXEq(j,k) = dVoldAlpha(j,k) * fzet(k) * gradZetaX(nThetaEquator,j,k) + &
               dVoldPsi(j,k) * f(j) * gradRhoX(nThetaEquator,j,k)
          dVoldYEq(j,k) = dVoldAlpha(j,k) * fzet(k) * gradZetaY(nThetaEquator,j,k) + &
               dVoldPsi(j,k) * f(j) * gradRhoY(nThetaEquator,j,k)
          dEntdXEq(j,k) = dEntdAlpha(j,k) * fzet(k) * gradZetaX(nThetaEquator,j,k) + &
               dEntdPsi(j,k) * f(j) * gradRhoX(nThetaEquator,j,k)
          dEntdYEq(j,k) = dEntdAlpha(j,k) * fzet(k) * gradZetaY(nThetaEquator,j,k) + &
               dEntdPsi(j,k) * f(j) * gradRhoY(nThetaEquator,j,k)
       END DO
    END DO
  
  
    DEALLOCATE(gradRhoX, STAT = ierr)
    DEALLOCATE(gradRhoY, STAT = ierr)
    DEALLOCATE(gradRhoZ, STAT = ierr)
    DEALLOCATE(gradZetaX, STAT = ierr)
    DEALLOCATE(gradZetaY, STAT = ierr)
    DEALLOCATE(gradZetaZ, STAT = ierr)
    DEALLOCATE(gradThetaX, STAT = ierr)
    DEALLOCATE(gradThetaY, STAT = ierr)
    DEALLOCATE(gradThetaZ, STAT = ierr)
  
    ALLOCATE(facVasGlobal(npsi,nzeta), stat = ierr)
    ALLOCATE(secondTermB(npsi,nzeta), stat = ierr)
  
    DO j = 1, npsi
       DO k = 1, nzeta
          facVasGlobal(j,k) = (dVoldPsi(j,k)*dPdAlpha(nThetaEquator,j,k) - &
               dVoldAlpha(j,k)*dPdPsi(nThetaEquator,j,k)) * bf(1,j,k)
          secondTermB(j,k) = jacobian(1,j,k)/bf(1,j,k)* &
               (f(j)*fzet(k)**2 * gradRhoGradTheta(1,j,k)*gradZetaSq(1,j,k)*dPdAlpha(1,j,k) + &
               f(j)**2 * fzet(k) * gradRhoGradTheta(1,j,k) * gradRhoGradZeta(1,j,k)*dPdPsi(1,j,k) - &
               f(j)*fzet(k)**2 * gradRhoGradZeta(1,j,k) * gradThetaGradZeta(1,j,k)*dPdAlpha(1,j,k) - &
               f(j)**2 * fzet(k) * gradRhoSq(1,j,k) * gradThetaGradZeta(1,j,k) * dPdPsi(1,j,k))
       END DO
    END DO
  
    DO k = 1, nzeta
       DO j = 1, npsi
          xx(1,j,k) = SQRT(x(1,j,k)**2 + y(1,j,k)**2)
          rr2 = x(1,j,k)**2 + y(1,j,k)**2 + z(1,j,k)**2
          rr1 = SQRT(rr2)
          ! thangle is the polar angle on the inner sphere delimiting the
          ! computational domain
          thangle = ASIN(z(1,j,k) / rr1)
          ! thangleOnEarth is the polar angle on Earth's surface
          thangleOnEarth = ACOS(SQRT((COS(thangle))**2/r0Start))
  
          ! Dipole field at the Earth's surface
          dipoleFactor = SQRT(1. + 3. * (SIN(thangleOnEarth))**2)
          dipoleFactor4RE = SQRT(1. + 3. * (SIN(thangle))**2)
          factorIncrease = dipoleFactor * r0Start**3 / dipoleFactor4RE
  
          IF (INT(r0Start) /= 1) facVasGlobal(j,k) = facVasGlobal(j,k) * factorIncrease
       END DO
    END DO

!    CALL cdf_open(ncdfId, 'entropy.cdf', 'w')

!    dimlens(1) = npsi
!    dimlens(2) = 1
!    CALL cdf_define(ncdfId, 'psival', dimlens, 'R8')

!    dimlens(1) = npsi
!    dimlens(2) = nzeta
!    CALL cdf_define(ncdfId, 'xEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'yEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'fluxVolume', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'dVoldXEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'dVoldYEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'entropy', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'dEntdXEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'dEntdYEq', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'facVasG', dimlens, 'R8')
!    CALL cdf_define(ncdfId, 'secondTermB', dimlens, 'R8')

!    CALL cdf_write(ncdfId, 'psival', psival(1:npsi))
!    CALL cdf_write(ncdfId, 'xEq', x(nThetaEquator,1:npsi,1:nzeta))
!    CALL cdf_write(ncdfId, 'yEq', y(nThetaEquator,1:npsi,1:nzeta))
!    CALL cdf_write(ncdfId, 'fluxVolume', vol_local(1:npsi, 1:nzeta))
!    CALL cdf_write(ncdfId, 'dVoldXEq', dVoldXEq(1:npsi, 1:nzeta))
!    CALL cdf_write(ncdfId, 'dVoldYEq', dVoldYEq(1:npsi, 1:nzeta))
!    CALL cdf_write(ncdfId, 'entropy', ent_local(1:npsi,1:nzeta))
!    CALL cdf_write(ncdfId, 'dEntdXEq', dEntdXEq(1:npsi, 1:nzeta))
!    CALL cdf_write(ncdfId, 'dEntdYEq', dEntdYEq(1:npsi, 1:nzeta))
!    CALL cdf_write(ncdfId, 'facVasG', facVasGlobal(1:npsi, 1:nzeta)*pjconst)
!    CALL cdf_write(ncdfId, 'secondTermB', secondTermB(1:npsi, 1:nzeta)*pjconst)

!    CALL cdf_close(ncdfId)
 
    DEALLOCATE(facVasGlobal, stat = idealerr)
    DEALLOCATE(secondTermB, stat = idealerr)
  
    DEALLOCATE(dVoldXEq, stat = idealerr)
    DEALLOCATE(dVoldYEq, stat = idealerr)
    DEALLOCATE(dEntdXEq, stat = idealerr)
    DEALLOCATE(dEntdYEq, stat = idealerr)
  
    ! Can de-allocate derivXRho etc.
  
  RETURN
  
  END SUBROUTINE entropy
  
!==============================================================================
  SUBROUTINE bounextp
    !!!! Module Variables
    USE ModScbGrids,     ONLY: nthe, npsi, npsim, nzetap
    use ModScbVariables, ONLY: bf, bsq, phij
    !!!! Module Subroutines/Functions
    use ModScbFunctions, ONLY: extap
    !!!! NR Modules
    use nrtype, ONLY: DP
    implicit none
  
    integer :: i,j,k
  
    DO k = 1,nzetap
       DO j = 2,npsim
          CALL extap(bj(4,j,k),bj(3,j,k),bj(2,j,k),bj(1,j,k))
          CALL extap(bj(nthe-3,j,k),bj(nthe-2,j,k),bj(nthe-1,j,k),bj(nthe,j,k))
  
          CALL extap(phij(4,j,k),phij(3,j,k),phij(2,j,k),phij(1,j,k))
          CALL extap(phij(nthe-3,j,k),phij(nthe-2,j,k),phij(nthe-1,j,k),phij(nthe,j,k))
  
!          CALL extap(bf(4,j,k),bf(3,j,k),bf(2,j,k),bf(1,j,k))
!          CALL extap(bf(nthe-3,j,k),bf(nthe-2,j,k),bf(nthe-1,j,k) ,bf(nthe,j,k))
!          bsq(1,j,k)=bf(1,j,k)**2
!          bsq(nthe,j,k)=bf(nthe,j,k)**2
       END DO
       DO i=1,nthe
          CALL extap(bj(i,4,k),bj(i,3,k),bj(i,2,k),bj(i,1,k))
          CALL extap(bj(i,npsi-3,k),bj(i,npsi-2,k),bj(i,npsi-1,k),bj(i,npsi,k))
          
          CALL extap(phij(i,4,k),phij(i,3,k),phij(i,2,k),phij(i,1,k))
          CALL extap(phij(i,npsi-3,k),phij(i,npsi-2,k),phij(i,npsi-1,k),phij(i,npsi,k))
  
!          CALL extap(bf(i,4,k),bf(i,3,k),bf(i,2,k),bf(i,1,k))
!          CALL extap(bf(i,npsi-3,k),bf(i,npsi-2,k),bf(i,npsi-1,k) ,bf(i,npsi,k))
!          IF (bf(i,npsi,k) < 0._dp) bf(i,npsi,k) = bf(i,npsi-1,k)
!          bsq(i,1,k) = bf(i,1,k)**2
!          bsq(i,npsi,k) = bf(i,npsi,k)**2
       END DO
    END DO
  
    RETURN
  
  END SUBROUTINE bounextp
  
!==============================================================================
  SUBROUTINE dps_general
    ! Calculates the depression on the Earth's surface - generalized
    ! Dessler-Parker-Sckopke relationship (Siscoe, 1970)
    ! Uses anisotropic pressure
    !!!! Module Variables
    USE ModScbMain,      ONLY: mu0, REarth, BEarth
    USE ModScbGrids,     ONLY: nzeta, npsi, nthe, dr, dt, dpPrime
    use ModScbVariables, ONLY: x, y, z, bsq, jacobian, pnormal, bnormal
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d

    IMPLICIT NONE
  
    INTEGER :: i, j, k, iplx
    REAL(DP) :: magneticEnergy(nthe), magneticEnergyInsideGeo(nthe), &
                magneticEnergyDipole, thermalEnergy(nthe), thermalEnergyInsideGeo(nthe), &
                rsq, totalEnergy, volumeTotal
  
    magneticEnergy = 0.0_dp
    magneticEnergyInsideGeo = 0.0_dp
    magneticEnergyDipole = 0.0_dp
    thermalEnergy = 0.0_dp
    thermalEnergyInsideGeo = 0.0_dp
    totalEnergy = 0.0_dp
    volumeTotal = 0.0_dp
  
    DO i = 2, nthe-1
       DO j = 2, npsi
          DO k = 2, nzeta
             rsq = x(i,j,k)**2 + y(i,j,k)**2
             magneticEnergy(i) = magneticEnergy(i) + jacobian(i,j,k) * dr * dpPrime * dt &
                               * (bsq(i,j,k)*bnormal**2) / (2._dp*mu0) * 1.E-18_dp*REarth**3 ! In Joules now
             IF (rsq < 6.6_dp**2) then
                magneticEnergyInsideGeo(i) = magneticEnergyInsideGeo(i) + jacobian(i,j,k) &
                                           * dr * dpPrime * dt * (bsq(i,j,k)*bnormal**2)  &
                                           / (2._dp*mu0) * 1.E-18_dp*REarth**3
             ENDIF
             thermalEnergy(i) = thermalEnergy(i) + jacobian(i,j,k) * dr * dpPrime * dt &
                              * (pper(i,j,k) + 0.5_dp*ppar(i,j,k))*pnormal * 1.E-9_dp  &
                              * REarth**3 ! In Joules now
             IF (rsq < 6.6_dp**2) then
                thermalEnergyInsideGeo(i) = thermalEnergyInsideGeo(i) + jacobian(i,j,k) &
                                          * dr * dpPrime * dt * (pper(i,j,k) &
                                          + 0.5_dp*ppar(i,j,k))*pnormal * 1.E-9_dp * REarth**3
             ENDIF
             volumeTotal = volumeTotal + jacobian(i,j,k) * dr * dpPrime * dt
          END DO
       END DO
    END DO
  
    magneticEnergyDipole = 4._dp*pi_d/(3._dp*mu0) * BEarth**2 * REarth**3 ! In J
  
    totalEnergy = SUM(magneticEnergy) + SUM(thermalEnergy)
    DstDPS = 1.3_dp * (-BEarth) * (2._dp*SUM(thermalEnergy))/(3._dp*magneticEnergyDipole) * 1.E9_dp
    DstDPSInsideGeo = 1.3_dp * (-BEarth) * (2._dp*SUM(thermalEnergyInsideGeo))/(3._dp*magneticEnergyDipole) * 1.E9_dp
    WRITE(*, '(A, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, A)') 'DstDPS, DstDPSGeo, DstBiot, DstBiotGeo = ', real(DstDPS), &
         real(DstDPSInsideGeo), real(DstBiot), real(DstBiotInsideGeo), ' nT' ! 1.3 factor due to currents induced in the Earth 
  
    RETURN
  
  END SUBROUTINE dps_general
  
!==============================================================================
  SUBROUTINE mapTheta
    !!!! Module Variables
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: diffmx, rjac, nisave,  x, y, z, sumb, sumdb, chiVal, &
                               chi
    !!!! Module Subroutines/Functions
    use ModScbSpline, ONLY: spline, splint
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d
  
    IMPLICIT NONE
  
    REAL(DP), DIMENSION(nthe) :: xOld, yOld, zOld, distance, chiValOld, chi2derivsX, &
                                 chi2derivsY, chi2derivsZ
    REAL(DP), DIMENSION(nthe,npsi,nzetap) :: chiPrev

    INTEGER :: i, j, k
  
    !     now move theta coordinates along each surface
    !     equal arc length along the i grids
  
    chiPrev(:,:,:) = chi(:,:,:)

    zetaloop: DO k = 2, nzeta
       fluxloop: DO j = 1, npsi
          distance(1) = 0._dp
          xOld(:) = x(1:nthe,j,k)
          yOld(:) = y(1:nthe,j,k)
          zOld(:) = z(1:nthe,j,k)
          chiValOld(1) = 0._dp
  
          DO i = 2, nthe
             distance(i) = distance(i-1) + SQRT((x(i,j,k)-x(i-1,j,k))**2 &
                  & +(y(i,j,k)-y(i-1,j,k))**2 +(z(i,j,k)-z(i-1,j,k))**2)
          END DO
  
          chiValOld = distance / distance(nthe) * pi_d
  
          ! "Natural" splines
          CALL spline(chiValOld, xOld, 1.E31_dp, 1.E31_dp, chi2derivsX)
          CALL spline(chiValOld, yOld, 1.E31_dp, 1.E31_dp, chi2derivsY)
          CALL spline(chiValOld, zOld, 1.E31_dp, 1.E31_dp, chi2derivsZ)
  
          DO i = 2, nthem
             x(i,j,k) = splint(chiValOld, xOld, chi2derivsX, chiVal(i))
             y(i,j,k) = splint(chiValOld, yOld, chi2derivsY, chiVal(i))
             z(i,j,k) = splint(chiValOld, zOld, chi2derivsZ, chiVal(i))
!             x(i,j,k) = splint(chiValOld, xOld, chi2derivsX, chiPrev(i,j,k))
!             y(i,j,k) = splint(chiValOld, yOld, chi2derivsY, chiPrev(i,j,k))
!             z(i,j,k) = splint(chiValOld, zOld, chi2derivsZ, chiPrev(i,j,k))
             chi(i,j,k) = chiVal(i)
          END DO
       END DO fluxloop
    END DO zetaloop
  
    !  periodic boundary conditions
    x(:,:,1) = x(:,:,nzeta)
    y(:,:,1) = y(:,:,nzeta)
    z(:,:,1) = z(:,:,nzeta)
    chi(:,:,1) = chi(:,:,nzeta)
    x(:,:,nzetap) = x(:,:,2)
    y(:,:,nzetap) = y(:,:,2)
    z(:,:,nzetap) = z(:,:,2)
    chi(:,:,nzetap) = chi(:,:,2)
  
  !C z(nThetaEquator,:,:) = 0.0_dp ! Symmetry
  
    RETURN
  
  END SUBROUTINE mapTheta
  
!==============================================================================
!******************************************************************************
SUBROUTINE pressure
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************
    !!!! Module Variables
    USE ModRamVariables, ONLY: PParH, PPerH, PParO, PPerO, PParHe, PPerHe, PParE, &
                               PPerE, PHI, LZ
    use ModRamParams,    ONLY: boundary
    use ModScbMain,      ONLY: iCountPressureCall, iSm2
    use ModScbParams,    ONLY: iLossCone, iOuterMethod, iReduceAnisotropy, isotropy
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap, nXRaw, nXRawExt, nYRaw, &
                               nAzimRAM
    use ModScbVariables, ONLY: x, y, z, xx, yy, bf, bsq, rhoVal, zetaVal, thetaVal, &
                               nZetaMidnight, nThetaEquator, pnormal, f, fzet, alfa, &
                               dela, azimRaw, radGrid, angleGrid, ratioEq, dPPerdRho, &
                               dPPerdZeta, dPPerdTheta, dBsqdRho, dBsqdZeta
    !!!! Module Subroutines/Functions
    USE ModScbSpline, ONLY: Spline_2D_derivs, Spline_2D_point, Spline_coord_Derivs
    USE ModScbFunctions, ONLY: SavGol7, pRoeRad, extap
    !!!! NR Modules
    use nrtype, ONLY: DP, SP, pi_d, twopi_d
    IMPLICIT NONE
  
    INTEGER :: i, iloopOut, ierflg, j, j1, k1, jSKBoundary, k, ierr, ierrDom, idealerr, m1, mstate, n1
  
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
    !C INTEGER, PARAMETER :: nXRawExt = 171, nAzimRAM = 49 ! For DMSP runs
    REAL(DP) :: xRaw(nXRaw,nYRaw), YRaw(nXRaw,nYRaw),pressProtonPerRaw(nXRaw,nYRaw), pressProtonParRaw(nXRaw,nYRaw), &
                pressOxygenPerRaw(nXRaw,nYRaw), pressOxygenParRaw(nXRaw,nYRaw), pressHeliumPerRaw(nXRaw,nYRaw), &
                pressHeliumParRaw(nXRaw,nYRaw), pressPerRaw(nXRaw,nYRaw), pressParRaw(nXRaw,nYRaw), &
                pressEleParRaw(nXRaw,nYRaw), pressElePerRaw(nXRaw,nYRaw), &     !Vania
                radRaw_local(nXRaw), ratioRaw(nXRaw,nYRaw), &
                radRoe(nXRoe), azimRoe(nYRoe), energRoe(0:nEnergRoe), PARoe(nPARoe), fluxRoe(nXRoe, nYRoe, nEnergRoe, 18), &
                pressProtonPerRoe(nXRoe, nYRoe), pressProtonParRoe(nXRoe, nYRoe), &
                pressPerRoe(nXRoe, nYRoe), pressParRoe(nXRoe, nYRoe), ratioRoe(nXRoe, nYRoe)
    REAL(DP) :: xRawExt(nXRawExt,nAzimRAM), YRawExt(nXRawExt,nAzimRAM),pressPerRawExt(nXRawExt,nAzimRAM), & 
         pressParRawExt(nXRawExt,nAzimRAM), radRawExt(nXRawExt), azimRawExt(nAzimRAM), ratioRawExt(nXRawExt,nAzimRAM)
    REAL(DP), PARAMETER :: l0 = 50._dp
    CHARACTER(len=93)  :: firstLine, secondLine
    CHARACTER(len=200) :: header
    !  INTEGER, PARAMETER :: ISLIM = nXRaw*nYRaw, NUMXOUT = npsi, NUMYOUT = nzeta-1
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
    REAL :: pressPerRawRowExt(nXRawExt*nAzimRAM), pressParRawRowExt(nXRawExt*nAzimRAM)
  
    REAL :: wrk(lwrk), wrk1(lwrk1), wrk2(lwrk2)
    REAL :: fpResids, fpResidsPer, fpResidsPar
    REAL :: smoothFactor, smoothFactorPer, smoothFactorPar  
    INTEGER, PARAMETER :: kwrk = 50000, kwrk1 = 5000!, nuest = nXRaw+7, nvest = nYRaw+7 !C nuest = nXRawExt + 7, nvest = nAzimRAM+7
    INTEGER, PARAMETER :: kx = 3, ky = 3  ! Must be 3 for polar, can vary for surfit
    !C  INTEGER, PARAMETER :: nxest = 24, nyest = 12, nmax = MAX(nxest, nyest)
    INTEGER, PARAMETER :: nxest = 15, nyest = 15, nmax = MAX(nxest, nyest) 
    ! For Roeder expansion, not wise to go for larger nx, ny as it might force an unnatural spline
    REAL :: coeff((nXRaw+7-4)*(nYRaw+7-4))
    REAL, SAVE :: coeff1((nxest-kx-1)*(nyest-ky-1)), coeff2((nxest-kx-1)*(nyest-ky-1))
    INTEGER :: iwrk(kwrk), iwrk1(kwrk1)
    REAL :: tu(nXRaw+7), tv(nYRaw+7)
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
    REAL,  ALLOCATABLE :: r(:), xSp(:), ySp(:), pValue(:), pValuePer(:), pValuePar(:), weight(:), &
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
    REAL(DP) :: xAr(48), yAr(48)
    REAL(DP) :: xSWMF(48,48), ySWMF(48,48), pressSWMF(48,48), rhoSWMF(48,48) 
    CHARACTER(len=4) :: ST3
    real(dp) :: drad

    integer :: mloc(2)
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
          dela(nzetap) = dela(2)
          xx(i,j,1) = xx(i,j,nzeta)
          xx(i,j,nzetap) = xx(i,j,2)
          yy(i,j,1) = alfa(i,j,1) - dela(1)
          yy(i,j,nzetap) = alfa(i,j,nzetap) - dela(nzetap)
       END DO
    END DO
  
  
    Isotropy_choice:  IF (isotropy == 1) THEN    ! isotropic case
! For now this takes the anisotropic pressure from RAM and used it for the SCB
! calculation. For some reason the isotropic case of RAM-SCB never worked, or at
! least never worked with the version of the code I started with. -ME
!          DO j1 = 1, nXRaw
!             DO k1 = 1, nYRaw
!                radRaw_local(j1) = LZ(j1+1)
!                azimRaw(k1) = PHI(k1)*12/pi_d
!                pressProtonPerRaw(j1,k1) = PPERH(j1+1,k1)
!                pressProtonParRaw(j1,k1) = PPARH(j1+1,k1)
!                pressOxygenPerRaw(j1,k1) = PPERO(j1+1,k1)
!                pressOxygenParRaw(j1,k1) = PPARO(j1+1,k1)
!                pressHeliumPerRaw(j1,k1) = PPERHE(j1+1,k1)
!                pressHeliumParRaw(j1,k1) = PPARHE(j1+1,k1)
!                pressElePerRaw(j1,k1)    = PPERE(j1+1,k1)
!                pressEleParRaw(j1,k1)    = PPARE(j1+1,k1)
!             END DO
!          END DO
!
!          azimRaw = azimRaw * 360./24 * pi_d / 180._dp ! In radians
!
!          pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw + pressElePerRaw) ! from keV/cm^3 to nPa
!          pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw + pressEleParRaw) ! from keV/cm^3 to nPa
!
!          ratioRaw = pressOxygenPerRaw/pressProtonPerRaw
!
!
!          radRawExt(1:nXRaw) = radRaw_local(1:nXRaw)
!          DO j1 = nXRaw+1, nXRawExt
!             radRawExt(j1) = radRaw_local(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw_local(nXRaw)-radRaw_local(1))/(REAL(nXRaw-1, DP))
!          END DO
!
!          azimRawExt(1:nAzimRAM) = azimRaw(1:nYRaw) ! nYRaw = nAzimRAM
!
!          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
!          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
!          ratioRawExt(1:nXRaw,:) = ratioRaw(1:nXRaw,:)
!
!          DO k1 = 1, nAzimRAM
!             DO j1 = nXRaw+1, nXRawExt
!                ! Alternatively, extrapolate the pressure assuming SK dependence
!                ! in regions where we don't know it instead of f(R)
!                ! extrapolation; 
!                ! or, decrease the order of the polynomial extrapolation
!                pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
!                     (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
!                pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
!                     (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
!                ratioRawExt(j1,k1) = pressOxygenPerRaw(nXRaw,k1)/pressProtonPerRaw(nXRaw,k1)
!             END DO
!          END DO
!
!          !  pressure functions
!          DO k = 2, nzeta
!             DO j = 1, npsi
!
!                radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
!                angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
!
!                IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
!                     angle = twopi_d - ASIN(y(nThetaEquator,j,k) / radius)
!
!                IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
!                     angle = - ASIN(y(nThetaEquator,j,k) / radius)
!
!                radGrid(j,k) = radius
!                angleGrid(j,k) = angle
!
!             END DO
!          END DO
!
!          IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing (possibly multiple) for the pressure
!             pressPerRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressPerRawExt(1:nXRawExt,1:nAzimRAM))
!             pressParRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressParRawExt(1:nXRawExt,1:nAzimRAM))
!          END IF
!
!          !Cubic spline interpolation
!          CALL Spline_2D_point(radRawExt**2, azimRawExt, pressPerRawExt, &
!               radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain)
!          CALL Spline_2D_point(radRawExt**2, azimRawExt, pressParRawExt, &
!               radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain)
!          CALL Spline_2D_point(radRawExt**2, azimRawExt, ratioRawExt, &
!               radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), ratioEq(1:npsi,2:nzeta), iDomain)
!          ratioEq(1:npsi,1) = ratioEq(1:npsi,nzeta)
!          IF (iDomain > 0) THEN
!             PRINT*, 'Stop; problem with pressure domain; iDomain = ', iDomain
!             STOP
!          END IF
!
!          ! Sometimes the interpolation can give very small negative values very 
!          ! near the Earth; inconsequential
!          WHERE(pperEq < 0.0) pperEq = 1e-2_dp
!          WHERE(pparEq < 0.0) pparEq = 1e-2_dp
!          DO k = 2, nzeta
!             DO j = 1, npsi
!                IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2 RE from Earth
!                   pperEq(j,k) = pressPerRaw(1,1)
!                   pparEq(j,k) = pressParRaw(1,1)
!                END IF
!             END DO
!          END DO
!
!         pperEq(:,nzeta+1) = pperEq(:,2)
!         pparEq(:,nzeta+1) = pparEq(:,2)
!         pperEq(:,1) = pperEq(:,nzeta)
!         pparEq(:,1) = pparEq(:,nzeta)
!
!         pperEq = pperEq/pnormal
!         pparEq = pparEq/pnormal
!
!       press = 1.0/3.0 * pparEq + 2.0/3.0 * pperEq
       do j=1,npsi
          do k=1,nzeta
             press(j,k) = pressureTsygMuk(x(nThetaEquator,j,k), y(nThetaEquator,j,k))
          enddo
       enddo

       CALL Spline_2D_derivs(rhoVal, zetaVal(2:nzeta), press(:,2:nzeta), &
            dPresdRho(:,2:nzeta), dPresdZeta(:,2:nzeta))
       press(:,nzetap) = press(:,2)
       press(:,1) = press(:,nzeta)
       dPresdRho(:,nzetap) = dPresdRho(:,2)
       dPresdRho(:,1) = dPresdRho(:,nzeta)
       dPresdZeta(:,nzetap) = dPresdZeta(:,2)
       dPresdZeta(:,1) = dPresdZeta(:,nzeta)
  
       CALL Spline_2D_derivs(rhoVal, zetaVal(2:nzeta), dPresdRho(:,2:nzeta), &
            dSqPresdRhoSq(:,2:nzeta), dSqPresdRhodZeta(:,2:nzeta))
       CALL Spline_2D_derivs(rhoVal, zetaVal(2:nzeta), dPresdZeta(:,2:nzeta), &
            dSqPresdRhodZeta(:,2:nzeta), dSqPresdZetaSq(:,2:nzeta))
       dSqPresdRhoSq(:,nzetap) = dSqPresdRhoSq(:,2)
       dSqPresdRhoSq(:,1) = dSqPresdRhoSq(:,nzeta)
       dSqPresdZetaSq(:,nzetap) = dSqPresdZetaSq(:,2)
       dSqPresdZetaSq(:,1) = dSqPresdZetaSq(:,nzeta)
  
       DO i = 1, nthe
          DO j = 1, npsi
             dpdPsi(i,j,1:nzetap) = 1._dp / f(j) * dPresdRho(j,1:nzetap)
             IF (iOuterMethod == 2) dSqPdPsiSq(i,j,1:nzetap) = 1._dp / f(j)**2 * dSqPresdRhoSq(j,1:nzetap)
          END DO
          DO k = 1, nzetap
             dpdAlpha(i,1:npsi, k) = dPresdZeta(1:npsi,k) / fzet(k)
             IF (iOuterMethod == 2) dSqPdAlphaSq(i,1:npsi,k) = dSqPresdZetaSq(1:npsi,k) / fzet(k)**2
          END DO
          pressure3D(i,1:npsi,1:nzetap) = press(1:npsi, 1:nzetap)
       END DO
  
  
    ELSE    ! Anisotropic pressure case
       IF ((boundary.eq.'LANL+').or.(boundary.eq.'SWMF')) THEN ! Calculation using RAM pressures
          ! print*, 'rank, 1st line: ', rank, firstLine; call flush(6)
          ! print*, 'rank, 2nd line:', rank, secondLine; call flush(6)
          ! print*, 'HERE'; call flush(6)
          DO j1 = 1, nXRaw
             DO k1 = 1, nYRaw
                radRaw_local(j1) = LZ(j1+1)
                azimRaw(k1) = PHI(k1)*12/pi_d
                pressProtonPerRaw(j1,k1) = PPERH(j1+1,k1)
                pressProtonParRaw(j1,k1) = PPARH(j1+1,k1)
                pressOxygenPerRaw(j1,k1) = PPERO(j1+1,k1)
                pressOxygenParRaw(j1,k1) = PPARO(j1+1,k1)
                pressHeliumPerRaw(j1,k1) = PPERHE(j1+1,k1)
                pressHeliumParRaw(j1,k1) = PPARHE(j1+1,k1)
                pressElePerRaw(j1,k1)    = PPERE(j1+1,k1)
                pressEleParRaw(j1,k1)    = PPARE(j1+1,k1)
             END DO
          END DO
  
          azimRaw = azimRaw * 360./24 * pi_d / 180._dp ! In radians
  
          pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw + pressElePerRaw) ! from keV/cm^3 to nPa
          pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw + pressEleParRaw) ! from keV/cm^3 to nPa
  
          ratioRaw = pressOxygenPerRaw/pressProtonPerRaw
  
  
          radRawExt(1:nXRaw) = radRaw_local(1:nXRaw)
          DO j1 = nXRaw+1, nXRawExt
             radRawExt(j1) = radRaw_local(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw_local(nXRaw)-radRaw_local(1))/(REAL(nXRaw-1, DP))
          END DO
  
          azimRawExt(1:nAzimRAM) = azimRaw(1:nYRaw) ! nYRaw = nAzimRAM
  
          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
          ratioRawExt(1:nXRaw,:) = ratioRaw(1:nXRaw,:)
          
          DO k1 = 1, nAzimRAM
             DO j1 = nXRaw+1, nXRawExt
                ! Alternatively, extrapolate the pressure assuming SK dependence in regions where we don't know it instead of f(R) extrapolation; 
                ! or, decrease the order of the polynomial extrapolation
                pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
                     (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
                pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) / &
                     (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53)) 
                ratioRawExt(j1,k1) = pressOxygenPerRaw(nXRaw,k1)/pressProtonPerRaw(nXRaw,k1)
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
             pressPerRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressPerRawExt(1:nXRawExt,1:nAzimRAM))
             pressParRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressParRawExt(1:nXRawExt,1:nAzimRAM))
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
                END DO
                press(j,k) = pEq
             END DO
          END DO
  
       ELSEIF (boundary == 'LANL') THEN ! RAM pressures & Roeder model extension; only protons for now in the 
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
  
          DO j1 = 1, nXRaw
             DO k1 = 1, nYRaw
                radRaw_local(j1) = LZ(j1)
                azimRaw(k1) = PHI(k1)*12/pi_d
                pressProtonPerRaw(j1,k1) = PPERH(j1,k1)
                pressProtonParRaw(j1,k1) = PPARH(j1,k1)
                pressOxygenPerRaw(j1,k1) = PPERO(j1,k1)
                pressOxygenParRaw(j1,k1) = PPARO(j1,k1)
                pressHeliumPerRaw(j1,k1) = PPERHE(j1,k1)
                pressHeliumParRaw(j1,k1) = PPARHE(j1,k1)
                pressElePerRaw(j1,k1)    = PPERE(j1,k1)
                pressEleParRaw(j1,k1)    = PPARE(j1,k1)
             END DO
          END DO
          azimRaw = azimRaw * 360./24 * pi_d / 180._dp
          pressPerRaw = 0.16_dp * (pressProtonPerRaw + 1.*pressOxygenPerRaw + 1.*pressHeliumPerRaw + 1.*pressElePerRaw)  ! from keV/cm^3 to nPa
          pressParRaw = 0.16_dp * (pressProtonParRaw + 1.*pressOxygenParRaw + 1.*pressHeliumParRaw + 1.*pressEleParRaw)  ! from keV/cm^3 to nPa
  
          radRawExt(1:nXRaw) = radRaw_local(1:nXRaw)
          DO j1 = nXRaw+1, nXRawExt
             radRawExt(j1) = radRaw_local(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw_local(nXRaw)-radRaw_local(1))/(REAL(nXRaw-1, DP))
          END DO
  
          azimRawExt(1:nAzimRAM) = azimRaw(1:nYRaw) ! nYRaw = nAzimRAM
          
          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
         
          DO k1 = 1, nAzimRAM
             DO j1 = nXRaw+1, nXRawExt
                ! Alternatively, extrapolate the pressure assuming SK dependence in regions where we don't know it instead of f(R) extrapolation; 
                ! or, decrease the order of the polynomial extrapolation
                pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
                pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) *  pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
             END DO
          END DO

          IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing (possibly multiple) for the pressure
             pressPerRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressPerRawExt(1:nXRawExt,1:nAzimRAM))
             pressParRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressParRawExt(1:nXRawExt,1:nAzimRAM))
          END IF

          ! Piecewise cubic spline interpolation; alternative - put all points scattered and do surfit
          CALL Spline_2D_point(radRawExt, azimRawExt, pressPerRawExt, &
               radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pperEq(1:npsi,2:nzeta), iDomain) 
          CALL Spline_2D_point(radRawExt, azimRawExt, pressParRawExt, &
               radGrid(1:npsi,2:nzeta), angleGrid(1:npsi,2:nzeta), pparEq(1:npsi,2:nzeta), iDomain) 
          IF (iDomain > 0) THEN
             PRINT*, 'Stop; problem with pressure domain; iDomain = ', iDomain
             STOP
          END IF
  
          IF (ALLOCATED(pValue)) DEALLOCATE(pValue, STAT=ierr)
  
          pperEq = pperEq / pnormal 
          pparEq = pparEq / pnormal
  
          ! Sometimes the interpolation can give very small negative values very 
          ! near the Earth; inconsequential
          WHERE(pperEq < 0.0) pperEq = MINVAL(pressPerRaw) ! 1e-1_dp/pnormal
          WHERE(pparEq < 0.0) pparEq = MINVAL(pressParRaw) ! 1e-1_dp/pnormal
          DO k = 1, nzeta
             DO j = npsi-3, 1, -1
                IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2RE from Earth
                   CALL extap(pperEq(j+3,k),pperEq(j+2,k),pperEq(j+1,k),pperEq(j,k))
                   CALL extap(pparEq(j+3,k),pparEq(j+2,k),pparEq(j+1,k),pparEq(j,k))
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
                      ppar(i,j,k) = pparN * (aN+1._dp)/(1._dp+aN*ratioB) * SQRT((rBI-1._dp)/(rBI-ratioB)) &
                                          * (1._dp-(1.+aN*ratioB)/(rBI+aN*ratioB))
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
                press(j,k) = pEq
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
  
    RETURN
  
  END SUBROUTINE pressure
 
FUNCTION pressureTsygMuk(xEqGsm, yEqGsm)

  ! Tsyganenko-Mukai statistical pressure model
  ! from Tsyganenko, N. and T. Mukai, Tail plasma sheet models derived from
  ! Geotail particle data
  ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 108, NO. A3, 1136,
  ! doi:10.1029/2002JA009707, 2003

  ! Table 1 
  !Average hT i = 3.795 keV hNi = 0.625 cm 3 hPi = 0.229 nPa
  !Model rms 
  ! deviation
  ! 1.422 0.342 0.0427
  ! Correlation 0.708 0.567 0.955
  ! Equation number (4)        (5)         (6)

  ! PSW=1.94E-6*xn_pd*vel**2 

  USE nrtype
  USE ModScbVariables, ONLY: pnormal
  use ModRamIndices,   ONLY: NameOmniFile
  use ModIOUnit,       ONLY: UNITTMP_
  IMPLICIT NONE

  REAL(DP) :: xEqGsm, yEqGsm, pressureTsygMuk
  REAL(DP), PARAMETER :: A1 = 0.057, A2 = 0.524, A3 = 0.0908, A4 = 0.527, &
       A5 = 0.078, A6 = -4.422, A7 = -1.533, A8 = -1.217, A9 = 2.54, &
       A10 = 0.32, A11 = 0.754, A12 = 1.048, A13 = -0.074, A14 = 1.015
  REAL(DP) :: Bperp, rhoStar, rhoStar10, PSWStar, F, FStar, phi, theta, presAt10RE
  integer :: aa, ab, ac, ad, i
  real(DP) :: ba, bb, bc, bd, be, pdynGlobal, byimfGlobal, bzimfGlobal

  pdynGlobal = 2.10
  byimfGlobal = 0.0
  bzimfGlobal = 1.0

  rhoStar = 0.1_dp * SQRT(xEqGsm**2+yEqGsm**2) ! 1/10RE * rho

  phi = ATAN2(xEqGsm, yEqGsm) ! azimuthal angle, atan(y/x)

  PSWStar = pdynGlobal/3._dp

  theta = ATAN2(byimfGlobal,bzimfGlobal) ! IMF clock angle
  IF (theta < 0.0) theta = theta + twopi_d
  Bperp = SQRT(byimfGlobal**2 + bzimfGlobal**2) ! BSW perpendicular to Sun-Earth axis
  F = ABS(Bperp)*SQRT(SIN(0.5_dp*theta))
  FStar = 0.2_dp*F

!C  print*, 'pdyn, theta, F = ', pdynGlobal, theta, F

  pressureTsygMuk = A1*rhoStar**A6 + A2*PSWStar**A11*rhoStar**A7 + A3*FStar**A12*rhoStar**A8 + &
       (A4*PSWStar**A13*EXP(-A9*rhoStar) + A5*FStar**A14*EXP(-A10*rhoStar))*(SIN(phi))**2
  pressureTsygMuk = pressureTsygMuk/pnormal !In computational units

  ! Continuation to R < 10 RE (rhoStar < 1.)
  rhoStar10 = 1._dp
  presAt10RE = A1*rhoStar10**A6 + A2*PSWStar**A11*rhoStar10**A7 + A3*FStar**A12*rhoStar10**A8 + &
       (A4*PSWStar**A13*EXP(-A9*rhoStar10) + A5*FStar**A14*EXP(-A10*rhoStar10))*(SIN(phi))**2
  presAt10RE = presAt10RE/pnormal
  IF (rhoStar < 1._dp) pressureTsygMuk = presAt10RE/pressureRad(10._dp)  * &
       pressureRad(10._dp*rhoStar)*(1._dp-exp(1.*rhoStar-1.)) + pressureTsygMuk*exp(1.*rhoStar-1)

  RETURN

END FUNCTION pressureTsygMuk

FUNCTION pressureRad(radius)

  USE nrtype
  use ModScbVariables, ONLY: pressurequot, pnormal

  IMPLICIT NONE

  integer :: iCorrectedPressure = 1
  REAL(DP) :: radius, LargeA, LA, dPdRRight, dPdR2Right, Acap, Bcap, Ccap, Dcap, &
       Acap2, Bcap2, Ccap2, pressureRad
  REAL(DP) :: m, n, pUp, pDown, x1, x2, delta1, delta2, pressureSK, pUp2, pDown2

  ! This subroutine outputs pressure as a function of radial distance r in the
  ! equatorial plane, midnight 

  ! Spence - Kivelson empirical formula
  IF (iCorrectedPressure /= 1) THEN
     !     IF (ABS(radius) > 3.5)
     pressureRad = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))

     LargeA = pressurequot / pnormal * 89.*EXP(-0.59*3.5) + 8.9 * 3.5**(-1.53) ! value at x = 3
     !     IF (ABS(radius) <= 3.5) pressureRad = pressurequot/
     !     pnormal*LargeA**(((ABS(radius) - 2.5)))
     ! Took this out for now

     RETURN
  END IF

  ! "corrected" Spence-Kivelson pressure for test runs, with a smaller grad p at large distances
  IF (iCorrectedPressure == 1) THEN
     ! This is the profile used for obtaining thin CS (FLR calculations)
     pressureSK = pressurequot / pnormal * (89.*EXP(-0.59*ABS(radius)) + 8.9*ABS(radius)**(-1.53))
     m = 4.
     n = 2.
     ! pUp = m * pressureSK
     pUp = 10. / pnormal
     pDown = n * pressureSK
     x1 = 7.5
     x2 = 9.25
     delta1 = 3.0_dp
     delta2 = 3.0_dp
     pUp2 = pUp / 2. * (1. + TANH((x1 - ABS(radius)) / delta1))
     pDown2 = pDown / 2. * (1. + TANH((ABS(radius) - x2)/delta2))
     pressureRad = pUp2 + pDown2
  END IF

  RETURN

END FUNCTION pressureRad

 
END MODULE ModScbRun
