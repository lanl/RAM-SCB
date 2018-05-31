!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

  MODULE ModScbRun
  ! Contains subroutines responsible for making the SCB calculations

    implicit none

    integer, save :: iteration

    contains
  
!==============================================================================
  SUBROUTINE scb_run(nIter)
    !!!! Module Variables
    use ModRamParams,    ONLY: boundary, electric, NameBoundMag, verbose
    use ModRamVariables, ONLY: KP
    use ModRamTiming,    ONLY: TimeRamNow
    USE ModScbMain,      ONLY: damp, iSm, nrelax, numit, relax, thresh
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, nAzimRAM, nXRawExt
    USE ModScbVariables, ONLY: alfa, alfaSav1, alfaSav2, psi, psiSav1, psiSav2, psiVal, &
                               alphaVal, blendAlpha, blendPsi, iAlphaMove, iPsiMove, &
                               decreaseConvAlpha, decreaseConvPsi, errorAlpha, errorPsi, &
                               diffmx, errorAlphaPrev, errorPsiPrev, x, y, z, sumb, &
                               sumdb, jacobian, xzero3, psiin, psiout, psitot, &
                               xpsiin, xpsiout, f, fp, fluxVolume, alfaPrev, nThetaEquator, &
                               constZ, fzet, fzetp, chiVal, chi, thetaVal, constTheta, &
                               normJxB, normGradP, SORFail, nFail, hICalc, normDiff, &
                               iConvGlobal, lconv, nisave, nitry
    use ModScbParams,    ONLY: decreaseConvAlphaMin, decreaseConvPsiMin, blendMin, &
                               decreaseConvAlphaMax, decreaseConvPsiMax, blendMax, &
                               blendAlphaInit, blendPsiInit, MinSCBIterations, &
                               iAMR, isEnergDetailNeeded, isFBDetailNeeded, &
                               method, isotropy
    !!!! Module Subroutine/Functions
    USE ModRamGSL,      ONLY: GSL_Interpolation_1D, GSL_Smooth_1D
    USE ModScbCompute,  ONLY: computeBandJacob, compute_convergence, metrics
    USE ModScbEuler,    ONLY: alfges, psiges, mapalpha, mappsi, directAlpha, &
                              iterateAlpha, directPsi, iteratePsi, psiFunctions, &
                              InterpolatePsiR, maptheta
    USE ModScbEquation, ONLY: newk, newj, metric, metrica ! LHS and RHS equations
    USE ModScbIO,       ONLY: Write_Convergence_Anisotropic, Update_Domain, Computational_Domain
    !!!! Share Modules
    use ModTimeConvert, ONLY: n_day_of_year
    USE ModIOUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, twopi_d, pi_d

    implicit none
  
    INTEGER, INTENT(IN) :: nIter
    INTEGER  :: iconv, nisave1, ierr, iCountEntropy, GSLerr
    INTEGER  :: i, j, k, SCBIterNeeded
    REAL(DP) :: outDistance, convDistance, blendInitial
    REAL(DP) :: sumdbconv, errorfirstalpha, diffmxfirstalpha, &
                errorfirstpsi, diffmxfirstpsi, dphi, phi, psis, &
                xpsitot, xpl
    REAL(DP) :: sumb1, sumdb1, diffmx1, normDiffPrev
    REAL(DP) :: entropyFixed(npsi,nzeta)
    REAL(DP), DIMENSION(500) :: psiSpline, xSpline, ySpline, zSpline
    REAL(DP), ALLOCATABLE, SAVE :: xPrev(:,:,:), yPrev(:,:,:), zPrev(:,:,:), &
                                   alphaPrev(:,:,:), psiPrev(:,:,:), xStart(:,:,:), &
                                   yStart(:,:,:), zStart(:,:,:), psiStart(:,:,:), &
                                   alphaStart(:,:,:), fStart(:)
    LOGICAL :: check

    ! Variables for timing
    integer :: time1, clock_rate, clock_max
    real(dp) :: starttime,stoptime
    clock_rate = 1000
    clock_max = 100000

    ALLOCATE(xStart(nthe,npsi,nzeta+1), yStart(nthe,npsi,nzeta+1), zStart(nthe,npsi,nzeta+1))
    xStart = 0.0; yStart = 0.0; zStart = 0.0
    ALLOCATE(psiStart(nthe,npsi,nzeta+1), alphaStart(nthe,npsi,nzeta+1), fStart(npsi))
    psiStart = 0.0; alphaStart = 0.0; fStart = 0.0
    IF (.NOT. ALLOCATED(xPrev)) ALLOCATE(xPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    xPrev = 0.0
    IF (.NOT. ALLOCATED(yPrev)) ALLOCATE(yPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    yPrev = 0.0
    IF (.NOT. ALLOCATED(zPrev)) ALLOCATE(zPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
    zPrev = 0.0
    IF (.NOT. ALLOCATED(alphaPrev)) ALLOCATE(alphaPrev(nthe,npsi,nzeta+1), STAT = ierr)
    alphaPrev = 0.0
    IF (.NOT. ALLOCATED(psiPrev)) ALLOCATE(psiPrev(nthe,npsi,nzeta+1), STAT = ierr)
    psiPrev = 0.0

    decreaseConvAlpha = decreaseConvAlphaMin + (decreaseConvAlphaMax - decreaseConvAlphaMin) &
                        *(MIN(Kp,6._dp))**2/36.
    decreaseConvPsi   = decreaseConvPsiMin + (decreaseConvPsiMax - decreaseConvPsiMin) &
                        *(MIN(Kp,6._dp))**2/36.

!!!!! Recalculate the SCB outerboundary
    if (nIter.ne.0) then
       call Update_Domain(check)
       IF ((iAMR == 1).and.(check)) THEN
          CALL InterpolatePsiR
          CALL mappsi
          CALL psiFunctions
          CALL maptheta
       ENDIF
    endif
!!!!!

    nFail = 0
    SORFail = .false.
    hICalc = .true.
    convDistance = 0.9

    sumb1 = 0._dp
    sumdb1 = 0._dp
    diffmx1 = 0._dp
    entropyFixed = 0._dp
    fluxVolume = 0._dp
  
    iteration = 0
    iCountEntropy = 1
    iconv = 0
    sumdbconv = 0.0_dp
    iConvGlobal = 0

    xStart = x
    yStart = y
    zStart = z
    psiStart = psi
    alphaStart = alfa
    fStart = f

    call system_clock(time1,clock_rate,clock_max)
    starttime=time1/real(clock_rate,dp)

    Outeriters: DO 

       call computeBandJacob
       CALL metrica

       ! Define the right-hand side of the betaEuler equation
       CALL pressure

       normDiffPrev = normDiff
       call compute_convergence
       IF (iteration == 0 .and. isFBDetailNeeded == 1 .and. method /= 3) THEN
          CALL Write_Convergence_Anisotropic('00')
       END IF
       !IF (iteration.ge.1) THEN
       !   if (normDiffPrev.lt.normDiff) exit OuterIters
       !ENDIF

       ! This is a new method for dynamically calculating the blending factors
       ! based on the initial state of the convergence. Note that if the beta
       ! check is changed in ModScbCompute then the ratio of normJxB to
       ! normGradP will also have to be adjusted to lay within similar
       ! boundaries. These blendings are conservative, and could be adjusted to
       ! be more exact, but for now are reasonable -ME
       if (iteration == 0) then
          if (normGradP.gt.80) then
             blendInitial = 0.15
             !convDistance = 0.2
          elseif (normGradP.gt.40) then
             blendInitial = 0.15
          elseif (normGradP.gt.30) then
             blendInitial = 0.15
          elseif (normJxB.ge.2.*normGradP) then
             blendInitial = 0.20
          elseif (normJxB.ge.1.5*normGradP) then
             blendInitial = 0.25
          elseif (normJxB.ge.normGradP) then
             blendInitial = 0.30
          else
             blendInitial = 0.20
          endif
          blendAlpha = blendInitial
          blendPsi = blendInitial

          ! The "distance" is the amount "traveled" between the current solution
          ! and the found solution. The value lies between 0 and 1 with 0 meaning
          ! you completely use the old solution and 1 meaning you completely use
          ! the new solution.
          !! The total "distance" you want to achieve
          !convDistance = 0.9 ! Moved to outside the loop so that it can be
                              ! modified in the loop as needed
          !! The number of iterations with the above blending parameter needed to
          !! achieve the wanted "distance"
          SCBIterNeeded = ceiling(log(1-convDistance)/log(1-blendInitial))
          !! The "distance" currently traveled, this is for when we implement a
          !! relaxation to the blending, but for now should just be set to > 1 to
          !! avoid it impacting the number of iterations performed.
          outDistance = 1.1
       endif

       if (iteration == 0) iteration = 1 ! Iteration 0 saves the full 3D pressure domain

       SCB_CALCULATION: IF (method /= 3) then

          CALL newk
  
          errorAlphaPrev = errorAlpha
  
          IF (iteration==1) alfaSav1 = alfa
          IF (iteration==1) psiSav1 = psi

          IF (MOD(iteration,nrelax) == 0) blendAlpha = relax * blendAlpha
          IF (MOD(iteration,nrelax) == 0) blendPsi = relax * blendPsi

          blendAlpha = MAX(blendAlpha,blendMin)
          blendAlpha = MIN(blendAlpha,blendMax)
          SELECT CASE (method)
            CASE(1)
               CALL directAlpha
            CASE(2)
               CALL iterateAlpha
          END SELECT

          if (SORFail) then
             x    = xStart
             y    = yStart
             z    = zStart
             psi  = psiStart
             alfa = alphaStart
             f = fStart
             hICalc = .false.
             exit OuterIters
          endif

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
  
          xPrev = x
          yPrev = y
          zPrev = z
          alphaPrev = alfa

          Move_points_in_alpha_theta: DO
             ! move zeta grid points along constant alphaEuler and theta lines
             CALL mapalpha
             ! move theta grid points along constant alphaEuler and zeta lines
             CALL maptheta
             CALL metrica
             IF (MINVAL(jacobian(2:nthe-1,2:npsi-1,:)) < 0._dp) THEN
                ! Revert to previous point configuration
                x = xPrev
                y = yPrev
                z = zPrev
                alfa = alphaPrev
                blendAlpha = damp * blendAlpha
                alfa(:,:,:) = alfa(:,:,:)*blendAlpha + (1.-blendAlpha)*alfaSav1(:,:,:)
                call metrica
                if (blendAlpha.lt.blendMin) then
                   !call CON_stop('Failed to converge Alpha potential with minimum blend size')
                   x    = xStart
                   y    = yStart
                   z    = zStart
                   psi  = psiStart
                   alfa = alphaStart
                   f = fStart
                   hICalc = .false.
                   exit OuterIters
                endif
                if (verbose) PRINT*, 'CE: Cycling alpha_theta pts, blendAlpha = ', blendAlpha
                CYCLE Move_points_in_alpha_theta
             END IF
             EXIT Move_points_in_alpha_theta
          END DO Move_points_in_alpha_theta
  
          IF (iAMR == 1) THEN
             CALL InterpolatePsiR
             CALL mappsi
             CALL psifunctions
             CALL maptheta
             psisav1 = psi
          ENDIF

          !IF (isotropy == 1) CALL entropy(entropyFixed, fluxVolume, iCountEntropy)
          IF (isFBDetailNeeded == 1) CALL Write_Convergence_Anisotropic('02')
 
          CALL metric
          IF (MINVAL(jacobian(2:nthe-1,2:npsi-1,:)) < 0._dp) then
             !call CON_STOP('CE: metric problem.')
             x = xStart
             y = yStart
             z = zStart
             psi = psiStart
             alfa = alphaStart
             f = fStart
             hICalc = .false.
             exit OuterIters
          ENDIF

          call computeBandJacob
          !CALL pressure
 
          !c  define the right-hand side of the alphaEuler equation
          CALL newj
  
          errorPsiPrev = errorPsi
          blendPsi = MAX(blendPsi,blendMin)
          blendPsi = MIN(blendPsi,blendMax)
          SELECT CASE (method)
            CASE(1)
              CALL directPsi
            CASE(2)
              CALL iteratePsi
          END SELECT

          if (SORFail) then
             x    = xStart
             y    = yStart
             z    = zStart
             psi  = psiStart
             alfa = alphaStart
             f = fStart
             hICalc = .false.
             exit OuterIters
          endif

          errorPsi = diffmx
          IF (iteration==1) errorPsiPrev = errorPsi
  
          xPrev = x
          yPrev = y
          zPrev = z
          psiPrev = psi

          Move_points_in_psi_theta: DO
             CALL mappsi
             CALL maptheta
             CALL metric
             IF (MINVAL(jacobian(2:nthe-1,2:npsi-1,:)) < 0._dp) THEN
                ! Revert to previous point configuration
                x = xPrev
                y = yPrev
                z = zPrev
                psi = psiPrev
                blendPsi = damp * blendPsi
                psi(:,:,:) = psi(:,:,:)*blendPsi + (1.-blendPsi)*psiSav1(:,:,:)
                if (blendPsi.lt.blendMin) then
                   !call CON_stop('Failed to converge Psi potential with minimum blend size')
                   x    = xStart
                   y    = yStart
                   z    = zStart
                   psi  = psiStart
                   alfa = alphaStart
                   f = fStart
                   hICalc = .false.
                   exit OuterIters
                endif
                if (verbose) PRINT*, 'CE: Cycling psi_theta pts, blendPsi = ', blendPsi
                CYCLE Move_points_in_psi_theta
             END IF
             EXIT Move_points_in_psi_theta
          END DO Move_points_in_psi_theta
  
          IF (iteration == 1) THEN
             errorfirstpsi = sumdb
             diffmxfirstpsi = diffmx
          END IF

          if (verbose) then
             WRITE(*,*) ' itout ',' blendAlpha ',' blendPsi ',' itAlpha ',' diffAlpha ',' errorAlpha ',&
                        ' itPsi ',' diffPsi ',' errorPsi '
             WRITE(*,*) iteration, blendAlpha, blendPsi, nisave1,sumdb1,errorAlpha/twopi_d, &
                        nisave,sumdb,errorPsi/MAXVAL(ABS(psival))
          endif
  
          IF (iAMR == 1) THEN
             CALL InterpolatePsiR
             CALL mappsi
             CALL psiFunctions
             CALL maptheta
             psisav1 = psi
          ENDIF

          ! Need to set an actual convergence criteria using the JxB and GradP
          ! values calculated in the Write_Convergence_Anisotropic subroutine.
          ! For now we will just use the residuals of the SOR calculations
          If ((errorAlpha.lt.decreaseConvAlpha).and.(errorPsi.lt.decreaseConvPsi)) then
             iConvGlobal = 1
          endif

          IF (((iteration.lt.numit).AND.(iConvGlobal.eq.0)) &
              .OR.(iteration.lt.MinSCBIterations) &
              .OR.(iteration.lt.SCBIterNeeded) &
              .OR.(outDistance.lt.convDistance)) THEN
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
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)

       !   The end of the iterative calculation
       PRINT*, iteration, "outer iterations performed in", stoptime-starttime, "seconds."
       IF (boundary /= 'SWMF') PRINT*, "End of calculation."
       PRINT*, ' '
       DEALLOCATE(xPrev, yPrev, zPrev, alphaPrev, psiPrev)
    end if
  
    IF (iteration > numit) lconv = 1
    nitry = nisave

    if (SORFail) then
       x = xStart
       y = yStart
       z = zStart
       psi = psiStart
       alfa = alphaStart
       f = fStart
       hICalc = .false.
    endif
  
    ! The following block should be uncommented for applications where equal-arc-length is needed or desirable
    ! Equal arc-length helps with the h and I integrals -ME
    !constTheta = 0.0_dp ! For equal-arc length
    !chiVal = (thetaVal + constTheta * SIN(2.*thetaVal)) 
    !CALL maptheta
    !CALL metrica
    ! This will compute the new Bfield on the new grid, to get the right pressure mapping

    CALL computeBandJacob
    CALL pressure
    CALL compute_convergence
    CALL entropy(entropyFixed, fluxVolume, iCountEntropy)

    ! Compute physical quantities: currents, field components etc..
    CALL metrics
    ! extrapolate to the fixed boundary
    CALL bounextp

    IF (isFBDetailNeeded == 1 .and. method /= 3) CALL Write_Convergence_Anisotropic('01')
    ! Computes energies and Dst from DPS relation, write to disk (+ Biot-Savart values) 
    ! Remove for speed
    IF (isotropy == 0 .AND. isEnergDetailNeeded == 1) CALL dps_general

    DEALLOCATE(xStart, yStart, zStart, psiStart, alphaStart, fStart)

    RETURN
  
  END SUBROUTINE scb_run
  
!==============================================================================
  SUBROUTINE energy
    !!!! Module Variables
    USE ModScbParams,    ONLY: isotropy
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, dr, dt, dpPrime
    use ModScbVariables, ONLY: pressure3D, x, jacobian, bsq
    !!!! NR Modules
    use nrtype, ONLY: DP

    implicit none
  
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
    use ModScbVariables, ONLY: x, y, z, xx, yy, jacobian, bf, nThetaEquator, &
                               f, fzet, rhoVal, thetaVal, zetaVal, psiVal, &
                               pjconst, r0Start, GradRhoSq, GradThetaSq, GradZetaSq, &
                               GradRhoGradTheta, GradRhoGradZeta, GradThetaGradZeta, &
                               derivXTheta, derivXRho, derivXZeta, &
                               derivYTheta, derivYRho, derivYZeta, &
                               derivZTheta, derivZRho, derivZZeta, &
                               gradRhoX, gradRhoY, gradRhoZ, dPdAlpha, &
                               gradZetaX, gradZetaY, gradZetaZ, dPdPsi, &
                               gradThetaX, gradThetaY, gradThetaZ, pressure3D

    !!!! Module Subroutines/Functions
    use ModRamGSL, ONLY: GSL_Derivs
    !!!! NR Modules
    use nrtype, ONLY: DP

    implicit none
  
    real(DP), INTENT(INOUT) :: ent_local(:,:), vol_local(:,:)
    integer, intent(IN) :: iteration_local
  
    INTEGER :: i, j, k, ierr, idealerr, ncdfId, GSLerr
    REAL(DP) :: yyp, phi, deltaPhi
    ! gradRhoSq, gradRhoGradZeta are global
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dVoldXEq, dVoldYEq, dVoldZeta, dVoldAlpha, &
                                             dVoldRho, dVoldPsi, dEntdXEq, dEntdYEq, &
                                             dEntdZeta, dEntdAlpha, dEntdRho, &
                                             dEntdPsi, facVasGlobal, secondTermB
    REAL(DP) :: delS
    REAL(DP) :: rr1, rr2, zangle, thangle, thangleOnEarth, rr, dza, dya
    REAL(DP) :: dipoleFactor, dipoleFactor4RE, factorIncrease
  
    ALLOCATE(dVoldZeta(npsi,nzeta), dVoldAlpha(npsi,nzeta), &
             dVoldRho(npsi,nzeta), dVoldPsi(npsi,nzeta), &
             dEntdRho(npsi,nzeta), dEntdPsi(npsi,nzeta), &
             dEntdZeta(npsi,nzeta), dEntdAlpha(npsi,nzeta))
    dVoldZeta = 0.0; dVoldAlpha = 0.0; dVoldRho = 0.0; dVoldPsi = 0.0; dEntdRho = 0.0
    dEntdPsi = 0.0; dEntdZeta = 0.0; dEntdAlpha = 0.0

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
  
    CALL GSL_Derivs(rhoVal, zetaVal(1:nzeta), vol_local(:,1:nzeta), &
                       dVoldRho(:,1:nzeta), dVoldZeta(:,1:nzeta), GSLerr)
    CALL GSL_Derivs(rhoVal, zetaVal(1:nzeta), ent_local(:,1:nzeta), &
                       dEntdRho(:,1:nzeta), dEntdZeta(:,1:nzeta), GSLerr)
  
    DO j = 1, npsi
       dVoldPsi(j,:) = 1._dp / f(j) * dVoldRho(j,:)
       dEntdPsi(j,:) = 1._dp / f(j) * dEntdRho(j,:)
    END DO
  
    DO k = 1, nzeta
       dVoldAlpha(:,k) = dVoldZeta(:,k) / fzet(k)
       dEntdAlpha(:,k) = dEntdZeta(:,k) / fzet(k)
    END DO
  
    ALLOCATE(dVoldXEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dVoldYEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dEntdXEq(npsi,nzeta), stat = ierr)
    ALLOCATE(dEntdYEq(npsi,nzeta), stat = ierr)
    dVoldXEq = 0.0; dVoldYEq = 0.0; dEntdXEq = 0.0; dEntdYEq = 0.0

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
  
    ALLOCATE(facVasGlobal(npsi,nzeta), stat = ierr)
    ALLOCATE(secondTermB(npsi,nzeta), stat = ierr)
    facVasGlobal = 0.0; secondTermB = 0.0

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

    DEALLOCATE(facVasGlobal, stat = idealerr)
    DEALLOCATE(secondTermB, stat = idealerr)
  
    DEALLOCATE(dVoldXEq, stat = idealerr)
    DEALLOCATE(dVoldYEq, stat = idealerr)
    DEALLOCATE(dEntdXEq, stat = idealerr)
    DEALLOCATE(dEntdYEq, stat = idealerr)
  
    ! Can de-allocate derivXRho etc.
    DEALLOCATE(dVoldZeta, dVoldAlpha, &
               dVoldRho, dVoldPsi, &
               dEntdRho, dEntdPsi, &
               dEntdZeta, dEntdAlpha)
 
  RETURN
  
  END SUBROUTINE entropy
  
!==============================================================================
  SUBROUTINE bounextp
    !!!! Module Variables
    USE ModScbGrids,     ONLY: nthe, npsi, npsim, nzetap
    use ModScbVariables, ONLY: bf, bsq, phij, bj
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
  
          CALL extap(bf(4,j,k),bf(3,j,k),bf(2,j,k),bf(1,j,k))
          CALL extap(bf(nthe-3,j,k),bf(nthe-2,j,k),bf(nthe-1,j,k) ,bf(nthe,j,k))
          bsq(1,j,k)=bf(1,j,k)**2
          bsq(nthe,j,k)=bf(nthe,j,k)**2
       END DO
       DO i=1,nthe
          CALL extap(bj(i,4,k),bj(i,3,k),bj(i,2,k),bj(i,1,k))
          CALL extap(bj(i,npsi-3,k),bj(i,npsi-2,k),bj(i,npsi-1,k),bj(i,npsi,k))
          
          CALL extap(phij(i,4,k),phij(i,3,k),phij(i,2,k),phij(i,1,k))
          CALL extap(phij(i,npsi-3,k),phij(i,npsi-2,k),phij(i,npsi-1,k),phij(i,npsi,k))
  
          CALL extap(bf(i,4,k),bf(i,3,k),bf(i,2,k),bf(i,1,k))
          CALL extap(bf(i,npsi-3,k),bf(i,npsi-2,k),bf(i,npsi-1,k) ,bf(i,npsi,k))
          IF (bf(i,npsi,k) < 0._dp) bf(i,npsi,k) = bf(i,npsi-1,k)
          bsq(i,1,k) = bf(i,1,k)**2
          bsq(i,npsi,k) = bf(i,npsi,k)**2
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
    use ModScbVariables, ONLY: x, y, z, bsq, jacobian, pnormal, bnormal, &
                               DstBiot, DstBiotInsideGeo, DstDPS, DstDPSInsideGeo, &
                               pper, ppar
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d

    implicit none
  
    INTEGER :: i, j, k, iplx
    REAL(DP) :: magneticEnergyDipole, rsq, totalEnergy, volumeTotal
    REAL(DP), ALLOCATABLE :: magneticEnergy(:), magneticEnergyInsideGeo(:), &
                             thermalEnergy(:), thermalEnergyInsideGeo(:)

    ALLOCATE(magneticEnergy(nthe), magneticEnergyInsideGeo(nthe), &
             thermalEnergy(nthe), thermalEnergyInsideGeo(nthe))
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

    DEALLOCATE(magneticEnergy,magneticEnergyInsideGeo,thermalEnergy,thermalEnergyInsideGeo)  
    RETURN
  
  END SUBROUTINE dps_general
  
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
    use ModScbMain,      ONLY: iCountPressureCall
    use ModScbParams,    ONLY: iLossCone, iOuterMethod, iReduceAnisotropy, isotropy, &
                               isPressureDetailNeeded, PressMode, Isotropic, iSm2
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap, nXRaw, nXRawExt, nYRaw, &
                               nAzimRAM
    use ModScbVariables, ONLY: x, y, z, xx, yy, bf, bsq, rhoVal, zetaVal, thetaVal, &
                               nZetaMidnight, nThetaEquator, pnormal, f, fzet, alfa, &
                               dela, azimRaw, radGrid, angleGrid, ratioEq, dPPerdRho, &
                               dPPerdZeta, dPPerdTheta, dBsqdRho, dBsqdZeta, dPdPsi, &
                               dSqPdPsiSq, dpdAlpha, dSqPdAlphaSq, pressure3D, ppar, &
                               pper, dPperdPsi, bsq, dBBdPsi, dBsqdPsi, dPperdAlpha, &
                               dBBdAlpha, dBsqdAlpha, dBsqdTheta, sigma, tau, &
                               gradRhoGradZeta, gradZetaSq, gradRhoGradTheta, &
                               gradThetaGradZeta, gradRhoSq
    !!!! Module Subroutines/Functions
    USE ModRamGSL,       ONLY: GSL_Derivs, GSL_Interpolation_2D, GSL_Interpolation_1D, &
                               GSL_Smooth_1D
    USE ModSCBIO,        ONLY: write_scb_pressure
    USE ModScbFunctions, ONLY: SavGol7, pRoeRad, extap
    !!!! NR Modules
    use nrtype, ONLY: DP, SP, pi_d, twopi_d

    implicit none

    INTEGER :: i, j, j1, k, k1, ierr, idealerr, GSLerr
    REAL(DP) :: radius, angle, bEqSq, aN, pperN, pparN, yyp, gParam, pEq, ratioB, rBI, &
                colatitudeMid, colatitudeNoo
    REAL(DP) :: aTemp(500), bTemp(500)
    REAL(DP), ALLOCATABLE :: press(:,:), dPresdRho(:,:), dPresdZeta(:,:), &
                             xEq(:,:), yEq(:,:), aratio(:,:), aratioOld(:,:), &
                             aLiemohn(:,:), dSqPresdRhoSq(:,:), dSqPresdZetaSq(:,:), &
                             dSqPresdRhodZeta(:,:), pperEq(:,:), pparEq(:,:), &
                             pperEqOld(:,:), pparEqOld(:,:), radGridEq(:,:), angleGridEq(:,:)

    REAL(DP), ALLOCATABLE :: dipoleFactorMid(:,:), dipoleFactorNoo(:,:)
    REAL(DP), ALLOCATABLE :: BigBracketPsi(:,:,:), BigBracketAlpha(:,:,:), dBBdRho(:,:,:), &
                             dBBdZeta(:,:,:), dummy1(:,:,:), dummy2(:,:,:)
    REAL(DP), ALLOCATABLE :: xRaw(:,:), YRaw(:,:),pressProtonPerRaw(:,:), pressProtonParRaw(:,:), &
                             pressOxygenPerRaw(:,:), pressOxygenParRaw(:,:), pressHeliumPerRaw(:,:), &
                             pressHeliumParRaw(:,:), pressPerRaw(:,:), pressParRaw(:,:), &
                             pressEleParRaw(:,:), pressElePerRaw(:,:), radRaw_local(:), &
                             ratioRaw(:,:)
    REAL(DP), ALLOCATABLE :: pressPerRawExt(:,:), pressParRawExt(:,:), &
                             radRawExt(:), azimRawExt(:)

    iCountPressureCall = iCountPressureCall + 1 ! global variable, counts how many times pressure is called

    ALLOCATE(press(npsi, nzeta+1), dPresdRho(npsi, nzeta+1), dPresdZeta(npsi, nzeta+1), &
             xEq(npsi, nzeta+1), yEq(npsi, nzeta+1), aratio(npsi, nzeta+1), &
             aratioOld(npsi, nzeta+1), aLiemohn(npsi, nzeta+1), dSqPresdRhoSq(npsi,nzeta+1), &
             dSqPresdZetaSq(npsi,nzeta+1), dSqPresdRhodZeta(npsi,nzeta+1), pperEq(npsi,nzeta+1), &
             pparEq(npsi,nzeta+1), pperEqOld(npsi,nzeta+1), pparEqOld(npsi,nzeta+1), &
             radGridEq(npsi, nzeta), angleGridEq(npsi,nzeta))
    ALLOCATE(xRaw(nXRaw,nYRaw), YRaw(nXRaw,nYRaw), pressProtonPerRaw(nXRaw,nYRaw), &
             pressProtonParRaw(nXRaw,nYRaw), pressOxygenPerRaw(nXRaw,nYRaw), &
             pressOxygenParRaw(nXRaw,nYRaw), pressHeliumPerRaw(nXRaw,nYRaw), &
             pressHeliumParRaw(nXRaw,nYRaw), pressPerRaw(nXRaw,nYRaw), pressParRaw(nXRaw,nYRaw), &
             pressEleParRaw(nXRaw,nYRaw), pressElePerRaw(nXRaw,nYRaw), &
             radRaw_local(nXRaw), ratioRaw(nXRaw,nYRaw))
    ALLOCATE(dipoleFactorMid(nthe,npsi),dipoleFactorNoo(nthe,npsi))
    ALLOCATE(pressPerRawExt(nXRawExt,nAzimRAM), pressParRawExt(nXRawExt,nAzimRAM))
    ALLOCATE(radRawExt(nXRawExt), azimRawExt(nAzimRAM))
    pressPerRawExt = 0.0; pressParRawExt = 0.0; radRawExt = 0.0; azimRawExt = 0.0
    press = 0.0; dPresdRho = 0.0; dPresdZeta = 0.0; xEq = 0.0; yEq = 0.0; aratio = 0.0; aratioOld = 0.0
    aLiemohn = 0.0; dSqPresdRhoSq = 0.0; dSqPresdZetaSq = 0.0; dSqPresdRhodZeta = 0.0; pperEq = 0.0
    pparEq = 0.0; pperEqOld = 0.0; pparEqOld = 0.0; radGridEq = 0.0; angleGridEq = 0.0
    xRaw = 0.0; YRaw = 0.0; pressProtonPerRaw = 0.0; pressProtonParRaw = 0.0; pressOxygenPerRaw = 0.0
    pressOxygenParRaw = 0.0; pressHeliumPerRaw = 0.0; pressHeliumParRaw = 0.0; pressPerRaw = 0.0
    pressParRaw = 0.0; pressEleParRaw = 0.0; pressElePerRaw = 0.0; radRaw_local = 0.0; ratioRaw = 0.0
    dipoleFactorMid = 0.0; dipoleFactorNoo = 0.0

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

    ! SCB Grid
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
 
    !if (iteration.eq.0) then
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

       radRawExt(1:nXRaw) = radRaw_local(1:nXRaw)
       DO j1 = nXRaw+1, nXRawExt
          radRawExt(j1) = radRaw_local(nXRaw) + REAL(j1-nXRaw, DP)*(radRaw_local(nXRaw)-radRaw_local(1))/(REAL(nXRaw-1, DP))
       END DO

       azimRawExt(1:nAzimRAM) = azimRaw(1:nYRaw) ! nYRaw = nAzimRAM
       IF (PressMode == 'SKD') then
          pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw) ! from keV/cm^3 to nPa
          pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw) ! from keV/cm^3 to nPa
          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
          DO k1 = 1, nAzimRAM
             DO j1 = nXRaw+1, nXRawExt
                pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) &
                   * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) &
                   /  (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
                pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) &
                   * (89.*EXP(-0.59*radRawExt(j1)) + 8.9*radRawExt(j1)**(-1.53)) &
                   /  (89.*EXP(-0.59*radRawExt(nXRaw)) + 8.9*radRawExt(nXRaw)**(-1.53))
             END DO
          END DO
       ELSEIF (PressMode == 'ROE') then
          pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw) ! from keV/cm^3 to nPa
          pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw) ! from keV/cm^3 to nPa
          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
          DO k1 = 1, nAzimRAM
             DO j1 = nXRaw+1, nXRawExt
                pressPerRawExt(j1,k1) = pressPerRawExt(nXRaw,k1) * pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
                pressParRawExt(j1,k1) = pressParRawExt(nXRaw,k1) * pRoeRad(radRawExt(j1))/pRoeRad(radRawExt(nXRaw))
             END DO
          END DO
       ELSEIF (PressMode == 'EXT') then
          pressPerRaw = 0.16_dp * (pressProtonPerRaw + pressOxygenPerRaw + pressHeliumPerRaw) ! from keV/cm^3 to nPa
          pressParRaw = 0.16_dp * (pressProtonParRaw + pressOxygenParRaw + pressHeliumParRaw) ! from keV/cm^3 to nPa
          pressPerRawExt(1:nXRaw,:) = pressPerRaw(1:nXRaw,:)
          pressParRawExt(1:nXRaw,:) = pressParRaw(1:nXRaw,:)
          DO k1 = 1, nAzimRAM
             DO j1 = nXRaw+1, nXRawExt
                pressPerRawExt(j1,k1) = pressPerRawExt(j1-1,k1) + (radRawExt(j1)-radRawExt(j1-1)) &
                                                                  /(radRawExt(j1-2)-radRawExt(j1-1)) &
                                                                  *(pressPerRawExt(j1-2,k1)-pressPerRawExt(j1-1,k1))
                pressParRawExt(j1,k1) = pressParRawExt(j1-1,k1) + (radRawExt(j1)-radRawExt(j1-1)) &
                                                                  /(radRawExt(j1-2)-radRawExt(j1-1)) &
                                                                  *(pressParRawExt(j1-2,k1)-pressParRawExt(j1-1,k1))
             END DO
          END DO

       ENDIF

       IF (iSm2 == 1) THEN ! Savitzky-Golay smoothing (possibly multiple) for the pressure
          pressPerRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressPerRawExt(1:nXRawExt,1:nAzimRAM))
          pressParRawExt(1:nXRawExt,1:nAzimRAM) = SavGol7(pressParRawExt(1:nXRawExt,1:nAzimRAM))
       ELSEIF (iSm2 == 2) THEN ! B-Spline Fit
          DO j = 1,nXRawExt
             CALL GSL_Smooth_1D(azimRawExt(1:nAzimRAM),pressPerRawExt(j,1:nAzimRAM),aTemp(1:500),bTemp(1:500),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',aTemp(1:500),bTemp(1:500),azimRawExt(1:nAzimRAM),pressPerRawExt(j,1:nAzimRAM),GSLerr)
             CALL GSL_Smooth_1D(azimRawExt(1:nAzimRAM),pressParRawExt(j,1:nAzimRAM),aTemp(1:500),bTemp(1:500),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',aTemp(1:500),bTemp(1:500),azimRawExt(1:nAzimRAM),pressParRawExt(j,1:nAzimRAM),GSLerr)
          ENDDO
          DO k = 1,nAzimRAM
             CALL GSL_Smooth_1D(radRawExt(1:nXRawExt),pressPerRawExt(1:nXRawExt,k),aTemp(1:500),bTemp(1:500),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',aTemp(1:500),bTemp(1:500),radRawExt(1:nXRawExt),pressPerRawExt(1:nXRawExt,k),GSLerr)
             CALL GSL_Smooth_1D(radRawExt(1:nXRawExt),pressParRawExt(1:nXRawExt,k),aTemp(1:500),bTemp(1:500),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',aTemp(1:500),bTemp(1:500),radRawExt(1:nXRawExt),pressParRawExt(1:nXRawExt,k),GSLerr)
          ENDDO
       !ELSEIF (iSm2 == 3) THEN ! Moving Average Filter
       ENDIF
    !endif

    Isotropy_choice:  IF (isotropy == 1) THEN    ! isotropic case
       Isotropic = 'RAM' ! For now just hard code the analytic isotropic pressure for testing -ME
       IF (Isotropic.eq.'RAM') THEN
          !Cubic interpolation
          CALL GSL_Interpolation_2D(radRawExt**2, azimRawExt, pressPerRawExt, &
                                    radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), &
                                    pperEq(1:npsi,2:nzeta), GSLerr)
          CALL GSL_Interpolation_2D(radRawExt**2, azimRawExt, pressParRawExt, &
                                    radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), &
                                    pparEq(1:npsi,2:nzeta), GSLerr)

          ! Sometimes the interpolation can give very small negative values very  near the Earth; inconsequential
          DO k = 1, nzeta
             DO j = 10, 1, -1
                IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2RE from Earth
                   CALL extap(pperEq(j+3,k),pperEq(j+2,k),pperEq(j+1,k),pperEq(j,k))
                   CALL extap(pparEq(j+3,k),pparEq(j+2,k),pparEq(j+1,k),pparEq(j,k))
                END IF
             END DO
          END DO
          WHERE(pperEq < 0.0) pperEq = MINVAL(pressPerRaw) ! 1e-1_dp/pnormal
          WHERE(pparEq < 0.0) pparEq = MINVAL(pressParRaw) ! 1e-1_dp/pnormal

          pperEq(:,nzeta+1) = pperEq(:,2)
          pparEq(:,nzeta+1) = pparEq(:,2)
          pperEq(:,1) = pperEq(:,nzeta)
          pparEq(:,1) = pparEq(:,nzeta)

          pperEq = pperEq/pnormal
          pparEq = pparEq/pnormal

          press = 1.0/3.0 * pparEq + 2.0/3.0 * pperEq

       ELSEIF (Isotropic.eq.'ANA') THEN
       ! Isotropic analytic pressure from TsygMuk and SK
          do j=1,npsi
             do k=1,nzeta
                press(j,k) = pressureTsygMuk(x(nThetaEquator,j,k), y(nThetaEquator,j,k))
             enddo
          enddo
       ENDIF

       press(:,nzetap) = press(:,2)
       press(:,1) = press(:,nzeta)
       CALL GSL_Derivs(rhoVal, zetaVal(1:nzeta), press(:,1:nzeta), &
                       dPresdRho(:,1:nzeta), dPresdZeta(:,1:nzeta), GSLerr)
       dPresdRho(:,nzetap) = dPresdRho(:,2)
       dPresdRho(:,1) = dPresdRho(:,nzeta)
       dPresdZeta(:,nzetap) = dPresdZeta(:,2)
       dPresdZeta(:,1) = dPresdZeta(:,nzeta)
  
       CALL GSL_Derivs(rhoVal, zetaVal(1:nzeta), dPresdRho(:,1:nzeta), &
                       dSqPresdRhoSq(:,1:nzeta), dSqPresdRhodZeta(:,1:nzeta), GSLerr)
       CALL GSL_Derivs(rhoVal, zetaVal(1:nzeta), dPresdZeta(:,1:nzeta), &
                       dSqPresdRhodZeta(:,1:nzeta), dSqPresdZetaSq(:,1:nzeta), GSLerr)
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

       ! These two are just needed for a couple of places
       pper = 2.0/3.0 * pressure3D
       ppar = 1.0/3.0 * pressure3D
    ELSE    ! Anisotropic pressure case
       IF ((boundary.eq.'LANL').or.(boundary.eq.'SWMF')) THEN ! Calculation using RAM pressures

          ! Cubic interpolation
          CALL GSL_Interpolation_2D(radRawExt**2, azimRawExt, pressPerRawExt, &
                                    radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), &
                                    pperEq(1:npsi,2:nzeta), GSLerr)
          CALL GSL_Interpolation_2D(radRawExt**2, azimRawExt, pressParRawExt, &
                                    radGrid(1:npsi,2:nzeta)**2, angleGrid(1:npsi,2:nzeta), &
                                    pparEq(1:npsi,2:nzeta), GSLerr)
  
          ! Sometimes the interpolation can give very small negative values very 
          ! near the Earth; inconsequential
          DO k = 1, nzeta
             DO j = 10, 1, -1
                IF (radGrid(j,k) < 2.0) THEN ! Extrapolation inside 2RE from Earth
                   CALL extap(pperEq(j+3,k),pperEq(j+2,k),pperEq(j+1,k),pperEq(j,k))
                   CALL extap(pparEq(j+3,k),pparEq(j+2,k),pparEq(j+1,k),pparEq(j,k))
                END IF
             END DO
          END DO
          WHERE(pperEq <= 0.0) pperEq = 1e-1_dp/pnormal
          WHERE(pparEq <= 0.0) pparEq = 1e-1_dp/pnormal
          pperEq(:,nzeta+1) = pperEq(:,2)
          pparEq(:,nzeta+1) = pparEq(:,2)
          pperEq(:,1) = pperEq(:,nzeta)
          pparEq(:,1) = pparEq(:,nzeta)

          ! Set pressures to normalized units
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
                   ELSEIF (iLossCone == 1) THEN
                      gParam = 1. / ((1. + aratio(j,k)*(1. - ratioB))**2)
                      ppar(i,j,k) = pEq * 1./(1.+2.*aratio(j,k)/3.) * SQRT(gParam)
                      pper(i,j,k) = pEq * (aratio(j,k)+1.)/(1.+2.*aratio(j,k)/3.) * gParam
                   END IF
                   sigma(i,j,k) = 1._dp + (pper(i,j,k)-ppar(i,j,k)) / bsq(i,j,k)
                   tau(i,j,k) = 1._dp - 2. * (pper(i,j,k) - ppar(i,j,k)) / bsq(i,j,k) * pper(i,j,k)/ppar(i,j,k)
                END DO
                press(j,k) = pEq
             END DO
          END DO
       ELSE     
          STOP 'PROBLEM in pressure.f90'
       END IF

       pressure3D = 1.0/3.0*ppar + 2.0/3.0*pper

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

       CALL GSL_Derivs(thetaVal, rhoVal, zetaVal, pper(1:nthe,1:npsi,1:nzeta), &
                       dPperdTheta, dPperdRho, dPperdZeta, GSLerr)
       CALL GSL_Derivs(thetaVal, rhoVal, zetaVal, bsq(1:nthe,1:npsi,1:nzeta), &
                       dBsqdTheta, dBsqdRho, dBsqdZeta, GSLerr)
  
       DO j = 1, npsi
          !DO k = 1, nzeta
          !   if (SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) > 7._dp) then
          !      dPPerdRho(:,j,k) = dPPerdRho(:,j-1,k)
          !   endif
          !ENDDO
          dPperdPsi(:,j,:) = 1./f(j) * dPperdRho(:,j,:)
          IF (iOuterMethod == 2) dBBdPsi(:,j,:) = dBBdRho(:,j,:) / f(j)
          dBsqdPsi(:,j,:) = 1./f(j) * dBsqdRho(:,j,:)
       END DO
  
       DO k = 1, nzeta
          !DO j = 1, npsi
          !   if (SQRT(x(nThetaEquator,j,k)**2+y(nThetaEquator,j,k)**2) > 7._dp) then
          !      dPPerdZeta(:,j,k) = dPPerdZeta(:,j-1,k)
          !   endif
          !END DO
          dPperdAlpha(:,:,k) = 1. / fzet(k) * dPperdZeta(:,:,k)
          IF (iOuterMethod == 2) dBBdAlpha(:,:,k) = dBBdZeta(:,:,k) / fzet(k)
          dBsqdAlpha(:,:,k) = 1. / fzet(k) * dBsqdZeta(:,:,k)
       END DO
  
       IF (iOuterMethod == 2) THEN ! If using the Newton method, need these
          ALLOCATE(BigBracketPsi(nthe,npsi,nzeta), stat = ierr)
          ALLOCATE(BigBracketAlpha(nthe,npsi,nzeta), stat = ierr)
          ALLOCATE(dBBdRho(nthe,npsi,nzeta), stat = ierr)
          ALLOCATE(dBBdZeta(nthe,npsi,nzeta), stat = ierr)
          ALLOCATE(dummy1(nthe,npsi,nzeta), stat = ierr)
          ALLOCATE(dummy2(nthe,npsi,nzeta), stat = ierr)
          BigBracketPsi = 0.0; BigBracketAlpha = 0.0; dBBdRho = 0.0
          dBBdZeta = 0.0; dummy1 = 0.0; dummy2 = 0.0
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
          CALL GSL_Derivs(thetaVal, rhoVal, zetaVal, &
                          BigBracketAlpha(1:nthe,1:npsi,1:nzeta), &
                          dummy1, dummy2, dBBdZeta, GSLerr)
          CALL GSL_Derivs(thetaVal, rhoVal, zetaVal, &
                          BigBracketPsi(1:nthe,1:npsi,1:nzeta), &
                          dummy1, dBBdRho, dummy2, GSLerr)

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

    IF (iteration==0.and.isPressureDetailNeeded==1) call write_scb_pressure

    DEALLOCATE(press, dPresdRho, dPresdZeta, xEq, yEq, aratio, aratioOld, &
               aLiemohn, dSqPresdRhoSq, dSqPresdZetaSq, dSqPresdRhodZeta, pperEq, &
               pparEq, pperEqOld, pparEqOld, radGridEq, angleGridEq)
    DEALLOCATE(xRaw, YRaw, pressProtonPerRaw, pressProtonParRaw, pressOxygenPerRaw, &
               pressOxygenParRaw, pressHeliumPerRaw, pressHeliumParRaw, pressPerRaw, &
               pressParRaw, pressEleParRaw, pressElePerRaw, radRaw_local, ratioRaw)
    DEALLOCATE(dipoleFactorMid, dipoleFactorNoo)
    DEALLOCATE(pressPerRawExt, pressParRawExt, radRawExt, azimRawExt)

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
  use ModIOUnit,       ONLY: UNITTMP_

  implicit none

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

  phi = ATAN2(yEqGsm, xEqGsm) ! azimuthal angle, atan(y/x)

  PSWStar = pdynGlobal/3._dp

  theta = ATAN2(bzimfGlobal,byimfGlobal) ! IMF clock angle
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


  implicit none

  integer, parameter :: iCorrectedPressure = 1
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
