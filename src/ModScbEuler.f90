!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbEuler
  ! Contains subroutines responsible for calculating the alpha (actually beta)
  ! and psi (actually alpha) Euler potentials
  
  use ModScbVariables, ONLY: psi, alfa, &
                             alphaVal, alphaValInitial, psiVal, blendAlpha, &
                             blendPsi, diffmx, rjac, decreaseConvAlpha, &
                             decreaseConvPsi, errorAlpha, errorAlphaPrev, &
                             errorPsi, errorPsiPrev, iAlphaMove, iPsiMove
  
  use ModScbFunctions, ONLY: extap
  

  implicit none
  
  contains

!==============================================================================
  SUBROUTINE mapTheta
    !!!! Module Variables
    USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: diffmx, rjac, nisave,  x, y, z, sumb, sumdb, chiVal

    !!!! Module Subroutines/Functions
    use ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d


    implicit none

    INTEGER :: i, j, k, i1, i2, GSLerr
    REAL(DP), DIMENSION(nthe) :: xOld, yOld, zOld, distance, chiValOld
    psiChange = 0
    theChange = 0

    ! Now move theta coordinates along each surface equal arc length along the i grids
    zetaloop: DO k = 2, nzeta
       fluxloop: DO j = 1+psiChange, npsi-psiChange
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

          i1 = 1 + theChange
          i2 = nthe - theChange
          CALL GSL_Interpolation_1D('Cubic',chiValOld,xOld,chiVal(i1:i2),x(i1:i2,j,k),GSLerr)
          CALL GSL_Interpolation_1D('Cubic',chiValOld,yOld,chiVal(i1:i2),y(i1:i2,j,k),GSLerr)
          CALL GSL_Interpolation_1D('Cubic',chiValOld,zOld,chiVal(i1:i2),z(i1:i2,j,k),GSLerr)
       END DO fluxloop
    END DO zetaloop

    !  periodic boundary conditions
    x(:,:,1) = x(:,:,nzeta)
    y(:,:,1) = y(:,:,nzeta)
    z(:,:,1) = z(:,:,nzeta)
    x(:,:,nzetap) = x(:,:,2)
    y(:,:,nzetap) = y(:,:,2)
    z(:,:,nzetap) = z(:,:,2)

    RETURN

  END SUBROUTINE mapTheta
  
!================================================!
!=========== Alpha Euler Potential ==============!
!================================================!
  SUBROUTINE alphaFunctions

    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: nzeta
    use ModScbVariables, ONLY: zetaVal, fzet, fzetp, alphaVal

    USE ModRamGSL, ONLY: GSL_Interpolation_1D


    implicit none

    integer :: GSLerr
    REAL(DP) :: alphaValue(nzeta)

    CALL GSL_Interpolation_1D('Cubic',zetaVal(1:nzeta), alphaval(1:nzeta), &
                              alphaValue, fzet(1:nzeta), GSLerr)
    CALL GSL_Interpolation_1D('Cubic',zetaVal(1:nzeta), fzet(1:nzeta), &
                              alphaValue, fzetp(1:nzeta), GSLerr)

    RETURN
  END SUBROUTINE alphaFunctions

!==============================================================================
  SUBROUTINE InterpolateAlphaPhi
    ! Interpolates the new values of psi at the new locations xnew on midnight 
    ! equator
    !!!! Module Variables
    USE ModScbParams,    ONLY: iAzimOffset
    USE ModScbGrids, ONLY: nzeta, npsi
    use ModScbVariables, ONLY: x, y, nThetaEquator, alphaVal, zetaVal
    !!!! Module Subroutines/Functions
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP, twopi_d


    implicit none

    REAL(DP), DIMENSION(nzeta+1) :: phiEqmid, alphaval1D
    INTEGER :: ialloc, ierr, k, j, jmax, GSLerr

    REAL(DP) :: distConsecFluxSqOld, distConsecFluxSq

    distConsecFluxSq = 0._dp
    distConsecFluxSqOld = 0._dp

    IF (iAzimOffset == 2) THEN
       DO j = 2, npsi
          phiEqmid = atan2(y(nThetaEquator,j,:),x(nThetaEquator,j,:))
          where (phiEqMid(3:nzeta).lt.0) phiEqMid(3:nzeta) = phiEqMid(3:nzeta) + twopi_d
          phiEqMid(1) = phiEqMid(nzeta) - twopi_d
          phiEqMid(nzeta+1) = phiEqMid(2) + twopi_d
          do k = 2, nzeta
             distConsecFluxSq = phiEqMid(k) - phiEqMid(k-1)
             IF (distConsecFluxSq > distConsecFluxSqOld) THEN
                distConsecFluxSqOld = distConsecFluxSq
                jMax = j
             END IF
          END DO
       end DO
    ELSE IF (iAzimOffset == 1) THEN
       jmax = npsi
    END IF

    phiEqmid = atan2(y(nThetaEquator,jMax,:),x(nThetaEquator,jMax,:))
    where (phiEqMid(3:nzeta).lt.0) phiEqMid(3:nzeta) = phiEqMid(3:nzeta) + twopi_d
    phiEqMid(1) = phiEqMid(nzeta) - twopi_d
    phiEqMid(nzeta+1) = phiEqMid(2) + twopi_d

    alphaval1D = alphaval

    CALL GSL_Interpolation_1D('Cubic',phiEqMid,alphaVal1D,zetaVal(1:nzeta),alphaVal(1:nzeta),GSLerr)

    alphaVal(1) = alphaVal(nzeta) - twopi_d
    alphaVal(nzeta+1) = alphaval(2) + twopi_d
    
    RETURN

  END SUBROUTINE InterpolateAlphaPhi

!==============================================================================
  SUBROUTINE alfges
    !   initial guess of alpha 

    use ModScbMain,  ONLY: DP  
    use ModScbGrids, ONLY: nthe, npsi, nzeta, nzetap, nthem
 
    use nrtype, ONLY: pi_d


    implicit none

    INTEGER :: i, j, k
  
    DO k = 1, nzetap
       alfa(1:nthe,1:npsi,k) = alphaval(k)
    END DO

  END SUBROUTINE alfges
                               
!==============================================================================
  SUBROUTINE mapAlpha(iSmoothMove)
    ! new cubic GSL interpolation, without involving linear distance calculation
    USE ModScbMain,      ONLY: DP
    USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, ny, nthem, nzetap
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, alfaPrev, &
                               left, right
  
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
  
    use nrtype, ONLY: pi_d
  

    implicit none
  
    INTEGER, INTENT(IN) :: iSmoothMove
    REAL(DP), DIMENSION(nzeta+1) :: xOld, yOld, zOld, alfaOld
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
    INTEGER :: i, j, k, ierr, GSLerr
    psiChange = 0
    theChange = 0

    blend = 0.1_dp**iAlphaMove
  
    IF (iSmoothMove /= 0 .AND. iAlphaMove > 1) THEN
       ! Add these in difficult equilibria
       PRINT*, 'mapalpha: blend = ', blend
       x = 1*blend*x + (1.-1*blend)*xPrev
       y = 1*blend*y + (1.-1*blend)*yPrev
       z = 1*blend*z + (1.-1*blend)*zPrev
    ELSE
       xPrev = x
       yPrev = y
       zPrev = z
  
       ierr = 0

       jloop : DO j = 1+psiChange, npsi-psiChange
          iloop: DO i = 1+theChange, nthe-theChange
             xOld(1:nzeta+1) = x(i,j,1:nzeta+1)
             yOld(1:nzeta+1) = y(i,j,1:nzeta+1)
             zOld(1:nzeta+1) = z(i,j,1:nzeta+1)
             alfaOld(1:nzeta+1) = alfa(i,j,1:nzeta+1)
             !DO k = 2,nzeta+1
             !   if (alfaOld(k).lt.alfaOld(k-1)) then
             !      alfaOld(k) = alfaOld(k-1)+1E-6
             !   endif
             !ENDDO

             CALL GSL_Interpolation_1D('Cubic',alfaOld,xOld,alphaVal(2:nzeta),x(i,j,2:nzeta),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',alfaOld,yOld,alphaVal(2:nzeta),y(i,j,2:nzeta),GSLerr)
             CALL GSL_Interpolation_1D('Cubic',alfaOld,zOld,alphaVal(2:nzeta),z(i,j,2:nzeta),GSLerr)
          END DO iloop
       END DO jloop

       ! Periodic boundary conditions
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzetap) = x(:,:,2)
       y(:,:,nzetap) = y(:,:,2)
       z(:,:,nzetap) = z(:,:,2)

       call alfges
    END IF
  
    RETURN
  
  END SUBROUTINE mapAlpha
  
  !------------------------------------------------------------------------------
  SUBROUTINE directAlpha

    CALL CON_STOP('Direct matric inversion is currently not working')  

    RETURN
  
  END SUBROUTINE directAlpha
  
!==============================================================================
  SUBROUTINE iterateAlpha
    !   iteration procedure to obtain matrix inversion for alfa Eqn.

    use ModRamParams,    ONLY: verbose  
    USE ModScbMain,      ONLY: nimax
    use ModScbParams,    ONLY: isSORDetailNeeded, iWantAlphaExtrapolation, &
                               InConAlpha
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecx, &
                               left, right, SORFail
  
    use nrtype, ONLY: DP, pi_d
  

    implicit none
 
    REAL(DP) :: anorm, diff, rjac, &
         & anormaverage, dyDummy, anormResid, anormError, anormf
    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:), alfaPrevTemp(:,:,:), om(:)
    REAL(DP) :: RESULT, omegaOpt
    REAL(DP), ALLOCATABLE :: angle(:,:,:), resid(:,:,:)
    INTEGER, ALLOCATABLE :: ni(:)
    INTEGER :: j, i, ict, jz, izmx, jzmx, kmx, ierr, nratio, &
         myPsiBegin, myPsiEnd, my_array_type2, psiRangeDiff, resultInt, loc(3)
    INTEGER, SAVE :: k, kp, km, iz, im, ip
    !$OMP THREADPRIVATE(k,kp,km,iz,im,ip)
  
    ALLOCATE(angle(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaprev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(ni(npsi), om(ny))
    angle = 0.0; alfaprev = 0.0; alfaPrevTemp = 0.0; resid = 0.0; ni = 0.0; om = 0.0

    DO k = 1, nzeta
       DO j = 1, npsi
          DO i = 1, nthe
             angle(i,j,k) = ACOS(x(i,j,k)/SQRT(x(i,j,k)**2+y(i,j,k)**2))  - pi_d
          END DO
       END DO
    END DO
  
    rjac = 1._dp - 2._dp*pi_d*pi_d / (REAL(nzeta,dp)*REAL(nzeta,dp) + REAL(nthe,dp)*REAL(nthe,dp))  ! Radius of convergence for the Jacobi method, can be used
    ! to find optimal SOR omega
  
    om = 1._dp
    !omegaOpt = 2._dp / (1._dp + pi_d/REAL(nzeta+nthe, DP))
    omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))
  
    alfaPrev = alfa
    ni = 0
    resid = 0

!$OMP PARALLEL DO
    psiloop: DO  jz = left+1, right-1
       ni(jz) = 1
  
       Iterations: DO WHILE (ni(jz) <= nimax)
          zetaloop: DO k = 2, nzeta
             kp = k + 1
             km = k - 1
             thetaloop: DO  iz = 2, nthe-1
                im = iz - 1
                ip = iz + 1
                ! Use natural rowwise ordering
                resid(iz,jz,k) = - vecd(iz,jz,k)*alfa(iz,jz,k)  &
                                 + vec1(iz,jz,k)*alfa(im,jz,km) &
                                 + vec2(iz,jz,k)*alfa(iz,jz,km) &
                                 + vec3(iz,jz,k)*alfa(ip,jz,km) &
                                 + vec4(iz,jz,k)*alfa(im,jz,k)  &
                                 + vec6(iz,jz,k)*alfa(ip,jz,k)  &
                                 + vec7(iz,jz,k)*alfa(im,jz,kp) &
                                 + vec8(iz,jz,k)*alfa(iz,jz,kp) &
                                 + vec9(iz,jz,k)*alfa(ip,jz,kp) &
                                 - vecx(iz,jz,k)
                alfa(iz,jz,k) = alfa(iz,jz,k) + om(jz) * (resid(iz,jz,k) / vecd(iz,jz,k))
                if (alfa(iz,jz,k).ne.alfa(iz,jz,k)) then
                   if (verbose) write(*,*) iz,jz,k,alfa(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'NaN encountered in ModScbEuler iterateAlpha'
                   alfa(iz,jz,k) = alfaprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                elseif (alfa(iz,jz,k)+1.0.eq.alfa(iz,jz,k)) then
                   if (verbose) write(*,*) iz,jz,k,alfa(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'Large number encountered in ModScbEuler iterateAlpha'
                   alfa(iz,jz,k) = alfaprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                endif
             END DO thetaloop
          END DO zetaloop

          om(jz) = omegaOpt
  
          ! Check for inner loop convergence of Alpha
          IF (maxval(abs(resid(2:nthem,jz,2:nzeta))).lt.InConAlpha) then
             EXIT Iterations
          END IF
  
          ni(jz)  =   ni(jz) + 1
  
       END DO Iterations
    END DO Psiloop
!$OMP END PARALLEL DO

    nisave = maxval(ni)
    sumdb = SUM(ABS(alfa(2:nthem,2:npsi-1,2:nzeta) - alfaprev(2:nthem,2:npsi-1,2:nzeta)))
    sumb = SUM(ABS(alfa(2:nthem,2:npsi-1,2:nzeta)))
    diffmx = maxval(abs(resid(2:nthem,2:npsi-1,2:nzeta)))

    !...  set "blending" in alpha for outer iteration loop
    DO j = 1, npsi
       DO i = 1, nthe
          DO k = 2, nzeta
             alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1._dp - blendAlpha)
          END DO
          alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
          alfa(i,j,nzetap) = alfa(i,j,2) + 2._dp * pi_d
       END DO
    END DO

    DO k = 2, nzeta
       DO i = 2,nthe-1
          CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsi-1,k),alfa(i,npsi,k))
          CALL extap(alfa(i,4,k),alfa(i,3,k),alfa(i,2,k),alfa(i,1,k))
       ENDDO
       !DO j = 1,npsi
       !   CALL extap(alfa(nthe-3,j,k),alfa(nthe-2,j,k),alfa(nthe-1,j,k),alfa(nthe,j,k))
       !   CALL extap(alfa(4,j,k),alfa(3,j,k),alfa(2,j,k),alfa(1,j,k))
       !ENDDO
    ENDDO

    IF (ALLOCATED(angle)) DEALLOCATE(angle, STAT = ierr)
    IF (ALLOCATED(alfaPrev)) DEALLOCATE(alfaPrev)
    IF (ALLOCATED(alfaPrevTemp)) DEALLOCATE(alfaPrevTemp)
    IF (ALLOCATED(resid)) DEALLOCATE(resid)
    DEALLOCATE(ni,om)
 
    RETURN
  
  END SUBROUTINE iterateAlpha
  
!================================================!
!============ Psi Euler Potential ===============!
!================================================!
  SUBROUTINE psiFunctions
  
    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: npsi
    use ModScbVariables, ONLY: rhoVal, f, fp
 
    use ModRamGSL, ONLY: GSL_Derivs 
  

    implicit none
  
    REAL(DP) :: psiValue(npsi)
    INTEGER :: j, GSLerr
  
    CALL GSL_Derivs(rhoVal(1:npsi), psival(1:npsi), f(1:npsi), GSLerr)
    CALL GSL_Derivs(rhoVal(1:npsi), f(1:npsi), fp(1:npsi), GSLerr)
  
    RETURN
  END SUBROUTINE psiFunctions

!==============================================================================
  SUBROUTINE InterpolatePsiR
    ! Interpolates the new values of psi at the new locations xnew on midnight equator
    !!!! Module Variables
    USE ModScbParams,    ONLY: iAzimOffset, psiChange
    USE ModScbGrids, ONLY: npsi, nzeta
    use ModScbVariables, ONLY: x, y, nThetaEquator, nZetaMidnight, psiVal, kmax, &
                               radEqmidNew
    !!!! Module Subroutines/Functions
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP


    implicit none

    REAL(DP), DIMENSION(npsi) :: radEqmid, psival1D
    INTEGER :: ialloc, ierr, j, k, GSLerr

    REAL(DP) :: deltaR, distConsecFluxSqOld, distConsecFluxSq
    REAL(DP) :: radius(npsi)

    distConsecFluxSq = 0._dp
    distConsecFluxSqOld = 0._dp

    IF (iAzimOffset == 2) THEN
       DO k = 2, nzeta
          distConsecFluxSq = x(nThetaEquator,npsi-psiChange,k)**2 &
                           + y(nThetaEquator,npsi-psiChange,k)**2
          IF (distConsecFluxSq > distConsecFluxSqOld) THEN
             distConsecFluxSqOld = distConsecFluxSq
             kMax = k 
          END IF
       end DO
    ELSE IF (iAzimOffset == 1) THEN
       kmax = kmax
    END IF

    DO j = 1, npsi
       radius(j) = SQRT(x(nThetaEquator,j,kMax)**2 + y(nThetaEquator,j,kMax)**2)
    END DO

    deltaR = (radius(npsi) - radius(1)) / REAL(npsi-1,dp)
    radEqMidNew(1) = radius(1)
    DO j = 2, npsi
       radEqMidNew(j) = radius(1) + REAL(j-1,dp)*deltaR
    END DO

    radEqmid = SQRT(x(nThetaEquator,:,kMax)**2 + y(nThetaEquator,:,kMax)**2)

    !do j = 2, npsi
    !   if (radEqmid(j).lt.radEqMid(j-1)) radEqMid(j) = radEqMid(j-1) + 1E-6
    !enddo

    psival1D = psival
    CALL GSL_Interpolation_1D('Cubic',radEqMid, psiVal1D, radEqMidNew(2:npsi), psiVal(2:npsi), GSLerr)

    RETURN

  END SUBROUTINE InterpolatePsiR

!==============================================================================
  SUBROUTINE psiges
    !   initial guess of poloidal flux
  
    USE ModScbGrids,     ONLY: npsi, nthe, nzetap
  
    INTEGER :: i, j, k
  
    !     initial guess for psi
    DO j=1,npsi
       psi(:,j,:)=psival(j)
    END DO
  
    RETURN
  
  END SUBROUTINE psiges
  
!==============================================================================
  SUBROUTINE mapPsi
    ! new psiMap, without computing linear distance
    ! first and last flux surface remain unchanged
    USE ModScbMain,      ONLY: DP
    USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, left, right
 
    use ModRamGSL, ONLY: GSL_Interpolation_1D 
  

    implicit none
  
    INTEGER :: iSmoothMove
    INTEGER :: ierr, k, i, j, GSLerr, i1, i2
    REAL(DP), DIMENSION(npsi) :: xOld, yOld, zOld, psiOld
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
    psiChange = 0
    theChange = 0

    iSmoothMove = 0
    blend = 0.1_dp**iPsiMove
    IF (iSmoothMove /= 0 .AND. iPsiMove > 1) THEN
       ! Add these in difficult equilibria
       PRINT*, 'mappsi: blend = ', blend
       x = 1*blend*x + (1.-1*blend)*xPrev
       y = 1*blend*y + (1.-1*blend)*yPrev
       z = 1*blend*z + (1.-1*blend)*zPrev
       !C z(nThetaEquator,:,:) = 0._dp ! Symmetry
    ELSE
       xPrev = x
       yPrev = y
       zPrev = z
       ierr = 0
  
       kloop: DO k = 2, nzeta
          iloop: DO i = 1+theChange, nthe-theChange
             xOld(1:npsi) = x(i,1:npsi,k)
             yOld(1:npsi) = y(i,1:npsi,k)
             zOld(1:npsi) = z(i,1:npsi,k)
             psiOld(1:npsi) = psi(i,1:npsi,k)
             !do j = 2,npsi
             !   if (psiOld(j).lt.psiOld(j-1)) then
             !      psiold(j) = psiOld(j-1) + 1E-6
             !   endif
             !enddo

             i1 = 1 + psiChange
             i2 = npsi - psiChange
             CALL GSL_Interpolation_1D('Cubic',psiOld, xOld, psiVal(i1:i2), x(i,i1:i2,k), GSLerr)
             CALL GSL_Interpolation_1D('Cubic',psiOld, yOld, psiVal(i1:i2), y(i,i1:i2,k), GSLerr)
             CALL GSL_Interpolation_1D('Cubic',psiOld, zOld, psiVal(i1:i2), z(i,i1:i2,k), GSLerr)
          END DO iloop
       END DO kloop
  
       ! periodic boundary condition in zeta
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzetap) = x(:,:,2)
       y(:,:,nzetap) = y(:,:,2)
       z(:,:,nzetap) = z(:,:,2)
  
       call psiges
    END IF
  
    RETURN
  
  END SUBROUTINE mapPsi
  
!==============================================================================
  SUBROUTINE directPsi

    CALL CON_Stop('Direct matrix inversion is not currently working')

    RETURN
  
  END SUBROUTINE directPsi
  
!==============================================================================
  SUBROUTINE iteratePsi
    !   iteration procedure to obtain matrix inversion

    use ModRamParams,    ONLY: verbose  
    USE ModScbMain,      ONLY: DP, nimax
    use ModScbParams,    ONLY: isSORDetailNeeded, InConPsi
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecr, &
                               left, right, SORFail
  
    use nrtype, ONLY: pi_d
  

    implicit none
 
    REAL(DP) :: omegaOpt
    REAL(DP) :: omc, anorm, anormf, diff, ano, sumbtest, anormResid, anormError
    REAL(DP), ALLOCATABLE :: psiPrev(:,:,:), psiPrevTemp(:,:,:), resid(:,:,:), om(:)
    INTEGER, ALLOCATABLE :: ni(:)
    INTEGER :: j, k, i, ierr, ict, nratio, myAlphaBegin, &
               myAlphaEnd, my_array_type_psi2, alphaRangeDiff, resultInt, loc(3)
    INTEGER, SAVE :: jz, jp, jm, iz, im, ip
    !$OMP THREADPRIVATE(jz,jp,jm,iz,im,ip)

    !     perform SOR iteration
    !..   choose for inner iteration loop
  
    ALLOCATE(psiPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(psiPrev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(ni(nzeta), om(na))
    psiPrevTemp = 0.0; psiPrev = 0.0; resid = 0.0; ni = 0.0; om = 0.0

    rjac = 1._dp - 2._dp*pi_d*pi_d / (REAL(nthe,dp)*REAL(nthe,dp) + REAL(npsi,dp)*REAL(npsi,dp)) ! Radius of conv. of Jacobi iteration,
    ! could be used to find omega optimal in SOR
    !rjac = (cos(pi_d/real(npsi,dp)) + (dr/dt)**2 * cos(pi_d/real(nthe,dp))) / &
    !     (1._dp +  (dr/dt)**2)
  
    omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))
    om = 1.0
    psiPrev = psi
    ni = 0
    resid = 0

!$OMP PARALLEL DO
    alphaLoop: DO k = 2, nzeta
       ni(k) = 1
  
       Iterations: DO WHILE (ni(k) <= nimax)
          !anormResid = 0._dp
          jLoop: DO jz = left+1, right-1
             jp = jz + 1
             jm = jz - 1
             iLoop: DO  iz = 2, nthe-1
                im = iz - 1
                ip = iz + 1
                resid(iz,jz,k) = &
                        - vecd(iz,jz,k)*psi(iz,jz,k) &
                        + vec1(iz,jz,k)*psi(im,jm,k) &
                        + vec2(iz,jz,k)*psi(iz,jm,k) &
                        + vec3(iz,jz,k)*psi(ip,jm,k) &
                        + vec4(iz,jz,k)*psi(im,jz,k) &
                        + vec6(iz,jz,k)*psi(ip,jz,k) &
                        + vec7(iz,jz,k)*psi(im,jp,k) &
                        + vec8(iz,jz,k)*psi(iz,jp,k) &
                        + vec9(iz,jz,k)*psi(ip,jp,k) &
                        - vecr(iz,jz,k)
                psi(iz,jz,k) = psi(iz,jz,k) + om(k)*resid(iz,jz,k)/vecd(iz,jz,k)
                if (psi(iz,jz,k).ne.psi(iz,jz,k)) then
                   if (verbose) write(*,*) iz,jz,k,psi(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'NaN encountered in ModScbEuler iteratePsi'
                   psi(iz,jz,k) = psiprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                elseif (psi(iz,jz,k)+1.0.eq.psi(iz,jz,k)) then
                   if (verbose) write(*,*) iz,jz,k,psi(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'Large number encountered in ModScbEuler iteratePsi'
                   psi(iz,jz,k) = psiprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                endif
             END DO iLoop
          END DO jLoop

          om(k) = omegaOpt
  
          ! Check for inner loop Psi convergence
          IF (maxval(abs(resid(2:nthem,2:npsi-1,k))).lt.InConPsi) then
             EXIT Iterations
          END IF
  
          ni(k) = ni(k) + 1
  
       END DO Iterations
    END DO AlphaLoop
!$OMP END PARALLEL DO

    nisave = maxval(ni)
    sumdb = SUM(ABS(psi(2:nthem,2:npsim,2:nzeta) - psiprev(2:nthem,2:npsim,2:nzeta)))
    sumb = SUM(ABS(psi(2:nthem,2:npsim,2:nzeta)))
    diffmx = maxval(abs(resid(2:nthem,2:npsi-1,2:nzeta)))

    ! Set blend for outer iteration loop, and apply periodic boundary conditions
    DO j = 1, npsi
       DO i = 1, nthe
          DO k = 2,nzeta
             psi(i,j,k) = psi(i,j,k) * blendPsi + psival(j) * (1._dp - blendPsi)
             !if ((j.gt.2).and.(psi(i,j,k).lt.psi(i,j-1,k))) then
             !   psi(i,j,k) = psi(i,j-1,k) + 1E-6
             !endif
          END DO
       END DO
    END DO

    DO k = 2, nzeta
       DO i = 2,nthe-1
          CALL extap(psi(i,npsi-3,k),psi(i,npsi-2,k),psi(i,npsi-1,k),psi(i,npsi,k))
          CALL extap(psi(i,4,k),psi(i,3,k),psi(i,2,k),psi(i,1,k))
       ENDDO
       !DO j = 1,npsi
       !   CALL extap(psi(nthe-3,j,k),psi(nthe-2,j,k),psi(nthe-1,j,k),psi(nthe,j,k))
       !   CALL extap(psi(4,j,k),psi(3,j,k),psi(2,j,k),psi(1,j,k))
       !ENDDO
    ENDDO

    psi(:,:,1) = psi(:,:,nzeta)
    psi(:,:,nzetap) = psi(:,:,2)


    IF (ALLOCATED(psiPrevTemp))  DEALLOCATE(psiPrevTemp)
    IF (ALLOCATED(psiPrev)) DEALLOCATE(psiPrev)
    IF (ALLOCATED(resid)) DEALLOCATE(resid)
    DEALLOCATE(ni, om)

    RETURN
  
  END SUBROUTINE iteratePsi
  
  END MODULE ModScbEuler
