!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbEuler
  ! Contains subroutines responsible for calculating the alpha (actually beta)
  ! and psi (actually alpha) Euler potentials
  
  implicit none
  
  contains

!==============================================================================
  SUBROUTINE mapTheta
    !!!! Module Variables
    USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: x, y, z, chiVal

    !!!! Module Subroutines/Functions
    use ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d

    implicit none

    INTEGER :: i, j, k, i1, i2, GSLerr
    REAL(DP), ALLOCATABLE :: xOld(:), yOld(:), zOld(:), distance(:), chiValOld(:)

    ALLOCATE(xOld(nthe),yOld(nthe),zOld(nthe),distance(nthe),chiValOld(nthe))
    xOld = 0.0; yOld = 0.0; zOld = 0.0; distance = 0.0; chiValOld = 0.0

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
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)

    DEALLOCATE(xOld,yOld,zOld,distance,chiValOld)
    RETURN

  END SUBROUTINE mapTheta
  
!================================================!
!=========== Alpha Euler Potential ==============!
!================================================!
!==============================================================================
  SUBROUTINE alfges
    !   initial guess of alpha 

    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables, ONLY: alfa, alphaVal
 
    implicit none

    INTEGER :: k
  
    DO k = 1, nzeta+1
       alfa(1:nthe,1:npsi,k) = alphaval(k)
    END DO

  END SUBROUTINE alfges
                               
!==============================================================================
  SUBROUTINE mapAlpha
    ! new cubic GSL interpolation, without involving linear distance calculation
    USE ModScbMain,      ONLY: DP
    USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: x, y, z, alfa, alphaVal
  
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
  
    implicit none
  
    INTEGER :: i, j, k, GSLerr, ii, jj
    REAL(DP) :: temp
    REAL(DP), ALLOCATABLE :: xOld(:), yOld(:), zOld(:), alfaOld(:)

    ALLOCATE(xOld(nzeta+1),yOld(nzeta+1),zOld(nzeta+1),alfaOld(nzeta+1))
    xOld = 0.0; yOld = 0.0; zOld = 0.0; alfaOld = 0.0

    jloop : DO j = 1+psiChange, npsi-psiChange
       iloop: DO i = 1+theChange, nthe-theChange
          xOld(1:nzeta+1) = x(i,j,1:nzeta+1)
          yOld(1:nzeta+1) = y(i,j,1:nzeta+1)
          zOld(1:nzeta+1) = z(i,j,1:nzeta+1)
          alfaOld(1:nzeta+1) = alfa(i,j,1:nzeta+1)
          !do ii = 1, nzeta+1
          !   do jj = ii, nzeta+1
          !      if (alfaOld(ii) > alfaOld(jj)) then
          !         temp = alfaOld(jj)
          !         alfaOld(jj)=alfaOld(ii)
          !         alfaOld(ii) = temp
          !      endif
          !   enddo
          !enddo
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
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)

    call alfges
 
    DEALLOCATE(xOld,yOld,zOld,alfaOld) 
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
    use ModScbParams,    ONLY: InConAlpha, psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: nisave, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecx, &
                               SORFail, alfa, alphaVal, diffmx, blendAlpha

    use ModScbFunctions, ONLY: extap
  
    use nrtype, ONLY: DP, pi_d
  
    implicit none
 
    INTEGER :: j, i, jz
    REAL(DP) :: rjac, omegaOpt

    INTEGER, ALLOCATABLE :: ni(:)
    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:), om(:), resid(:,:,:)
    
    INTEGER, SAVE :: k, kp, km, iz, im, ip
    !$OMP THREADPRIVATE(k,kp,km,iz,im,ip)
  
    ALLOCATE(alfaprev(nthe,npsi,nzeta+1))
    ALLOCATE(resid(nthe,npsi,nzeta+1))
    ALLOCATE(ni(npsi), om(ny))
    alfaprev = 0._dp; resid = 0._dp; ni = 0.0; om = 0._dp

    rjac = 1._dp - 2._dp*pi_d*pi_d / (REAL(nzeta,dp)*REAL(nzeta,dp) + REAL(nthe,dp)*REAL(nthe,dp))  
    ! Radius of convergence for the Jacobi method, can be used to find optimal SOR omega
  
    om = 1._dp
    !omegaOpt = 2._dp / (1._dp + pi_d/REAL(nzeta+nthe, DP))
    omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))

    alfaPrev(:,:,:) = alfa(:,:,:)
    ni = 0
    resid = 0

!$OMP PARALLEL DO
    psiloop: DO  jz = 2, npsi-1
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
                if (isnan(alfa(iz,jz,k))) then
                   if (verbose) write(*,*) iz,jz,k,alfa(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'NaN encountered in ModScbEuler iterateAlpha'
                   alfa(iz,jz,k) = alfaprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                elseif (alfa(iz,jz,k).ge.1e10) then
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
    sumdb = SUM(ABS(alfa(2:nthe-1,2:npsi-1,2:nzeta) - alfaprev(2:nthe-1,2:npsi-1,2:nzeta)))
    sumb = SUM(ABS(alfa(2:nthe-1,2:npsi-1,2:nzeta)))
    diffmx = maxval(abs(resid(2:nthe-1,2:npsi-1,2:nzeta)))

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
       if (psiChange == 0) then
          DO i = 2,nthe-1
             CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsi-1,k),alfa(i,npsi,k))
             CALL extap(alfa(i,4,k),alfa(i,3,k),alfa(i,2,k),alfa(i,1,k))
          ENDDO
       endif
       if (theChange == 0) then
          DO j = 1,npsi
             CALL extap(alfa(nthe-3,j,k),alfa(nthe-2,j,k),alfa(nthe-1,j,k),alfa(nthe,j,k))
             CALL extap(alfa(4,j,k),alfa(3,j,k),alfa(2,j,k),alfa(1,j,k))
          ENDDO
       endif
    ENDDO

    DEALLOCATE(alfaPrev,resid)
    DEALLOCATE(ni,om)
 
    RETURN
  
  END SUBROUTINE iterateAlpha
  
!================================================!
!============ Psi Euler Potential ===============!
!================================================!
  SUBROUTINE psiFunctions
  
    USE ModScbGrids,     ONLY: npsi
    use ModScbVariables, ONLY: rhoVal, f, fp, psival
 
    use ModRamGSL, ONLY: GSL_Derivs 
  
    implicit none
  
    INTEGER :: GSLerr
  
    CALL GSL_Derivs(rhoVal(1:npsi), psival(1:npsi), f(1:npsi), GSLerr)
    CALL GSL_Derivs(rhoVal(1:npsi), f(1:npsi), fp(1:npsi), GSLerr)
  
    RETURN
  END SUBROUTINE psiFunctions

!==============================================================================
  SUBROUTINE InterpolatePsiR
    ! Interpolates the new values of psi at the new locations xnew on midnight equator
    !!!! Module Variables
    USE ModScbParams,    ONLY: iAzimOffset, psiChange
    USE ModScbGrids,     ONLY: npsi, nzeta
    use ModScbVariables, ONLY: x, y, nThetaEquator, psiVal, kmax, radEqmidNew
    !!!! Module Subroutines/Functions
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP

    implicit none

    INTEGER :: j, k, GSLerr
    REAL(DP) :: deltaR, distConsecFluxSqOld, distConsecFluxSq
    REAL(DP), ALLOCATABLE :: radEqmid(:), psiVal1D(:), radius(:)

    ALLOCATE(radEqMid(npsi), psiVal1D(npsi), radius(npsi))
    radEqMid = 0.0; psiVal1D = 0.0; radius = 0.0

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

    psival1D = psival
    CALL GSL_Interpolation_1D('Cubic',radEqMid, psiVal1D, radEqMidNew(2:npsi), psiVal(2:npsi), GSLerr)

    DEALLOCATE(radEqMid,psiVal1D,radius)
    RETURN

  END SUBROUTINE InterpolatePsiR

!==============================================================================
  SUBROUTINE psiges
    USE ModScbGrids,     ONLY: npsi
    use ModScbVariables, ONLY: psival, psi
  
    INTEGER :: j
  
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
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: x, y, z, psi, psiVal
 
    use ModRamGSL, ONLY: GSL_Interpolation_1D 
  
    implicit none
  
    INTEGER :: k, i, j, GSLerr, i1, i2, ii, jj
    REAL(DP) :: temp
    REAL(DP), ALLOCATABLE :: xOld(:), yOld(:), zOld(:), psiOld(:)

    ALLOCATE(xOld(npsi), yOld(npsi), zOld(npsi), psiOld(npsi))
    xOld = 0.0; yOld = 0.0; zOld = 0.0; psiOld = 0.0

    kloop: DO k = 2, nzeta
       iloop: DO i = 1+theChange, nthe-theChange
          xOld(1:npsi) = x(i,1:npsi,k)
          yOld(1:npsi) = y(i,1:npsi,k)
          zOld(1:npsi) = z(i,1:npsi,k)
          psiOld(1:npsi) = psi(i,1:npsi,k)
          !do ii = 1, npsi
          !   do jj = ii, npsi
          !      if (psiOld(ii) > psiOld(jj)) then
          !         temp = psiOld(jj)
          !         psiOld(jj) = psiOld(ii)
          !         psiOld(ii) = temp
          !      endif
          !   enddo
          !enddo
          !DO j = 2,npsi
          !   if (psiOld(j).lt.psiOld(j-1)) then
          !      psiOld(j) = psiOld(j-1)+1E-6
          !   endif
          !ENDDO

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
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)
  
    call psiges

    DEALLOCATE(xOld,yOld,zOld,psiOld)  
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
    use ModScbParams,    ONLY: InConPsi, psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecr, &
                               SORFail, psi, psiVal, blendPsi, diffmx

    use ModScbFunctions, ONLY: extap
  
    use nrtype, ONLY: pi_d
  
    implicit none
 
    REAL(DP) :: omegaOpt, rjac
    REAL(DP), ALLOCATABLE :: psiPrev(:,:,:), resid(:,:,:), om(:)
    INTEGER, ALLOCATABLE :: ni(:)
    INTEGER :: j, k, i

    INTEGER, SAVE :: jz, jp, jm, iz, im, ip
    !$OMP THREADPRIVATE(jz,jp,jm,iz,im,ip)

    !     perform SOR iteration
    !..   choose for inner iteration loop  
    ALLOCATE(psiPrev(nthe,npsi,nzeta+1))
    ALLOCATE(resid(nthe,npsi,nzeta+1))
    ALLOCATE(ni(nzeta), om(na))
    psiPrev = 0.0; resid = 0.0; ni = 0.0; om = 0.0

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
          jLoop: DO jz = 2, npsi-1
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
                if (isnan(psi(iz,jz,k))) then
                   if (verbose) write(*,*) iz,jz,k,psi(iz,jz,k), resid(iz,jz,k), vecd(iz,jz,k)
                   if (verbose) write(*,*) 'NaN encountered in ModScbEuler iteratePsi'
                   psi(iz,jz,k) = psiprev(iz,jz,k)
                   resid(iz,jz,k) = 0._dp
                   SORFail = .true.
                   EXIT Iterations
                elseif (psi(iz,jz,k).ge.1e10) then
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
          END DO
       END DO
    END DO
    psi(:,:,1) = psi(:,:,nzeta)
    psi(:,:,nzetap) = psi(:,:,2)

    DO k = 2, nzeta
       if (psiChange == 0) then
          DO i = 2,nthe-1
             CALL extap(psi(i,npsi-3,k),psi(i,npsi-2,k),psi(i,npsi-1,k),psi(i,npsi,k))
             CALL extap(psi(i,4,k),psi(i,3,k),psi(i,2,k),psi(i,1,k))
          ENDDO
       endif
       if (theChange == 0) then
          DO j = 1,npsi
             CALL extap(psi(nthe-3,j,k),psi(nthe-2,j,k),psi(nthe-1,j,k),psi(nthe,j,k))
             CALL extap(psi(4,j,k),psi(3,j,k),psi(2,j,k),psi(1,j,k))
          ENDDO
       endif
    ENDDO

    IF (ALLOCATED(psiPrev)) DEALLOCATE(psiPrev)
    IF (ALLOCATED(resid)) DEALLOCATE(resid)
    DEALLOCATE(ni, om)

    RETURN
  
  END SUBROUTINE iteratePsi
  
  END MODULE ModScbEuler
