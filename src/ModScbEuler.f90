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
    !USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: diffmx, rjac, nisave,  x, y, z, sumb, sumdb, chiVal

    !!!! Module Subroutines/Functions
    use ModRamGSL, ONLY: GSL_Interpolation_1D
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d

    IMPLICIT NONE

    INTEGER :: i, j, k, i1, i2, GSLerr
    REAL(DP), DIMENSION(nthe) :: xOld, yOld, zOld, distance, chiValOld
    INTEGER :: psiChange = 0, theChange = 0

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

    IMPLICIT NONE

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

    IMPLICIT NONE

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
 
    INTEGER :: i, j, k
  
    DO k = 1, nzetap
       alfa(1:nthe,1:npsi,k) = alphaval(k)
    END DO

  END SUBROUTINE alfges
                               
!==============================================================================
  SUBROUTINE mapAlpha(iSmoothMove)
    ! new cubic GSL interpolation, without involving linear distance calculation
    USE ModScbMain,      ONLY: DP
    !USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, ny, nthem, nzetap
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, alfaPrev, &
                               left, right
  
    USE ModRamGSL, ONLY: GSL_Interpolation_1D
  
    use nrtype, ONLY: pi_d
  
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: iSmoothMove
    REAL(DP), DIMENSION(nzeta+1) :: xOld, yOld, zOld, alfaOld
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
    INTEGER :: i, j, k, ierr, GSLerr
    integer :: psiChange = 0, theChange = 0

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
!    USE ModScbMain,      ONLY: DP
!    use ModScbGrids,     ONLY: nzeta, nzetap, nthe, nthem, npsi, npsim
!    use ModScbVariables, ONLY: sumb, sumdb, vecd, vec1, vec2, &
!                               vec3, vec4, vec6, vec7, vec8, vec9, vecx, &
!                               left, right
!  
!    use nrtype, ONLY: twopi_d, pi_d
!    use nrmod, ONLY: gaussj
! 
!    IMPLICIT NONE
!  
!    ! define the (NSIZE-1) vectors Fk and [(NSIZE-1)x(NSIZE-1)] matrices
!    ! Ek, & 
!    ! which will be kept till end
!  
!    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:)
!    REAL(DP), DIMENSION(:, :), ALLOCATABLE :: FMATRIXFULL, UVECTORFULL
!    REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: EMATRIXFULL
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: AMATRIX, BMATRIX, CMATRIX, EMATRIX, &
!         MATRIX1, MATVER, UNITMAT, RESID, MATRIXINV
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: DMATRIX, FMATRIX, MATRIX2, UVECTOR, &
!         UVECTORPLUS
!    REAL(DP), ALLOCATABLE :: AUXMATR(:,:)            !!! :: AUXMATR(nzeta-1,1)
!    REAL(DP) :: tpi2, tpi3
!  
!    INTEGER :: jz, ierrffull, ierrefull, ierra, ierrb, k, ierrc, ierrd, ierre, ierrf, &
!         ierrmat1, ierrmatver, ierrmatrixinversion, ierresid, ierrunit, i, j, ierrdamat, ierrdbmat, &
!         ierrdcmat, ierrddmat, ierrdemat, ierrdfmat, ierrdm1mat, ierrmdmatver, ierrdmatrixinversion, &
!         ierrdunitmat, ierrdm2mat, ierruv, ierruvp, ierr, izmx, jzmx, kmx, ierrdmatver, ierraux, &
!         ierrdaux, idealerr
!  
!    INTEGER :: maximum_location(3)
!  
!    ! Start here the loop for different flux surfaces
!  
!    ALLOCATE(alfaPrev(SIZE(alfa,1), SIZE(alfa,2), SIZE(alfa,3)), stat = ierr)
!  
!    alfaprev(:,:,:) = alfa(:,:,:)
!  
!    sumdb = 0._dp
!    sumb = 0._dp
!    diffmx = 0._dp
!  
!    fluxloop: DO jz = left+1, right-1   ! jz represents the flux surface
!  
!       ALLOCATE(FMATRIXFULL(nthe, nzeta-1), STAT=ierrffull)
!       ALLOCATE(EMATRIXFULL(nthe, nzeta-1, nzeta-1), STAT=ierrefull)
!       ALLOCATE(UVECTORFULL(nthe, nzeta-1), STAT=ierrffull)
!  
!       EMATRIXFULL(1,:,:) = 0._dp
!  
!       UVECTORFULL(1,1:nzeta-1) = alphaVal(1:nzeta-1)     ! BC on theta ; here first index 
!       ! is theta
!       UVECTORFULL(nthe,1:nzeta-1) = alphaVal(1:nzeta-1) ! BC on theta
!  
!       FMATRIXFULL(1,:) = UVECTORFULL(1,:)
!  
!  
!       ! for k from 2 to NSIZE-1, define Ak, Bk, Ck, Dk
!       ! after each "step" k, these can be de-allocated, and only Ekfull, Fkfull
!       ! kept
!       ! simple ek, fk will also be de-allocated
!  
!       k = 2  ! this is theta index
!  
!       kloop: DO
!  
!          ALLOCATE(AMATRIX(nzeta-1,nzeta-1), STAT=ierra)
!          ALLOCATE(BMATRIX(nzeta-1,nzeta-1), STAT=ierrb)
!          ALLOCATE(CMATRIX(nzeta-1,nzeta-1), STAT=ierrc)
!          ALLOCATE(DMATRIX(nzeta-1), STAT=ierrd)
!          ALLOCATE(EMATRIX(nzeta-1,nzeta-1), STAT=ierre)
!          ALLOCATE(FMATRIX(nzeta-1), STAT=ierrf)
!          ALLOCATE(MATRIX1(nzeta-1,nzeta-1), STAT=ierrmat1)
!          ALLOCATE(MATRIX2(nzeta-1), STAT=ierrmat1)
!          ALLOCATE(MATVER(nzeta-1,nzeta-1), STAT=ierrmatver)
!          ALLOCATE(MATRIXINV(nzeta-1,nzeta-1), STAT=ierrmatrixInversion)
!          ALLOCATE(RESID(nzeta-1,nzeta-1), STAT=ierresid)
!          ALLOCATE(UNITMAT(nzeta-1,nzeta-1), STAT=ierrunit)
!          ALLOCATE(AUXMATR(nzeta-1, 1), STAT=ierraux)
!  
!          AMATRIX = 0.
!          BMATRIX = 0.
!          CMATRIX = 0.
!  
!          IF ((ierra /= 0) .OR. (ierrb /=0) .OR. (ierrc /= 0) .OR. (ierrd /= 0)) &
!               STOP "A,B,C,D : Allocation failed"
!  
!          DO i = 1, nzeta-3
!             DMATRIX(i) = vecx(k, jz, i+1)
!          END DO
!  
!          tpi2 = twopi_d ! yy(k,jz,nzeta) - yy(k,jz,1)
!          tpi3 = twopi_d ! yy(k,jz,nzetap) - yy(k,jz,2)
!  
!          DMATRIX(nzeta-2) = vecx(k, jz, nzeta-1) -tpi2*vec7(k, jz, nzeta-1) &
!               -tpi2*vec8(k,jz,nzeta-1)-tpi2*vec9(k,jz,nzeta-1)
!          DMATRIX(nzeta-1) = vecx(k, jz, nzeta) + tpi2*vecd(k, jz, nzeta) &
!               - tpi2*vec4(k, jz, nzeta) - &
!               tpi2*vec6(k, jz, nzeta) -tpi3*vec7(k, jz, nzeta) &
!               -tpi3*vec8(k, jz, nzeta) &
!               -tpi3*vec9(k, jz, nzeta)
!  
!          FMATRIX = FMATRIXFULL(k-1,:)
!          EMATRIX = EMATRIXFULL(k-1,:,:)
!  
!          jloop: DO j = 1, nzeta-1
!             iloop: DO i = 1, nzeta-1
!                IF (i == j)  THEN
!                   AMATRIX(i,j) = vec3(k, jz, i+1)
!                   BMATRIX(i,j) = -vec2(k, jz, i+1)
!                   CMATRIX(i,j) = vec1(k, jz, i+1)
!                END IF
!                IF (j == i+1) THEN
!                   AMATRIX(i,j) = vec6(k, jz, i+1)
!                   BMATRIX(i,j) = vecd(k, jz, i+1)
!                   CMATRIX(i,j) = vec4(k, jz, i+1)
!                END IF
!                IF (j == i+2) THEN
!                   AMATRIX(i,j) = vec9(k, jz, i+1)
!                   BMATRIX(i,j) = -vec8(k, jz, i+1)
!                   CMATRIX(i,j) = vec7(k, jz, i+1)
!                END IF
!                IF (i == nzeta-2) THEN
!                   AMATRIX(i,1) = vec9(k, jz, nzeta-1)
!                   BMATRIX(i,1) = -vec8(k, jz, nzeta-1)
!                   CMATRIX(i,1) = vec7(k, jz, nzeta-1)
!                END IF
!                IF(i == nzeta-1) THEN
!                   AMATRIX(i,1) = vec6(k, jz, nzeta)
!                   AMATRIX(i,2) = vec9(k, jz, nzeta)
!                   BMATRIX(i,1) = vecd(k, jz, nzeta)
!                   BMATRIX(i,2) = -vec8(k, jz, nzeta)
!                   CMATRIX(i,1) = vec4(k, jz, nzeta)
!                   CMATRIX(i,2) = vec7(k, jz, nzeta)
!                END IF
!             END DO iloop
!          END DO jloop
!  
!          ! AMATRIX, BMATRIX, CMATRIX, DMATRIX, EMATRIXk, FMATRIXk  &
!          ! are initialized at this point
!  
!          !Time to do the inverse of (Bk - Ck*Ek-1)
!  
!          MATRIX1 = BMATRIX - MATMUL(CMATRIX, EMATRIX)
!          MATRIX2 = MATMUL(CMATRIX,FMATRIX) - DMATRIX
!  
!          AUXMATR = 100._dp
!  
!          MATVER = MATRIX1
!          UNITMAT = 0._dp
!  
!          DO i = 1, nzeta-1
!             UNITMAT(i,i) = 1._dp
!          END DO
!  
!          ! Now time to invert MATRIX1
!  
!          ! CALL gaussj(MATRIX1, nzeta-1, nzeta-1, AUXMATR, 1, 1, nmaximum) 
!          CALL gaussj(MATRIX1, AUXMATR)  ! f90 routine (NR) with full pivoting
!  
!          ! CALL matrixInversion(nzeta-1, 0, MATRIX1, auxmatr, det)  ! French
!          ! routine (fast)
!  
!          ! Now MATRIX1 is the inverse, do iterative improvement of it 
!  
!          RESID = MATMUL(MATVER, MATRIX1)
!          MATRIXINV = MATRIX1 + MATMUL(MATRIX1, UNITMAT-RESID)
!  
!          ! MATRIXINV is the iterative improvement of inv. of MATRIX1 
!  
!          EMATRIX = MATMUL(MATRIXINV, AMATRIX)
!          FMATRIX = MATMUL(MATRIXINV, MATRIX2)
!  
!          ! Now we define the full E, F for this order (k)
!  
!          FMATRIXFULL(k,:) = FMATRIX
!          EMATRIXFULL(k,:,:) = EMATRIX
!  
!          ! Now we have what we need, E and F for this order (k), 
!          ! So we can de-allocate the matrices we defined here.
!  
!          IF (ALLOCATED(AMATRIX)) DEALLOCATE(AMATRIX,STAT=ierrdamat)
!          IF (ALLOCATED(BMATRIX)) DEALLOCATE(BMATRIX,STAT=ierrdbmat)
!          IF (ALLOCATED(CMATRIX)) DEALLOCATE(CMATRIX,STAT=ierrdcmat)
!          IF (ALLOCATED(DMATRIX)) DEALLOCATE(DMATRIX,STAT=ierrddmat)
!          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierrdemat)
!          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierrdfmat)
!          IF (ALLOCATED(MATRIX1)) DEALLOCATE(MATRIX1,STAT=ierrdm1mat)
!          IF (ALLOCATED(MATVER)) DEALLOCATE(MATVER,STAT=ierrdmatver)
!          IF (ALLOCATED(MATRIXINV)) DEALLOCATE(MATRIXINV,STAT=ierrdmatrixInversion)
!          IF (ALLOCATED(RESID)) DEALLOCATE(RESID,STAT=ierresid)
!          IF (ALLOCATED(UNITMAT)) DEALLOCATE(UNITMAT,STAT=ierrdunitmat)
!          IF (ALLOCATED(MATRIX2)) DEALLOCATE(MATRIX2,STAT=ierrdm2mat)
!          IF (ALLOCATED(AUXMATR)) DEALLOCATE(AUXMATR, STAT=ierrdaux)
!  
!          IF (k == nthe-1) EXIT kloop   ! We only need E, F up to the order N-1
!  
!          k = k + 1
!  
!       END DO kloop
!  
!       ! Now it's time to start calculating the U vectors, keeping in mind that 
!       ! U(N) is known from the BCs
!  
!       k = nthe - 1
!  
!       kloop2: DO
!  
!          ALLOCATE(UVECTOR(nzeta-1), STAT=ierruv)
!          ALLOCATE(UVECTORPLUS(nzeta-1), STAT=ierruvp)
!          ALLOCATE(EMATRIX(nzeta-1, nzeta-1), STAT = ierre)
!          ALLOCATE(FMATRIX(nzeta-1), STAT = ierrf)
!  
!          FMATRIX = FMATRIXFULL(k,:)
!          UVECTORPLUS = UVECTORFULL(k+1,:)
!          EMATRIX = EMATRIXFULL(k,:,:)
!  
!          UVECTOR = MATMUL(EMATRIX, UVECTORPLUS) + FMATRIX
!  
!          UVECTORFULL(k,:) = UVECTOR
!  
!          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierr)
!          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierr)
!          IF (ALLOCATED(UVECTOR)) DEALLOCATE(UVECTOR,STAT=ierr)
!          IF (ALLOCATED(UVECTORPLUS)) DEALLOCATE(UVECTORPLUS,STAT=ierr)
!  
!          IF (k == 2) EXIT kloop2
!  
!          k = k - 1
!  
!       END DO kloop2
!  
!       ! Now all data is in the matrix UVECTORFULL
!  
!       ! De-allocate Efull, Ffull
!  
!       IF (ALLOCATED(EMATRIXFULL)) DEALLOCATE(EMATRIXFULL,STAT=ierr)
!       IF (ALLOCATED(FMATRIXFULL)) DEALLOCATE(FMATRIXFULL,STAT=ierr)
!  
!       alfa(1:nthe,jz,1:nzeta-1) = UVECTORFULL(1:nthe,1:nzeta-1)
!  
!       DO k = 1, nthe
!          tpi2 = twopi_d
!          tpi3 = twopi_d
!  
!          alfa(k, jz, nzeta) = alfa(k, jz, 1) + tpi2
!          alfa(k, jz, nzeta+1) = alfa(k, jz, 2) + tpi3
!       END DO
!  
!       IF (ALLOCATED(UVECTORFULL)) DEALLOCATE(UVECTORFULL,STAT=ierr)
!  
!    END DO fluxloop
!  
!    DO k = 2,nzeta
!       DO i = 2,nthem
!          CALL extap(alfa(i,4,k),alfa(i,3,k),alfa(i,2,k),alfa(i,1,k))
!          CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsim,k), &
!                     alfa(i,npsi,k))
!       END DO
!    END DO
!  
!    DO  j = left,right
!       DO  i = 1,nthe
!          DO k = 2, nzeta
!             alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1. - blendAlpha)
!          END DO
!  
!          !  periodic boundary condition in delta, alfa = phi + delta
!          alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
!          alfa(i,j,nzetap) = alfa(i,j,2) + 2._dp * pi_d
!       END DO
!    END DO
!  
!    sumdb = SUM(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
!    sumb = SUM(ABS(alfa(2:nthem,1:npsi,2:nzeta)))
!    diffmx = MAXVAL(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
!  
!    maximum_location =  MAXLOC(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
!  
!    izmx = maximum_location(1)
!    jzmx = maximum_location(2)
!    kmx = maximum_location(3)
!  
!    ierr = 2
!  
!    DEALLOCATE(alfaPrev, stat = idealerr)
!  
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
  
    IMPLICIT NONE
    save
 
    REAL(DP) :: anorm, diff, rjac, &
         & anormaverage, dyDummy, anormResid, anormError, anormf
    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:), alfaPrevTemp(:,:,:), om(:)
    REAL(DP) :: RESULT, omegaOpt
    REAL(DP), ALLOCATABLE :: angle(:,:,:), resid(:,:,:)
    INTEGER, ALLOCATABLE :: ni(:)
    INTEGER :: j, i, ict, jz, k, kp, km, iz, im, ip, izmx, jzmx, kmx, ierr, nratio, &
         myPsiBegin, myPsiEnd, my_array_type2, psiRangeDiff, resultInt, loc(3)
    !$OMP THREADPRIVATE(k,kp,km,iz,im,ip)
  
    ALLOCATE(angle(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaprev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(ni(npsi), om(ny))

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
  
    IMPLICIT NONE
  
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

    IMPLICIT NONE

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
    !USE ModScbParams,    ONLY: psiChange, theChange
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, left, right
 
    use ModRamGSL, ONLY: GSL_Interpolation_1D 
  
    IMPLICIT NONE
  
    INTEGER :: iSmoothMove = 0
    INTEGER :: ierr, k, i, j, GSLerr, i1, i2
    REAL(DP), DIMENSION(npsi) :: xOld, yOld, zOld, psiOld
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
    integer :: psiChange = 0, theChange = 0

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
!    USE ModScbMain,      ONLY: DP
!    USE ModScbGrids,     ONLY: nzeta, nzetap, nthe, nthem, npsi, npsim
!    use ModScbVariables, ONLY: sumb, sumdb, vecd, vec1, vec2, &
!                               vec3, vec4, vec6, vec7, vec8, vec9, vecr, &
!                               left, right
! 
!    use nrmod, ONLY: gaussj 
!    IMPLICIT NONE
!  
!    REAL(DP), ALLOCATABLE :: psiprev(:,:,:)
!    INTEGER :: ierra, ierrb, ierrc, ierrd, i, j, k
!    REAL(DP), DIMENSION(:, :), ALLOCATABLE :: FMATRIXFULL, UVECTORFULL
!    REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: EMATRIXFULL
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: AMATRIX, BMATRIX, CMATRIX, EMATRIX, &
!         MATRIX1, MATVER, MATRIXINV, RESID
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: DMATRIX, FMATRIX, MATRIX2, UVECTOR, &
!         UVECTORPLUS
!  
!    REAL(DP), ALLOCATABLE :: AUXMATR(:,:)  !!! DIMENSION(nthe-2,1) :: AUXMATR
!    REAL(DP), ALLOCATABLE :: UNITMAT(:,:)  !!! :: UNITMAT(nthe-2,nthe-2)
!    REAL(DP) :: diff
!  
!    INTEGER, DIMENSION(3) :: maximum_location
!  
!    INTEGER :: jz, ierrffull, ierrefull,  ierre, ierrf, &
!         ierrmat1, ierrmatver, ierrmatrixinversion, ierresid, ierrdamat, ierrdbmat, &
!         ierrdcmat, ierrddmat, ierrdemat, ierrdfmat, ierrdm1mat, ierrmdmatver, ierrdmatrixinversion, &
!         ierrdunitmat, ierrdm2mat, ierruv, ierruvp, ierr, izmx, jzmx, kmx, ierrdmatver, ierrmatv, kzeta, &
!         ierrmatres, ierrdmver, ierrdresid, iz, ierraux, ierrunit, ierrdaux, ierrdunit, idealerr
!  
!    ALLOCATE(psiPrev(SIZE(psi,1), SIZE(psi,2), SIZE(psi,3)), stat = ierr)
!  
!    psiprev(:,:,:) = psi(:,:,:)
!  
!    sumdb = 0._dp
!    sumb = 0._dp
!    diffmx = 0._dp
!  
!    ! The subroutine for the psi equation is a bit different from that for alpha;
!  
!    ! Start here the loop for different alpha surfaces
!  
!    zetaloop: DO kzeta = 2, nzeta
!  
!       ! PRINT*, 'Alpha surface:', kzeta
!  
!       ALLOCATE(FMATRIXFULL(npsi, nthe-2), STAT=ierrffull)
!       ALLOCATE(EMATRIXFULL(npsi, nthe-2, nthe-2), STAT=ierrefull)
!       ALLOCATE(UVECTORFULL(npsi, nthe), STAT=ierrffull)
!  
!       EMATRIXFULL(1,:,:) = 0.   ! first index is the flux surface
!  
!       UVECTORFULL(1,:) = psival(1)
!       UVECTORFULL(npsi,:) = psival(npsi)
!  
!       FMATRIXFULL(1,:) = psival(1)
!  
!       k = 2    ! k are the flux surfaces here
!  
!       kloop: DO 
!  
!          ALLOCATE(AMATRIX(nthe-2,nthe-2), STAT=ierra)
!          ALLOCATE(BMATRIX(nthe-2,nthe-2), STAT=ierrb)
!          ALLOCATE(CMATRIX(nthe-2,nthe-2), STAT=ierrc)
!          ALLOCATE(DMATRIX(nthe-2), STAT=ierrd)
!          ALLOCATE(EMATRIX(nthe-2,nthe-2), STAT=ierre)
!          ALLOCATE(FMATRIX(nthe-2), STAT=ierrf)
!          ALLOCATE(MATRIX1(nthe-2,nthe-2), STAT=ierrmat1)
!          ALLOCATE(MATVER(nthe-2,nthe-2), STAT=ierrmatv)
!          ALLOCATE(MATRIXINV(nthe-2,nthe-2), STAT=ierrmatrixInversion)
!          ALLOCATE(RESID(nthe-2,nthe-2), STAT=ierrmatres)
!          ALLOCATE(MATRIX2(nthe-2), STAT=ierrmat1)
!          ALLOCATE(AUXMATR(nthe-2, 1), STAT=ierraux)
!          ALLOCATE(UNITMAT(nthe-2, nthe-2), STAT=ierrunit)
!  
!          AMATRIX = 0.
!          BMATRIX = 0.
!          CMATRIX = 0.
!  
!          IF ((ierra /= 0) .OR. (ierrb /=0) .OR. (ierrc /= 0) .OR. (ierrd /= 0)) &
!               STOP "A,B,C,D : Allocation failed"
!  
!          DO i = 1, nthe-2
!             DMATRIX(i) = vecr(i+1, k, kzeta)
!  
!             IF (i==1) THEN
!                DMATRIX(i) = vecr(i+1, k, kzeta) - &
!                     vec1(i+1, k, kzeta)*psi(i,k-1,kzeta) &
!                     -vec4(i+1, k, kzeta)*psi(i,k,kzeta) &
!                     -vec7(i+1, k, kzeta)*psi(i,k+1,kzeta)
!             END IF
!  
!             IF (i==nthe-2) THEN
!                DMATRIX(i) = vecr(i+1, k, kzeta) - &
!                     vec3(i+1, k, kzeta)*psi(i,k-1,kzeta) &
!                     -vec6(i+1, k, kzeta)*psi(i,k,kzeta) &
!                     -vec9(i+1, k, kzeta)*psi(i,k+1,kzeta)
!             END IF
!  
!          END DO
!  
!          EMATRIX(:,:) = EMATRIXFULL(k-1,:,:)
!          FMATRIX(:) = FMATRIXFULL(k-1,:) ! FMATRIX of order k-1 
!  
!          jloop: DO j = 1, nthe-2
!  
!             iloop: DO i = 1, nthe-2
!  
!                IF (i == j)  THEN
!                   AMATRIX(i,j) = vec8(i+1, k, kzeta)
!                   BMATRIX(i,j) = vecd(i+1, k, kzeta)
!                   CMATRIX(i,j) = vec2(i+1, k, kzeta)
!                END IF
!  
!                IF (j == i+1) THEN
!                   AMATRIX(i,j) = vec9(i+1, k, kzeta)
!                   BMATRIX(i,j) = -vec6(i+1, k, kzeta)
!                   CMATRIX(i,j) = vec3(i+1, k, kzeta)
!                END IF
!  
!                IF (j == i-1) THEN
!                   AMATRIX(i,j) = vec7(i+1, k, kzeta)
!                   BMATRIX(i,j) = -vec4(i+1, k, kzeta)
!                   CMATRIX(i,j) = vec1(i+1, k, kzeta)
!                END IF
!  
!             END DO iloop
!          END DO jloop
!  
!  
!          ! AMATRIX, BMATRIX, CMATRIX, DMATRIX, EMATRIXk, FMATRIXk  are initialized &
!          ! at this point
!  
!          !Time to do the inverse of (Bk - Ck*Ek-1)
!  
!          MATRIX1 = BMATRIX - MATMUL(CMATRIX, EMATRIX)
!          MATRIX2 = MATMUL(CMATRIX,FMATRIX) - DMATRIX
!  
!          ! Now time to invert MATRIX1
!  
!          AUXMATR = 100.   ! just a test value, only the inverse of matrix1 interests us
!  
!          UNITMAT = 0.
!  
!          DO i = 1,nthe-2
!             UNITMAT(i,i) = 1.
!          END DO
!  
!          MATVER = MATRIX1  ! MATVER keeps matrix1
!  
!          ! CALL gaussj(MATRIX1, nthe-2, nthe-2, AUXMATR, 1, 1, nmaximum)    
!  
!          CALL gaussj(MATRIX1,AUXMATR)   ! f90 routine (NR) with pivoting
!  
!          ! CALL matrixInversion(nthe-2, 0, MATRIX1, auxmatr, det)  ! French routine
!  
!          ! Now on output MATRIX1 should be the inverse
!  
!          ! Start iterative improvement for the solution
!  
!          RESID = MATMUL(MATVER, MATRIX1)
!  
!          MATRIXINV = MATRIX1 + MATMUL(MATRIX1, UNITMAT-RESID)
!  
!          ! MATRIXINV is the corrected value for the inverse
!  
!          EMATRIX = MATMUL(MATRIXINV, AMATRIX)
!          FMATRIX = MATMUL(MATRIXINV, MATRIX2)
!  
!          ! Now we define the full E, F for this order (k)
!  
!          FMATRIXFULL(k, :) = FMATRIX(:)
!          EMATRIXFULL(k, :, :) = EMATRIX(:, :)
!  
!          ! Now we have what we need, E and F for this order (k), 
!          ! So we can de-allocate the matrices we defined here.
!  
!          IF (ALLOCATED(AMATRIX)) DEALLOCATE(AMATRIX,STAT=ierrdamat)
!          IF (ALLOCATED(BMATRIX)) DEALLOCATE(BMATRIX,STAT=ierrdbmat)
!          IF (ALLOCATED(CMATRIX)) DEALLOCATE(CMATRIX,STAT=ierrdcmat)
!          IF (ALLOCATED(DMATRIX)) DEALLOCATE(DMATRIX,STAT=ierrddmat)
!          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierrdemat)
!          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierrdfmat)
!          IF (ALLOCATED(MATRIX1)) DEALLOCATE(MATRIX1,STAT=ierrdm1mat)
!          IF (ALLOCATED(MATVER)) DEALLOCATE(MATVER,STAT=ierrdmver)
!          IF (ALLOCATED(MATRIXINV)) DEALLOCATE(MATRIXINV,STAT=ierrmatrixInversion)
!          IF (ALLOCATED(RESID)) DEALLOCATE(MATVER,STAT=ierrdresid)
!          IF (ALLOCATED(MATRIX2)) DEALLOCATE(MATRIX2,STAT=ierrdm2mat)
!          IF (ALLOCATED(AUXMATR)) DEALLOCATE(AUXMATR, STAT=ierrdaux)
!          IF (ALLOCATED(UNITMAT)) DEALLOCATE(UNITMAT, STAT=ierrdunit)
!  
!          IF (k == npsi - 1) EXIT kloop   ! We only need E, F up to the order N-1
!  
!          k = k + 1
!  
!       END DO kloop
!  
!       ! Now it's time to start calculating the U vectors, keeping in mind that 
!       ! U(N) is known from the BCs
!  
!       k = npsi - 1
!  
!       kloop2: DO
!  
!          ALLOCATE(UVECTOR(nthe-2), STAT=ierruv)
!          ALLOCATE(UVECTORPLUS(nthe-2), STAT=ierruvp)
!          ALLOCATE(EMATRIX(nthe-2, nthe-2), STAT = ierre)
!          ALLOCATE(FMATRIX(nthe-2), STAT = ierrf)
!  
!          FMATRIX(:) = FMATRIXFULL(k,:)
!          EMATRIX(:,:) = EMATRIXFULL(k,:,:)
!  
!          DO i = 1, nthe-2
!             UVECTORPLUS(i) = UVECTORFULL(k+1,i+1)
!          END DO
!  
!          UVECTOR = MATMUL(EMATRIX, UVECTORPLUS) + FMATRIX
!  
!          DO i = 1, nthe-2
!             UVECTORFULL(k, i+1) = UVECTOR(i)
!          END DO
!  
!          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierr)
!          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierr)
!          IF (ALLOCATED(UVECTOR)) DEALLOCATE(UVECTOR,STAT=ierr)
!          IF (ALLOCATED(UVECTORPLUS)) DEALLOCATE(UVECTORPLUS,STAT=ierr)
!  
!          IF (k == 2) EXIT kloop2
!  
!          k = k-1
!  
!       END DO kloop2
!  
!       ! Now all data is in the matrix UVECTORFULL
!       ! De-allocate efull, ffull
!  
!       IF (ALLOCATED(EMATRIXFULL)) DEALLOCATE(EMATRIXFULL,STAT=ierr)
!       IF (ALLOCATED(FMATRIXFULL)) DEALLOCATE(FMATRIXFULL,STAT=ierr)
!  
!       DO k = 1, npsi
!          DO i = 2, nthe-1
!             psi(i, k, kzeta) = UVECTORFULL(k,i)
!          END DO
!       END DO
!  
!       ! the values for theta = 1 and NSIZE are the same as defined in psiges
!       ! Now psi is completely determined for zeta surface kzeta = kzeta
!  
!       IF (ALLOCATED(UVECTORFULL)) DEALLOCATE(UVECTORFULL,STAT=ierr)
!  
!    END DO zetaloop
!  
!    !...  set blend for outer iteration loop
!  
!    DO  j = 1,npsi
!       DO  i = 1,nthe
!          DO k = 2, nzeta
!             psi(i,j,k) = psi(i,j,k) * blendPsi + psival(j) * (1. - blendPsi)  
!             !  blend does the "blending" between computed psi at (n+1) and at (n)
!          END DO
!          psi(i,j,1) = psi(i,j,nzeta)
!          psi(i,j,nzetap) = psi(i,j,2)
!       END DO
!    END DO
!  
!    DO k = 2, nzeta
!       DO jz = 2, npsim
!          DO iz = 2, nthem
!             diff = ABS(psi(iz,jz,k)-psiprev(iz,jz,k))
!             IF(diff.GT.diffmx) diffmx = diff
!             sumdb = sumdb + ABS(psi(iz,jz,k)-psiprev(iz,jz,k))
!             sumb = sumb + ABS(psi(iz,jz,k))
!          END DO
!       END DO
!    END DO
!  
!    DEALLOCATE(psiPrev, stat = idealerr)
!  
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
  
    IMPLICIT NONE
    save
 
    REAL(DP) :: omegaOpt
    REAL(DP) :: omc, anorm, anormf, diff, ano, sumbtest, anormResid, anormError
    REAL(DP), ALLOCATABLE :: psiPrev(:,:,:), psiPrevTemp(:,:,:), resid(:,:,:), om(:)
    INTEGER, ALLOCATABLE :: ni(:)
    INTEGER :: j, k, i, ierr, ict, jz, jp, jm, iz, im, ip, nratio, myAlphaBegin, &
               myAlphaEnd, my_array_type_psi2, alphaRangeDiff, resultInt, loc(3)
    !$OMP THREADPRIVATE(jz,jp,jm,iz,im,ip)

    !     perform SOR iteration
    !..   choose for inner iteration loop
  
    ALLOCATE(psiPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(psiPrev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(ni(nzeta), om(na))

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
