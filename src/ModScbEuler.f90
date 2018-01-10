MODULE ModScbEuler
  ! Contains subroutines responsible for calculating the alpha (actually beta)
  ! and psi (actually alpha) Euler potentials
  
  use ModScbVariables, ONLY: psi, psiSav1, psiSav2, alfa, alfaSav1, alfaSav2, &
                             alphaVal, alphaValInitial, psiVal, blendAlpha, &
                             blendPsi, diffmx, rjac, decreaseConvAlpha, &
                             decreaseConvPsi, errorAlpha, errorAlphaPrev, &
                             errorPsi, errorPsiPrev, iAlphaMove, iPsiMove
  
  use ModScbFunctions, ONLY: extap
  
  implicit none
  
  contains
  
!================================================!
!=========== Alpha Euler Potential ==============!
!================================================!
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
    ! new cubic spline interpolation, without involving linear distance
    ! calculation
  
    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, ny, nthem, nzetap
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, alfaPrev, &
                               left, right
  
    USE ModScbSpline, ONLY: spline, splint
  
    use nrtype, ONLY: pi_d
  
    IMPLICIT NONE
  
    !  ignore i = 1 & nthe lines since the boundary condition is delta = 0 there
    !  (on Earth surface)
  
    INTEGER, INTENT(IN) :: iSmoothMove
    REAL(DP), DIMENSION(nzeta+1) :: xOld, yOld, zOld, alfaOld, x2derivs, y2derivs, &
                                    z2derivs
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
    INTEGER :: i, j, k, ierr
  
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
  
       jloop : DO j = 1, npsi
          iloop: DO i = 2, nthem
             xOld(1:nzetap) = x(i,j,1:nzetap)
             yOld(1:nzetap) = y(i,j,1:nzetap)
             zOld(1:nzetap) = z(i,j,1:nzetap)
             alfaOld(1:nzetap) = alfa(i,j,1:nzetap)

             CALL spline(alfaOld, xOld, 1.E31_dp, 1.E31_dp, x2derivs)
             CALL spline(alfaOld, yOld, 1.E31_dp, 1.E31_dp, y2derivs)
             CALL spline(alfaOld, zOld, 1.E31_dp, 1.E31_dp, z2derivs)
  
             DO k = 2,nzeta
                x(i,j,k) = splint(alfaOld, xOld, x2derivs, alphaval(k))
                y(i,j,k) = splint(alfaOld, yOld, y2derivs, alphaval(k))
                z(i,j,k) = splint(alfaOld, zOld, z2derivs, alphaval(k))
                alfa(i,j,k) = alphaval(k)
             END DO
  
             ! Periodic boundary conditions
             x(i,j,1) = x(i,j,nzeta)
             y(i,j,1) = y(i,j,nzeta)
             z(i,j,1) = z(i,j,nzeta)
             alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp*pi_d
             x(i,j,nzetap) = x(i,j,2)
             y(i,j,nzetap) = y(i,j,2)
             z(i,j,nzetap) = z(i,j,2)
             alfa(i,j,nzetap) = alfa(i,j,2) + 2._dp*pi_d
  
          END DO iloop
       END DO jloop
    END IF
  
    RETURN
  
  END SUBROUTINE mapAlpha
  
  !------------------------------------------------------------------------------
  SUBROUTINE directAlpha
  
    USE ModScbMain,      ONLY: DP
    use ModScbGrids,     ONLY: nzeta, nzetap, nthe, nthem, npsi, npsim
    use ModScbVariables, ONLY: sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecx, &
                               left, right
  
    use nrtype, ONLY: twopi_d, pi_d
    use nrmod, ONLY: gaussj
 
    IMPLICIT NONE
  
    ! define the (NSIZE-1) vectors Fk and [(NSIZE-1)x(NSIZE-1)] matrices
    ! Ek, & 
    ! which will be kept till end
  
    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:)
    REAL(DP), DIMENSION(:, :), ALLOCATABLE :: FMATRIXFULL, UVECTORFULL
    REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: EMATRIXFULL
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: AMATRIX, BMATRIX, CMATRIX, EMATRIX, &
         MATRIX1, MATVER, UNITMAT, RESID, MATRIXINV
    REAL(DP), DIMENSION(:), ALLOCATABLE :: DMATRIX, FMATRIX, MATRIX2, UVECTOR, &
         UVECTORPLUS
    REAL(DP), ALLOCATABLE :: AUXMATR(:,:)            !!! :: AUXMATR(nzeta-1,1)
    REAL(DP) :: tpi2, tpi3
  
    INTEGER :: jz, ierrffull, ierrefull, ierra, ierrb, k, ierrc, ierrd, ierre, ierrf, &
         ierrmat1, ierrmatver, ierrmatrixinversion, ierresid, ierrunit, i, j, ierrdamat, ierrdbmat, &
         ierrdcmat, ierrddmat, ierrdemat, ierrdfmat, ierrdm1mat, ierrmdmatver, ierrdmatrixinversion, &
         ierrdunitmat, ierrdm2mat, ierruv, ierruvp, ierr, izmx, jzmx, kmx, ierrdmatver, ierraux, &
         ierrdaux, idealerr
  
    INTEGER :: maximum_location(3)
  
    ! Start here the loop for different flux surfaces
  
    ALLOCATE(alfaPrev(SIZE(alfa,1), SIZE(alfa,2), SIZE(alfa,3)), stat = ierr)
  
    alfaprev(:,:,:) = alfa(:,:,:)
  
    sumdb = 0._dp
    sumb = 0._dp
    diffmx = 0._dp
  
    fluxloop: DO jz = left+1, right-1   ! jz represents the flux surface
  
       ALLOCATE(FMATRIXFULL(nthe, nzeta-1), STAT=ierrffull)
       ALLOCATE(EMATRIXFULL(nthe, nzeta-1, nzeta-1), STAT=ierrefull)
       ALLOCATE(UVECTORFULL(nthe, nzeta-1), STAT=ierrffull)
  
       EMATRIXFULL(1,:,:) = 0._dp
  
       UVECTORFULL(1,1:nzeta-1) = alphaVal(1:nzeta-1)     ! BC on theta ; here first index 
       ! is theta
       UVECTORFULL(nthe,1:nzeta-1) = alphaVal(1:nzeta-1) ! BC on theta
  
       FMATRIXFULL(1,:) = UVECTORFULL(1,:)
  
  
       ! for k from 2 to NSIZE-1, define Ak, Bk, Ck, Dk
       ! after each "step" k, these can be de-allocated, and only Ekfull, Fkfull
       ! kept
       ! simple ek, fk will also be de-allocated
  
       k = 2  ! this is theta index
  
       kloop: DO
  
          ALLOCATE(AMATRIX(nzeta-1,nzeta-1), STAT=ierra)
          ALLOCATE(BMATRIX(nzeta-1,nzeta-1), STAT=ierrb)
          ALLOCATE(CMATRIX(nzeta-1,nzeta-1), STAT=ierrc)
          ALLOCATE(DMATRIX(nzeta-1), STAT=ierrd)
          ALLOCATE(EMATRIX(nzeta-1,nzeta-1), STAT=ierre)
          ALLOCATE(FMATRIX(nzeta-1), STAT=ierrf)
          ALLOCATE(MATRIX1(nzeta-1,nzeta-1), STAT=ierrmat1)
          ALLOCATE(MATRIX2(nzeta-1), STAT=ierrmat1)
          ALLOCATE(MATVER(nzeta-1,nzeta-1), STAT=ierrmatver)
          ALLOCATE(MATRIXINV(nzeta-1,nzeta-1), STAT=ierrmatrixInversion)
          ALLOCATE(RESID(nzeta-1,nzeta-1), STAT=ierresid)
          ALLOCATE(UNITMAT(nzeta-1,nzeta-1), STAT=ierrunit)
          ALLOCATE(AUXMATR(nzeta-1, 1), STAT=ierraux)
  
          AMATRIX = 0.
          BMATRIX = 0.
          CMATRIX = 0.
  
          IF ((ierra /= 0) .OR. (ierrb /=0) .OR. (ierrc /= 0) .OR. (ierrd /= 0)) &
               STOP "A,B,C,D : Allocation failed"
  
          DO i = 1, nzeta-3
             DMATRIX(i) = vecx(k, jz, i+1)
          END DO
  
          tpi2 = twopi_d ! yy(k,jz,nzeta) - yy(k,jz,1)
          tpi3 = twopi_d ! yy(k,jz,nzetap) - yy(k,jz,2)
  
          DMATRIX(nzeta-2) = vecx(k, jz, nzeta-1) -tpi2*vec7(k, jz, nzeta-1) &
               -tpi2*vec8(k,jz,nzeta-1)-tpi2*vec9(k,jz,nzeta-1)
          DMATRIX(nzeta-1) = vecx(k, jz, nzeta) + tpi2*vecd(k, jz, nzeta) &
               - tpi2*vec4(k, jz, nzeta) - &
               tpi2*vec6(k, jz, nzeta) -tpi3*vec7(k, jz, nzeta) &
               -tpi3*vec8(k, jz, nzeta) &
               -tpi3*vec9(k, jz, nzeta)
  
          FMATRIX = FMATRIXFULL(k-1,:)
          EMATRIX = EMATRIXFULL(k-1,:,:)
  
          jloop: DO j = 1, nzeta-1
             iloop: DO i = 1, nzeta-1
                IF (i == j)  THEN
                   AMATRIX(i,j) = vec3(k, jz, i+1)
                   BMATRIX(i,j) = -vec2(k, jz, i+1)
                   CMATRIX(i,j) = vec1(k, jz, i+1)
                END IF
                IF (j == i+1) THEN
                   AMATRIX(i,j) = vec6(k, jz, i+1)
                   BMATRIX(i,j) = vecd(k, jz, i+1)
                   CMATRIX(i,j) = vec4(k, jz, i+1)
                END IF
                IF (j == i+2) THEN
                   AMATRIX(i,j) = vec9(k, jz, i+1)
                   BMATRIX(i,j) = -vec8(k, jz, i+1)
                   CMATRIX(i,j) = vec7(k, jz, i+1)
                END IF
                IF (i == nzeta-2) THEN
                   AMATRIX(i,1) = vec9(k, jz, nzeta-1)
                   BMATRIX(i,1) = -vec8(k, jz, nzeta-1)
                   CMATRIX(i,1) = vec7(k, jz, nzeta-1)
                END IF
                IF(i == nzeta-1) THEN
                   AMATRIX(i,1) = vec6(k, jz, nzeta)
                   AMATRIX(i,2) = vec9(k, jz, nzeta)
                   BMATRIX(i,1) = vecd(k, jz, nzeta)
                   BMATRIX(i,2) = -vec8(k, jz, nzeta)
                   CMATRIX(i,1) = vec4(k, jz, nzeta)
                   CMATRIX(i,2) = vec7(k, jz, nzeta)
                END IF
             END DO iloop
          END DO jloop
  
          ! AMATRIX, BMATRIX, CMATRIX, DMATRIX, EMATRIXk, FMATRIXk  &
          ! are initialized at this point
  
          !Time to do the inverse of (Bk - Ck*Ek-1)
  
          MATRIX1 = BMATRIX - MATMUL(CMATRIX, EMATRIX)
          MATRIX2 = MATMUL(CMATRIX,FMATRIX) - DMATRIX
  
          AUXMATR = 100._dp
  
          MATVER = MATRIX1
          UNITMAT = 0._dp
  
          DO i = 1, nzeta-1
             UNITMAT(i,i) = 1._dp
          END DO
  
          ! Now time to invert MATRIX1
  
          ! CALL gaussj(MATRIX1, nzeta-1, nzeta-1, AUXMATR, 1, 1, nmaximum) 
          CALL gaussj(MATRIX1, AUXMATR)  ! f90 routine (NR) with full pivoting
  
          ! CALL matrixInversion(nzeta-1, 0, MATRIX1, auxmatr, det)  ! French
          ! routine (fast)
  
          ! Now MATRIX1 is the inverse, do iterative improvement of it 
  
          RESID = MATMUL(MATVER, MATRIX1)
          MATRIXINV = MATRIX1 + MATMUL(MATRIX1, UNITMAT-RESID)
  
          ! MATRIXINV is the iterative improvement of inv. of MATRIX1 
  
          EMATRIX = MATMUL(MATRIXINV, AMATRIX)
          FMATRIX = MATMUL(MATRIXINV, MATRIX2)
  
          ! Now we define the full E, F for this order (k)
  
          FMATRIXFULL(k,:) = FMATRIX
          EMATRIXFULL(k,:,:) = EMATRIX
  
          ! Now we have what we need, E and F for this order (k), 
          ! So we can de-allocate the matrices we defined here.
  
          IF (ALLOCATED(AMATRIX)) DEALLOCATE(AMATRIX,STAT=ierrdamat)
          IF (ALLOCATED(BMATRIX)) DEALLOCATE(BMATRIX,STAT=ierrdbmat)
          IF (ALLOCATED(CMATRIX)) DEALLOCATE(CMATRIX,STAT=ierrdcmat)
          IF (ALLOCATED(DMATRIX)) DEALLOCATE(DMATRIX,STAT=ierrddmat)
          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierrdemat)
          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierrdfmat)
          IF (ALLOCATED(MATRIX1)) DEALLOCATE(MATRIX1,STAT=ierrdm1mat)
          IF (ALLOCATED(MATVER)) DEALLOCATE(MATVER,STAT=ierrdmatver)
          IF (ALLOCATED(MATRIXINV)) DEALLOCATE(MATRIXINV,STAT=ierrdmatrixInversion)
          IF (ALLOCATED(RESID)) DEALLOCATE(RESID,STAT=ierresid)
          IF (ALLOCATED(UNITMAT)) DEALLOCATE(UNITMAT,STAT=ierrdunitmat)
          IF (ALLOCATED(MATRIX2)) DEALLOCATE(MATRIX2,STAT=ierrdm2mat)
          IF (ALLOCATED(AUXMATR)) DEALLOCATE(AUXMATR, STAT=ierrdaux)
  
          IF (k == nthe-1) EXIT kloop   ! We only need E, F up to the order N-1
  
          k = k + 1
  
       END DO kloop
  
       ! Now it's time to start calculating the U vectors, keeping in mind that 
       ! U(N) is known from the BCs
  
       k = nthe - 1
  
       kloop2: DO
  
          ALLOCATE(UVECTOR(nzeta-1), STAT=ierruv)
          ALLOCATE(UVECTORPLUS(nzeta-1), STAT=ierruvp)
          ALLOCATE(EMATRIX(nzeta-1, nzeta-1), STAT = ierre)
          ALLOCATE(FMATRIX(nzeta-1), STAT = ierrf)
  
          FMATRIX = FMATRIXFULL(k,:)
          UVECTORPLUS = UVECTORFULL(k+1,:)
          EMATRIX = EMATRIXFULL(k,:,:)
  
          UVECTOR = MATMUL(EMATRIX, UVECTORPLUS) + FMATRIX
  
          UVECTORFULL(k,:) = UVECTOR
  
          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierr)
          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierr)
          IF (ALLOCATED(UVECTOR)) DEALLOCATE(UVECTOR,STAT=ierr)
          IF (ALLOCATED(UVECTORPLUS)) DEALLOCATE(UVECTORPLUS,STAT=ierr)
  
          IF (k == 2) EXIT kloop2
  
          k = k - 1
  
       END DO kloop2
  
       ! Now all data is in the matrix UVECTORFULL
  
       ! De-allocate Efull, Ffull
  
       IF (ALLOCATED(EMATRIXFULL)) DEALLOCATE(EMATRIXFULL,STAT=ierr)
       IF (ALLOCATED(FMATRIXFULL)) DEALLOCATE(FMATRIXFULL,STAT=ierr)
  
       alfa(1:nthe,jz,1:nzeta-1) = UVECTORFULL(1:nthe,1:nzeta-1)
  
       DO k = 1, nthe
          tpi2 = twopi_d
          tpi3 = twopi_d
  
          alfa(k, jz, nzeta) = alfa(k, jz, 1) + tpi2
          alfa(k, jz, nzeta+1) = alfa(k, jz, 2) + tpi3
       END DO
  
       IF (ALLOCATED(UVECTORFULL)) DEALLOCATE(UVECTORFULL,STAT=ierr)
  
    END DO fluxloop
  
    DO k = 2,nzeta
       DO i = 2,nthem
          CALL extap(alfa(i,4,k),alfa(i,3,k),alfa(i,2,k),alfa(i,1,k))
          CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsim,k), &
                     alfa(i,npsi,k))
       END DO
    END DO
  
    DO  j = left,right
       DO  i = 1,nthe
          DO k = 2, nzeta
             alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1. - blendAlpha)
          END DO
  
          !  periodic boundary condition in delta, alfa = phi + delta
          alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
          alfa(i,j,nzetap) = alfa(i,j,2) + 2._dp * pi_d
       END DO
    END DO
  
    sumdb = SUM(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
    sumb = SUM(ABS(alfa(2:nthem,1:npsi,2:nzeta)))
    diffmx = MAXVAL(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
  
    maximum_location =  MAXLOC(ABS(alfa(2:nthem,1:npsi,2:nzeta) - alfaprev(2:nthem,1:npsi,2:nzeta)))
  
    izmx = maximum_location(1)
    jzmx = maximum_location(2)
    kmx = maximum_location(3)
  
    ierr = 2
  
    DEALLOCATE(alfaPrev, stat = idealerr)
  
    RETURN
  
  END SUBROUTINE directAlpha
  
!==============================================================================
  SUBROUTINE iterateAlpha
    !   iteration procedure to obtain matrix inversion for alfa Eqn.
  
    USE ModScbMain,      ONLY: nimax
    use ModScbParams,    ONLY: isSORDetailNeeded, iWantAlphaExtrapolation, &
                               InConAlpha
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, nzeta, nzetap, ny
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecx, &
                               left, right
  
    use nrmod,  ONLY: polint
    use nrtype, ONLY: DP, pi_d
  
    IMPLICIT NONE
  
    REAL(DP) :: om(ny), anorm, diff, rjac, &
         & anormaverage, dyDummy, anormResid, anormError, anormf
    REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:), alfaPrevTemp(:,:,:)
    REAL(DP) :: RESULT, omegaOpt
    REAL(DP), ALLOCATABLE :: angle(:,:,:), resid(:,:,:)
    INTEGER :: j, i, ni, ict, jz, k, kp, km, iz, im, ip, izmx, jzmx, kmx, ierr, nratio, &
         myPsiBegin, myPsiEnd, my_array_type2, psiRangeDiff, resultInt, loc(3)
  !  LOGICAL isnand ! Intrinsic in PGF
  
    ALLOCATE(angle(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaprev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(alfaPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)
 
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
    nisave = 0
    resid = 0

    psiloop: DO  jz = left+1, right-1
       ni = 1
       anormf = 0._dp
       DO k = 2, nzeta
          DO iz = 2, nthem
             anormf = anormf + ABS(vecx(iz,jz,k))
          END DO
       END DO
  
       Iterations: DO WHILE (ni <= nimax)
  
          anormResid = 0._dp

          zetaloop: DO k = 2, nzeta
             kp = k + 1
             km = k - 1
             thetaloop: DO  iz = 2, nthem
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
                anormResid = anormResid + ABS(resid(iz,jz,k)) ! Residue will decrease with iterations
                alfa(iz,jz,k) = alfa(iz,jz,k) + om(jz) * (resid(iz,jz,k) / vecd(iz,jz,k))
                if (alfa(iz,jz,k)+1.0.eq.alfa(iz,jz,k)) then
                   write(*,*) iz, jz, k
                   write(*,*) x(iz,jz,k), y(iz,jz,k), z(iz,jz,k)
                   call CON_stop('NaN encountered in ModScbEuler iterateAlpha')
                endif
             END DO thetaloop
          END DO zetaloop

          om(jz) = omegaOpt
          IF (isSorDetailNeeded == 1) THEN
             IF (MOD(ni, 50) == 1 .AND. jz == 2) THEN
                WRITE(*, '(A, 1X, I2, 1X, I4, A1, 1X, 3(E12.3))') &
                     'ni, Residue for alpha, RHS, residue/RHS at iteration ni = ', ni, ':',  &
                     anormResid, anormf, anormResid/anormf
             END IF
          END IF
  
          ! Check for inner loop convergence of Alpha
          IF (maxval(abs(resid(2:nthem,jz,2:nzeta))).lt.InConAlpha) then
             EXIT Iterations
          END IF
  
          ni  =   ni + 1
  
       END DO Iterations
       nisave = MAX(nisave, ni) ! We want to know the maximum number of iterations for each processor
    END DO Psiloop
 
    sumdb = SUM(ABS(alfa(2:nthem,2:npsi-1,2:nzeta) - alfaprev(2:nthem,2:npsi-1,2:nzeta)))
    sumb = SUM(ABS(alfa(2:nthem,2:npsi-1,2:nzeta)))
    diffmx = maxval(abs(resid(2:nthem,2:npsi-1,2:nzeta)))
  
    !  boundary condition at themin and themax delta = 0 and phi = alphaVal
    DO k = 2, nzeta
       DO j = 1,npsi
          alfa(1,j,k) = alphaVal(k)
          alfa(nthe,j,k) = alphaVal(k)
       END DO
       !  extrapolate alfa to the j = 1 & npsi surfaces
       DO i = 1, nthe
          do j = right,npsi
             IF (iWantAlphaExtrapolation == 0) THEN
                alfa(i,j,k) = alfa(i,j-1,k) ! If the extrapolation is problematic
             ELSE
                CALL extap(alfa(i,j-3,k),alfa(i,j-2,k),alfa(i,j-1,k),alfa(i,j,k))
             END IF
          enddo
          do j = left,1
             CALL extap(alfa(i,j+3,k),alfa(i,j+2,k),alfa(i,j+1,k),alfa(i,j,k)) ! This is never a problem - very close to Earth
          enddo
       END DO
    END DO

    !...  set "blending" in alpha for outer iteration loop
    DO  j = 1, npsi
       DO  i = 1, nthe
          DO  k = 2, nzeta
             alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1._dp - blendAlpha)
          END DO
          alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
          alfa(i,j,nzetap) = alfa(i,j,2) + 2._dp * pi_d
       END DO
    END DO

    IF (ALLOCATED(angle)) DEALLOCATE(angle, STAT = ierr)
    IF (ALLOCATED(alfaPrev)) DEALLOCATE(alfaPrev)
    IF (ALLOCATED(alfaPrevTemp)) DEALLOCATE(alfaPrevTemp)
    IF (ALLOCATED(resid)) DEALLOCATE(resid)
 
    RETURN
  
  END SUBROUTINE iterateAlpha
  
!================================================!
!============ Psi Euler Potential ===============!
!================================================!
  SUBROUTINE psiFunctions
  
    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: npsi
    use ModScbVariables, ONLY: rhoVal, f, fp
  
    USE ModScbSpline, ONLY: Spline_derivs_1D
  
    IMPLICIT NONE
  
    REAL(DP) :: psiValue(npsi)
    INTEGER :: j
  
    CALL Spline_derivs_1D(rhoVal(1:npsi), psival(1:npsi), psiValue, f(1:npsi))
    CALL Spline_derivs_1D(rhoVal(1:npsi), f(1:npsi), psiValue, fp(1:npsi))
  
    RETURN
  END SUBROUTINE psiFunctions

!==============================================================================
  SUBROUTINE InterpolatePsiR
    ! Interpolates the new values of psi at the new locations xnew on midnight 
    ! equator
    !!!! Module Variables
    USE ModScbParams,  ONLY: iAzimOffset
    USE ModScbGrids, ONLY: npsi
    use ModScbVariables, ONLY: x, y, nThetaEquator, nZetaMidnight, psiVal, kmax, &
                               radEqmidNew
    !!!! Module Subroutines/Functions
    USE ModScbSpline, ONLY: spline, splint
    !!!! NR Modules
    use nrtype, ONLY: DP

    IMPLICIT NONE

    REAL(DP), DIMENSION(npsi) :: radEqmid, psival1D, psi2Deriv
    INTEGER :: ialloc, ierr, j

    IF (iAzimOffset == 1) THEN
       radEqmid(1:npsi) = SQRT(x(nThetaEquator, 1:npsi, nZetaMidnight)**2 + &
            y(nThetaEquator, 1:npsi, nZetaMidnight)**2)
    ELSE IF (iAzimOffset == 2) THEN
       radEqmid(1:npsi) = SQRT(x(nThetaEquator, 1:npsi, kMax)**2 + &
            y(nThetaEquator, 1:npsi, kMax)**2)
    END IF

    psival1D(1:npsi) = psival(1:npsi)

    CALL spline(radEqmid, psival1D, 1.E31_dp, 1.E31_dp, psi2Deriv)

    DO j = 2, npsi-1
       psival(j) = splint(radEqmid, psival1D, psi2Deriv, radEqmidNew(j))
    END DO

    RETURN

  END SUBROUTINE InterpolatePsiR

!==============================================================================
  SUBROUTINE psiges
    !   initial guess of poloidal flux
  
    USE ModScbGrids,     ONLY: npsi, nthe, nzetap
  
    INTEGER :: i, j, k
  
    !     initial guess for psi
    DO  k=1,nzetap
       DO  j=1,npsi
          DO  i=1,nthe
             psi(i,j,k)=psival(j)
          END DO
       END DO
    END DO
  
    RETURN
  
  END SUBROUTINE psiges
  
!==============================================================================
  SUBROUTINE mapPsi(iSmoothMove)
    ! new psiMap, without computing linear distance
    ! first and last flux surface remain unchanged
  
    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, psiPrev, &
                               left, right
  
    USE ModScbSpline, ONLY: spline, splint
  
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: iSmoothMove
    INTEGER :: ierr, k, i, j
    REAL(DP) :: distance1deriv1, distance1derivn
    REAL(DP), DIMENSION(npsi) :: xOld, yOld, zOld, psiOld, distanceOld, x2derivs, y2derivs, &
                                 z2derivs
    REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
    REAL(DP) :: blend
  
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
          iloop: DO i = 1, nthe
             xOld(1:npsi) = x(i,1:npsi,k)
             yOld(1:npsi) = y(i,1:npsi,k)
             zOld(1:npsi) = z(i,1:npsi,k)
             psiOld(1:npsi) = psi(i,1:npsi,k)
  
             CALL spline(psiold, xOld, 1.E31_dp, 1.E31_dp, x2derivs)
             CALL spline(psiold, yOld, 1.E31_dp, 1.E31_dp, y2derivs)
             CALL spline(psiold, zOld, 1.E31_dp, 1.E31_dp, z2derivs)
  
             DO j = 2, npsi-1
                x(i,j,k) = splint(psiold, xOld, x2derivs, psival(j))
                y(i,j,k) = splint(psiold, yOld, y2derivs, psival(j))
                z(i,j,k) = splint(psiold, zOld, z2derivs, psival(j))
                psi(i,j,k) = psival(j)
                !x(i,j,k) = splint(psiold, xOld, x2derivs, psiPrev(i,j,k))
                !y(i,j,k) = splint(psiold, yOld, y2derivs, psiPrev(i,j,k))
                !z(i,j,k) = splint(psiold, zOld, z2derivs, psiPrev(i,j,k))
             END DO
          END DO iloop
       END DO kloop
  
       ! periodic boundary condition in zeta
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzetap) = x(:,:,2)
       y(:,:,nzetap) = y(:,:,2)
       z(:,:,nzetap) = z(:,:,2)
  
       DO j = 1, npsi
          psi(:,j, 1) = psival(j)
          psi(:,j,nzetap) = psival(j)
       END DO
  
    END IF
  
    RETURN
  
  END SUBROUTINE mapPsi
  
!==============================================================================
  SUBROUTINE directPsi
  
    USE ModScbMain,      ONLY: DP
    USE ModScbGrids,     ONLY: nzeta, nzetap, nthe, nthem, npsi, npsim
    use ModScbVariables, ONLY: sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecr, &
                               left, right
 
    use nrmod, ONLY: gaussj 
    IMPLICIT NONE
  
    REAL(DP), ALLOCATABLE :: psiprev(:,:,:)
    INTEGER :: ierra, ierrb, ierrc, ierrd, i, j, k
    REAL(DP), DIMENSION(:, :), ALLOCATABLE :: FMATRIXFULL, UVECTORFULL
    REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: EMATRIXFULL
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: AMATRIX, BMATRIX, CMATRIX, EMATRIX, &
         MATRIX1, MATVER, MATRIXINV, RESID
    REAL(DP), DIMENSION(:), ALLOCATABLE :: DMATRIX, FMATRIX, MATRIX2, UVECTOR, &
         UVECTORPLUS
  
    REAL(DP), ALLOCATABLE :: AUXMATR(:,:)  !!! DIMENSION(nthe-2,1) :: AUXMATR
    REAL(DP), ALLOCATABLE :: UNITMAT(:,:)  !!! :: UNITMAT(nthe-2,nthe-2)
    REAL(DP) :: diff
  
    INTEGER, DIMENSION(3) :: maximum_location
  
    INTEGER :: jz, ierrffull, ierrefull,  ierre, ierrf, &
         ierrmat1, ierrmatver, ierrmatrixinversion, ierresid, ierrdamat, ierrdbmat, &
         ierrdcmat, ierrddmat, ierrdemat, ierrdfmat, ierrdm1mat, ierrmdmatver, ierrdmatrixinversion, &
         ierrdunitmat, ierrdm2mat, ierruv, ierruvp, ierr, izmx, jzmx, kmx, ierrdmatver, ierrmatv, kzeta, &
         ierrmatres, ierrdmver, ierrdresid, iz, ierraux, ierrunit, ierrdaux, ierrdunit, idealerr
  
    ALLOCATE(psiPrev(SIZE(psi,1), SIZE(psi,2), SIZE(psi,3)), stat = ierr)
  
    psiprev(:,:,:) = psi(:,:,:)
  
    sumdb = 0._dp
    sumb = 0._dp
    diffmx = 0._dp
  
    ! The subroutine for the psi equation is a bit different from that for alpha;
  
    ! Start here the loop for different alpha surfaces
  
    zetaloop: DO kzeta = 2, nzeta
  
       ! PRINT*, 'Alpha surface:', kzeta
  
       ALLOCATE(FMATRIXFULL(npsi, nthe-2), STAT=ierrffull)
       ALLOCATE(EMATRIXFULL(npsi, nthe-2, nthe-2), STAT=ierrefull)
       ALLOCATE(UVECTORFULL(npsi, nthe), STAT=ierrffull)
  
       EMATRIXFULL(1,:,:) = 0.   ! first index is the flux surface
  
       UVECTORFULL(1,:) = psival(1)
       UVECTORFULL(npsi,:) = psival(npsi)
  
       FMATRIXFULL(1,:) = psival(1)
  
       k = 2    ! k are the flux surfaces here
  
       kloop: DO 
  
          ALLOCATE(AMATRIX(nthe-2,nthe-2), STAT=ierra)
          ALLOCATE(BMATRIX(nthe-2,nthe-2), STAT=ierrb)
          ALLOCATE(CMATRIX(nthe-2,nthe-2), STAT=ierrc)
          ALLOCATE(DMATRIX(nthe-2), STAT=ierrd)
          ALLOCATE(EMATRIX(nthe-2,nthe-2), STAT=ierre)
          ALLOCATE(FMATRIX(nthe-2), STAT=ierrf)
          ALLOCATE(MATRIX1(nthe-2,nthe-2), STAT=ierrmat1)
          ALLOCATE(MATVER(nthe-2,nthe-2), STAT=ierrmatv)
          ALLOCATE(MATRIXINV(nthe-2,nthe-2), STAT=ierrmatrixInversion)
          ALLOCATE(RESID(nthe-2,nthe-2), STAT=ierrmatres)
          ALLOCATE(MATRIX2(nthe-2), STAT=ierrmat1)
          ALLOCATE(AUXMATR(nthe-2, 1), STAT=ierraux)
          ALLOCATE(UNITMAT(nthe-2, nthe-2), STAT=ierrunit)
  
          AMATRIX = 0.
          BMATRIX = 0.
          CMATRIX = 0.
  
          IF ((ierra /= 0) .OR. (ierrb /=0) .OR. (ierrc /= 0) .OR. (ierrd /= 0)) &
               STOP "A,B,C,D : Allocation failed"
  
          DO i = 1, nthe-2
             DMATRIX(i) = vecr(i+1, k, kzeta)
  
             IF (i==1) THEN
                DMATRIX(i) = vecr(i+1, k, kzeta) - &
                     vec1(i+1, k, kzeta)*psi(i,k-1,kzeta) &
                     -vec4(i+1, k, kzeta)*psi(i,k,kzeta) &
                     -vec7(i+1, k, kzeta)*psi(i,k+1,kzeta)
             END IF
  
             IF (i==nthe-2) THEN
                DMATRIX(i) = vecr(i+1, k, kzeta) - &
                     vec3(i+1, k, kzeta)*psi(i,k-1,kzeta) &
                     -vec6(i+1, k, kzeta)*psi(i,k,kzeta) &
                     -vec9(i+1, k, kzeta)*psi(i,k+1,kzeta)
             END IF
  
          END DO
  
          EMATRIX(:,:) = EMATRIXFULL(k-1,:,:)
          FMATRIX(:) = FMATRIXFULL(k-1,:) ! FMATRIX of order k-1 
  
          jloop: DO j = 1, nthe-2
  
             iloop: DO i = 1, nthe-2
  
                IF (i == j)  THEN
                   AMATRIX(i,j) = vec8(i+1, k, kzeta)
                   BMATRIX(i,j) = vecd(i+1, k, kzeta)
                   CMATRIX(i,j) = vec2(i+1, k, kzeta)
                END IF
  
                IF (j == i+1) THEN
                   AMATRIX(i,j) = vec9(i+1, k, kzeta)
                   BMATRIX(i,j) = -vec6(i+1, k, kzeta)
                   CMATRIX(i,j) = vec3(i+1, k, kzeta)
                END IF
  
                IF (j == i-1) THEN
                   AMATRIX(i,j) = vec7(i+1, k, kzeta)
                   BMATRIX(i,j) = -vec4(i+1, k, kzeta)
                   CMATRIX(i,j) = vec1(i+1, k, kzeta)
                END IF
  
             END DO iloop
          END DO jloop
  
  
          ! AMATRIX, BMATRIX, CMATRIX, DMATRIX, EMATRIXk, FMATRIXk  are initialized &
          ! at this point
  
          !Time to do the inverse of (Bk - Ck*Ek-1)
  
          MATRIX1 = BMATRIX - MATMUL(CMATRIX, EMATRIX)
          MATRIX2 = MATMUL(CMATRIX,FMATRIX) - DMATRIX
  
          ! Now time to invert MATRIX1
  
          AUXMATR = 100.   ! just a test value, only the inverse of matrix1 interests us
  
          UNITMAT = 0.
  
          DO i = 1,nthe-2
             UNITMAT(i,i) = 1.
          END DO
  
          MATVER = MATRIX1  ! MATVER keeps matrix1
  
          ! CALL gaussj(MATRIX1, nthe-2, nthe-2, AUXMATR, 1, 1, nmaximum)    
  
          CALL gaussj(MATRIX1,AUXMATR)   ! f90 routine (NR) with pivoting
  
          ! CALL matrixInversion(nthe-2, 0, MATRIX1, auxmatr, det)  ! French routine
  
          ! Now on output MATRIX1 should be the inverse
  
          ! Start iterative improvement for the solution
  
          RESID = MATMUL(MATVER, MATRIX1)
  
          MATRIXINV = MATRIX1 + MATMUL(MATRIX1, UNITMAT-RESID)
  
          ! MATRIXINV is the corrected value for the inverse
  
          EMATRIX = MATMUL(MATRIXINV, AMATRIX)
          FMATRIX = MATMUL(MATRIXINV, MATRIX2)
  
          ! Now we define the full E, F for this order (k)
  
          FMATRIXFULL(k, :) = FMATRIX(:)
          EMATRIXFULL(k, :, :) = EMATRIX(:, :)
  
          ! Now we have what we need, E and F for this order (k), 
          ! So we can de-allocate the matrices we defined here.
  
          IF (ALLOCATED(AMATRIX)) DEALLOCATE(AMATRIX,STAT=ierrdamat)
          IF (ALLOCATED(BMATRIX)) DEALLOCATE(BMATRIX,STAT=ierrdbmat)
          IF (ALLOCATED(CMATRIX)) DEALLOCATE(CMATRIX,STAT=ierrdcmat)
          IF (ALLOCATED(DMATRIX)) DEALLOCATE(DMATRIX,STAT=ierrddmat)
          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierrdemat)
          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierrdfmat)
          IF (ALLOCATED(MATRIX1)) DEALLOCATE(MATRIX1,STAT=ierrdm1mat)
          IF (ALLOCATED(MATVER)) DEALLOCATE(MATVER,STAT=ierrdmver)
          IF (ALLOCATED(MATRIXINV)) DEALLOCATE(MATRIXINV,STAT=ierrmatrixInversion)
          IF (ALLOCATED(RESID)) DEALLOCATE(MATVER,STAT=ierrdresid)
          IF (ALLOCATED(MATRIX2)) DEALLOCATE(MATRIX2,STAT=ierrdm2mat)
          IF (ALLOCATED(AUXMATR)) DEALLOCATE(AUXMATR, STAT=ierrdaux)
          IF (ALLOCATED(UNITMAT)) DEALLOCATE(UNITMAT, STAT=ierrdunit)
  
          IF (k == npsi - 1) EXIT kloop   ! We only need E, F up to the order N-1
  
          k = k + 1
  
       END DO kloop
  
       ! Now it's time to start calculating the U vectors, keeping in mind that 
       ! U(N) is known from the BCs
  
       k = npsi - 1
  
       kloop2: DO
  
          ALLOCATE(UVECTOR(nthe-2), STAT=ierruv)
          ALLOCATE(UVECTORPLUS(nthe-2), STAT=ierruvp)
          ALLOCATE(EMATRIX(nthe-2, nthe-2), STAT = ierre)
          ALLOCATE(FMATRIX(nthe-2), STAT = ierrf)
  
          FMATRIX(:) = FMATRIXFULL(k,:)
          EMATRIX(:,:) = EMATRIXFULL(k,:,:)
  
          DO i = 1, nthe-2
             UVECTORPLUS(i) = UVECTORFULL(k+1,i+1)
          END DO
  
          UVECTOR = MATMUL(EMATRIX, UVECTORPLUS) + FMATRIX
  
          DO i = 1, nthe-2
             UVECTORFULL(k, i+1) = UVECTOR(i)
          END DO
  
          IF (ALLOCATED(EMATRIX)) DEALLOCATE(EMATRIX,STAT=ierr)
          IF (ALLOCATED(FMATRIX)) DEALLOCATE(FMATRIX,STAT=ierr)
          IF (ALLOCATED(UVECTOR)) DEALLOCATE(UVECTOR,STAT=ierr)
          IF (ALLOCATED(UVECTORPLUS)) DEALLOCATE(UVECTORPLUS,STAT=ierr)
  
          IF (k == 2) EXIT kloop2
  
          k = k-1
  
       END DO kloop2
  
       ! Now all data is in the matrix UVECTORFULL
       ! De-allocate efull, ffull
  
       IF (ALLOCATED(EMATRIXFULL)) DEALLOCATE(EMATRIXFULL,STAT=ierr)
       IF (ALLOCATED(FMATRIXFULL)) DEALLOCATE(FMATRIXFULL,STAT=ierr)
  
       DO k = 1, npsi
          DO i = 2, nthe-1
             psi(i, k, kzeta) = UVECTORFULL(k,i)
          END DO
       END DO
  
       ! the values for theta = 1 and NSIZE are the same as defined in psiges
       ! Now psi is completely determined for zeta surface kzeta = kzeta
  
       IF (ALLOCATED(UVECTORFULL)) DEALLOCATE(UVECTORFULL,STAT=ierr)
  
    END DO zetaloop
  
    !...  set blend for outer iteration loop
  
    DO  j = 1,npsi
       DO  i = 1,nthe
          DO k = 2, nzeta
             psi(i,j,k) = psi(i,j,k) * blendPsi + psival(j) * (1. - blendPsi)  
             !  blend does the "blending" between computed psi at (n+1) and at (n)
          END DO
          psi(i,j,1) = psi(i,j,nzeta)
          psi(i,j,nzetap) = psi(i,j,2)
       END DO
    END DO
  
    DO k = 2, nzeta
       DO jz = 2, npsim
          DO iz = 2, nthem
             diff = ABS(psi(iz,jz,k)-psiprev(iz,jz,k))
             IF(diff.GT.diffmx) diffmx = diff
             sumdb = sumdb + ABS(psi(iz,jz,k)-psiprev(iz,jz,k))
             sumb = sumb + ABS(psi(iz,jz,k))
          END DO
       END DO
    END DO
  
    DEALLOCATE(psiPrev, stat = idealerr)
  
    RETURN
  
  END SUBROUTINE directPsi
  
!==============================================================================
  SUBROUTINE iteratePsi
    !   iteration procedure to obtain matrix inversion
  
    USE ModScbMain,      ONLY: DP, nimax
    use ModScbParams,    ONLY: isSORDetailNeeded, InConPsi
    USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, na
    USE ModScbVariables, ONLY: nisave, x, y, z, sumb, sumdb, vecd, vec1, vec2, &
                               vec3, vec4, vec6, vec7, vec8, vec9, vecr, &
                               left, right
  
    use nrtype, ONLY: pi_d
  
    IMPLICIT NONE
  
    REAL(DP) :: om(na)
    REAL(DP) :: omegaOpt
    REAL(DP) :: omc, anorm, anormf, diff, ano, sumbtest, anormResid, anormError
    REAL(DP), ALLOCATABLE :: psiPrev(:,:,:), psiPrevTemp(:,:,:), resid(:,:,:)
    INTEGER :: j, k, i, ierr, ni, ict, jz, jp, jm, iz, im, ip, nratio, myAlphaBegin, myAlphaEnd, &
         my_array_type_psi2, alphaRangeDiff, resultInt, loc(3)
  
    !     perform SOR iteration
    !..   choose for inner iteration loop
  
    ALLOCATE(psiPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(psiPrev(nthe,npsi,nzeta+1), STAT = ierr)
    ALLOCATE(resid(nthe,npsi,nzeta+1), STAT = ierr)

    rjac = 1._dp - 2._dp*pi_d*pi_d / (REAL(nthe,dp)*REAL(nthe,dp) + REAL(npsi,dp)*REAL(npsi,dp)) ! Radius of conv. of Jacobi iteration,
    ! could be used to find omega optimal in SOR
    !rjac = (cos(pi_d/real(npsi,dp)) + (dr/dt)**2 * cos(pi_d/real(nthe,dp))) / &
    !     (1._dp +  (dr/dt)**2)
  
    om = 1._dp
    omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))
    psiPrev = psi
    nisave = 0
    resid = 0

    alphaLoop: DO k = 2, nzeta
       ni = 1
       anormf = 0._dp
       DO jz = 2, npsi-1
          DO iz = 2, nthem
             anormf = anormf + ABS(vecr(iz,jz,k))
          END DO
       END DO
  
       Iterations: DO WHILE (ni <= nimax)
          anormResid = 0._dp
          jLoop: DO jz = left+1, right-1
             jp = jz + 1
             jm = jz - 1
             iLoop: DO  iz = 2, nthem
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
                anormResid = anormResid + ABS(resid(iz,jz,k))
                psi(iz,jz,k) = psi(iz,jz,k) + om(k)*resid(iz,jz,k)/vecd(iz,jz,k)
                if (psi(iz,jz,k)+1.0.eq.psi(iz,jz,k)) then
                   call CON_stop('NaN encountered in ModScbEuler iteratePsi')
                endif
             END DO iLoop
          END DO jLoop

          om(k) = omegaOpt
  
          IF (isSorDetailNeeded == 1) THEN
             IF(MOD(ni,50) == 1 .AND. k == 2) THEN
                WRITE(*, '(A, 1X, I3, 1X, I4, A1, 1X, 3(E12.3))') 'ni, Residue for psi, RHS , residue/RHS at ni = ', &
                     ni, ':', anormResid, anormf, anormResid/anormf
             END IF
          END IF

          ! Check for inner loop Psi convergence
          IF (maxval(abs(resid(2:nthem,2:npsi-1,k))).lt.InConPsi) then
             EXIT Iterations
          END IF
  
          ni = ni + 1
  
       END DO Iterations
       ! print*, 'k, ni = ', k, ni
       nisave = MAX(nisave, ni)
    END DO AlphaLoop

    sumdb = SUM(ABS(psi(2:nthem,2:npsim,2:nzeta) - psiprev(2:nthem,2:npsim,2:nzeta)))
    sumb = SUM(ABS(psi(2:nthem,2:npsim,2:nzeta)))
    diffmx = maxval(abs(resid(2:nthem,2:npsi-1,2:nzeta)))
  
    ! Set blend for outer iteration loop, and apply periodic boundary conditions
    DO  j = 1, npsi
       DO  i = 1, nthe
          DO  k = 2,nzeta
             psi(i,j,k) = psi(i,j,k) * blendPsi + psival(j) * (1._dp - blendPsi)
          END DO
          psi(i,j,1) = psi(i,j,nzeta)
          psi(i,j,nzetap) = psi(i,j,2)
       END DO
    END DO

    IF (ALLOCATED(psiPrevTemp))  DEALLOCATE(psiPrevTemp)
    IF (ALLOCATED(psiPrev)) DEALLOCATE(psiPrev)
    IF (ALLOCATED(resid)) DEALLOCATE(resid)

    RETURN
  
  END SUBROUTINE iteratePsi
  
  END MODULE ModScbEuler
