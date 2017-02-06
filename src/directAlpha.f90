SUBROUTINE directAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)

  USE nrtype
  USE Module1 

  IMPLICIT NONE

  INTERFACE gauss
     SUBROUTINE gaussj(a,b)
       USE nrtype 
       USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b  ! modification from original form, kind=8
     END SUBROUTINE gaussj
  END INTERFACE

  ! define the (NSIZE-1) vectors Fk and [(NSIZE-1)x(NSIZE-1)] matrices Ek, & 
  ! which will be kept till end

  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd, vec1, vec2, vec3, vec4, vec6, vec7, &
       vec8, vec9, vecx

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

  fluxloop: DO jz = 2, npsim   ! jz represents the flux surface

     ALLOCATE(FMATRIXFULL(nthe, nzeta-1), STAT=ierrffull)
     ALLOCATE(EMATRIXFULL(nthe, nzeta-1, nzeta-1), STAT=ierrefull)
     ALLOCATE(UVECTORFULL(nthe, nzeta-1), STAT=ierrffull)

     EMATRIXFULL(1,:,:) = 0._dp 

     UVECTORFULL(1,1:nzeta-1) = alphaVal(1:nzeta-1)     ! BC on theta    ; here first index 
     ! is theta
     UVECTORFULL(nthe,1:nzeta-1) = alphaVal(1:nzeta-1) ! BC on theta

     FMATRIXFULL(1,:) = UVECTORFULL(1,:)


     ! for k from 2 to NSIZE-1, define Ak, Bk, Ck, Dk
     ! after each "step" k, these can be de-allocated, and only Ekfull, Fkfull kept
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
        tpi3 = twopi_d ! yy(k,jz,nzetp) - yy(k,jz,2)

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

        ! CALL matrixInversion(nzeta-1, 0, MATRIX1, auxmatr, det)  ! French routine (fast)

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
        CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsim,k)  &
             ,alfa(i,npsi,k))
     END DO
  END DO

  DO  j = 1,npsi
     DO  i = 1,nthe
        DO k = 2, nzeta
           alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1. - blendAlpha)
        END DO

        !  periodic boundary condition in delta, alfa = phi + delta
        alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
        alfa(i,j,nzetp) = alfa(i,j,2) + 2._dp * pi_d
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
