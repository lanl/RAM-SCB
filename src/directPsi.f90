SUBROUTINE directPsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTERFACE gauss
     SUBROUTINE gaussj(a,b)
       USE nrtype 
       USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b  ! modification from original form, kind=8
       INTEGER, DIMENSION(SIZE(a,1)) :: ipiv,indxr,indxc
       LOGICAL(lgt), DIMENSION(SIZE(a,1)) :: lpiv
       REAL(DP) :: pivinv
       REAL(DP), DIMENSION(SIZE(a,1)) :: dumc
       INTEGER, TARGET :: irc(2)
       INTEGER :: i,l,n
       INTEGER,  POINTER :: irow,icol
     END SUBROUTINE gaussj
  END INTERFACE

  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd,vec1,vec2,vec3,vec4, &
       & vec6,vec7,vec8,vec9,vecr

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
        psi(i,j,nzetp) = psi(i,j,2)
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

  ! sumdb = sum(abs(psi(2:nthem,2:npsim,2:nzeta) - psiprev(2:nthem,2:npsim,2:nzeta)))
  ! sumb = sum(abs(psi(2:nthem,2:npsim,2:nzeta)))
  ! diffmx = maxval(abs(psi(2:nthem,2:npsim,2:nzeta)))
  ! maximum_location =  maxloc(abs(psi(2:nthem,2:npsim,2:nzeta) - psiprev(2:nthem,2:npsim,2:nzeta)))

  ! izmx = maximum_location(1)
  ! jzmx = maximum_location(2)
  ! kmx = maximum_location(3) 


  RETURN

END SUBROUTINE directPsi
