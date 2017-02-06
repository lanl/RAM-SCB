!******************************************************************************
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

SUBROUTINE iterateAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)
  !   iteration procedure to obtain matrix inversion for alfa Eqn.

  USE nrtype
  USE Module1
  USE mpi
  USE ModRamMpi, ONLY: iComm

  IMPLICIT NONE

  INTERFACE polinomial_interpolation  ! (and extrapolation)
     SUBROUTINE polint(xa,ya,x,y,dy)
       USE nrtype; USE nrutil, ONLY : assert_eq, iminloc, nrerror
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa, ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y, dy
     END SUBROUTINE polint
  END INTERFACE

  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vec1, vec2, vec3, vec4, vec6, vec7, &
       vec8, vec9, vecx
  REAL(DP), DIMENSION(:,:,:), INTENT(IN OUT) :: vecd
  REAL(DP) :: om(ny), anorm, resid, diff, &
       & anormaverage, dyDummy, anormResid, anormError, anormf
  REAL(DP), ALLOCATABLE :: alfaPrev(:,:,:), alfaPrevTemp(:,:,:), radTot(:,:,:)
  REAL(DP) :: RESULT, omegaOpt
  REAL(DP), ALLOCATABLE :: angle(:,:,:)
  INTEGER :: psiBeg(0:63), psiEnd(0:63), param(0:63), my_array_type(0:63)  ! Assume no more than 64 processors
  INTEGER :: sizes(3), subsizes(3), starts(3)
  INTEGER :: request(3, 0:63)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: j, i, ni, ict, jz, k, kp, km, iz, im, ip, izmx, jzmx, kmx, ierr, nratio, &
       myPsiBegin, myPsiEnd, my_array_type2, psiRangeDiff, resultInt
!  LOGICAL isnand ! Intrinsic in PGF

  ALLOCATE(radTot(nthe, npsi, nzeta+1), STAT = ierr)
  ALLOCATE(angle(nthe,npsi,nzeta+1), STAT = ierr)
  xx = SQRT(x**2 + y**2)
  radTot = SQRT(xx**2+z**2)

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
  ! omegaOpt = 2._dp / (1._dp + 2._dp*pi_d/REAL(nzeta+nthe, DP))
  omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))

  ALLOCATE(alfaprev(nthe,npsi,nzeta+1), STAT = ierr)
  ALLOCATE(alfaPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)

  alfaPrev = alfa

  ! Give each processor a chunk of the total number of psi surfaces to compute (psis - jz - are not connected here)
  nratio = CEILING(REAL(npsi-2, dp)/REAL(numProc,dp))

  IF (numProc /= 1) THEN
     IF (rank == 0) THEN
        psiBeg(rank) = 2
        psiEnd(rank) = nratio + 1
     ELSE IF (rank /= numProc-1) THEN
        psiBeg(rank) = nratio * rank + 2
        psiEnd(rank) = MIN(npsi-1, nratio * (rank+1) + 1)
     ELSEIF (rank == numProc-1) THEN
        psiBeg(rank) = MIN(npsi-1, nratio * (numProc-1) + 2)
        psiEnd(rank) = npsi-1
     END IF
  ELSE
     psiBeg(rank) = 2
     psiEnd(rank) = npsi-1
  END IF

  ! Make master know all psi ranges
  CALL MPI_BARRIER(iComm, ierr)
  myPsiBegin = psiBeg(rank)
  CALL MPI_GATHER(myPsiBegin, 1, MPI_INTEGER, param, 1, MPI_INTEGER, &
       0, iComm, ierr)

  IF (rank == 0)  psiBeg(0:numProc-1) = param(0:numProc-1)

  myPsiEnd = psiEnd(rank)
  CALL MPI_GATHER(myPsiEnd,1,MPI_INTEGER,param,1,MPI_INTEGER, &
       0, iComm, ierr)

  IF (rank == 0) psiEnd(0:numProc-1) = param(0:numProc-1)

  nisave = 0

  psiloop: DO  jz = psiBeg(rank), psiEnd(rank)
     
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
              resid = - vecd(iz,jz,k)*alfa(iz,jz,k) + &
                   vec1(iz,jz,k)*alfa(im,jz,km) + & 
                   vec2(iz,jz,k)*alfa(iz,jz,km) + &
                   vec3(iz,jz,k)*alfa(ip,jz,km) + & 
                   vec4(iz,jz,k)*alfa(im,jz,k) + &
                   vec6(iz,jz,k)*alfa(ip,jz,k) + & 
                   vec7(iz,jz,k)*alfa(im,jz,kp) + &
                   vec8(iz,jz,k)*alfa(iz,jz,kp) + & 
                   vec9(iz,jz,k)*alfa(ip,jz,kp) - vecx(iz,jz,k)
       !       if (isSorDetailNeeded) write(*, '(A, 3I4, 12F11.3)') 'itAlpha: jz, k, iz, vecs, alfa, resid = ', jz, k, iz, &
       !            vecd(iz,jz,k), vec1(iz,jz,k), vec2(iz,jz,k), vec3(iz,jz,k), vec4(iz,jz,k), &
       !            vec6(iz,jz,k), vec7(iz,jz,k), vec8(iz,jz,k), vec9(iz,jz,k), vecx(iz,jz,k), alfa(iz,jz,k), resid
       !C       IF (isnand(resid)) STOP 'itAlpha: NaNs.'
              anormResid = anormResid + ABS(resid) ! Residue will decrease with iterations
              alfa(iz,jz,k) = alfa(iz,jz,k) + om(jz) * resid / vecd(iz,jz,k)
           END DO thetaloop
        END DO zetaloop

        om(jz) = omegaOpt

        !  periodic boundary condition in delta, alfa = phi + delta
        !  DO i = 2, nthem
        !    alfa(i,jz,1) = alfa(i,jz,nzeta) - twopi_d
        !    alfa(i,jz,nzetp) = alfa(i,jz,2) + twopi_d
        !  END DO

        IF (isSorDetailNeeded == 1) THEN
           IF (MOD(ni, 50) == 1 .AND. jz == psiBeg(rank)) THEN
              WRITE(*, '(A, 1X, I2, 1X, I4, A1, 1X, 3(E12.3))') &
                   'RANK, ni, Residue for alpha, RHS, residue/RHS at iteration ni = ', rank, ni, ':',  &
                   anormResid, anormf, anormResid/anormf
              CALL FLUSH(6) ! Latter ratio is important for convergence
           END IF
        END IF

        IF((ni > 1) .AND. (anormResid < 1e-3_dp*anormf)) THEN  ! SOR iteration has converged
           EXIT Iterations
        END IF

        ni  =   ni + 1

     END DO Iterations
     nisave = MAX(nisave, ni) ! We want to know the maximum number of iterations for each processor
  END DO Psiloop

  sumdb = SUM(ABS(alfa(2:nthem,psiBeg(rank):psiEnd(rank),2:nzeta) - alfaprev(2:nthem,psiBeg(rank):psiEnd(rank),2:nzeta)))
  sumb = SUM(ABS(alfa(2:nthem,psiBeg(rank):psiEnd(rank),2:nzeta)))
  diffmx = MAXVAL(ABS(alfa(2:nthem,psiBeg(rank):psiEnd(rank),2:nzeta) - alfaprev(2:nthem,psiBeg(rank):psiEnd(rank),2:nzeta)))

  IF (ALLOCATED(alfaPrev)) DEALLOCATE(alfaPrev)
  IF (ALLOCATED(alfaPrevTemp)) DEALLOCATE(alfaPrevTemp)

  CALL MPI_BARRIER(iComm, ierr)

  CALL MPI_REDUCE(sumdb, RESULT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       0, iComm, ierr)
  IF (rank == 0) sumdb = RESULT
  CALL MPI_REDUCE(sumb, RESULT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       0, iComm, ierr)
  IF (rank == 0) sumb = RESULT
  CALL MPI_REDUCE(diffmx, RESULT, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, iComm, ierr)
  IF (rank == 0) diffmx = RESULT
  CALL MPI_REDUCE(nisave, resultInt, 1, MPI_INTEGER, MPI_MAX, &
       0, iComm, ierr)
  IF (rank == 0) nisave = resultInt ! Now we have the maximum number of iterations among all processors/flux surfaces

800 CONTINUE

  ! Have to assemble alpha on the master node from the calculations done by all processors
  CALL MPI_BARRIER(iComm, ierr)  ! necessary to make sure all data is newly computed
  ! by all processors

  sizes(1) = nthe
  sizes(2) = npsi
  sizes(3) = nzeta+1
  subsizes(1) = nthe
  subsizes(3) = nzeta+1
  starts(1) = 0
  starts(3) = 0

  IF (rank /= 0) THEN
     subsizes(2) = psiEnd(rank) - psiBeg(rank)+1
     starts(2) = psiBeg(rank) -1 ! because it starts from 0
     CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
          my_array_type(rank), ierr)
     CALL MPI_TYPE_COMMIT(my_array_type(rank), ierr)
     CALL MPI_SEND(alfa, 1, my_array_type(rank), 0, 17+3*rank, iComm, ierr)
  END IF

  IF (rank == 0 .AND. numProc > 1) THEN
     DO i = 1, numProc-1
        subsizes(2) = psiEnd(i) - psiBeg(i) + 1
        starts(2) = psiBeg(i)-1
        CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
             my_array_type(i), ierr)
        CALL MPI_TYPE_COMMIT(my_array_type(i), ierr)
        CALL MPI_IRECV(alfa, 1, my_array_type(i), i, 17+3*i, iComm, request(1, i), ierr)
     END DO
     DO i = 1, numProc-1
        DO j = 1, 1
           CALL MPI_WAIT(request(j, i), status, ierr) !C Completes the non-blocking receives above
        END DO
     END DO
  END IF
999 CONTINUE

  CALL MPI_BARRIER(iComm, ierr)

  IF (rank == 0) THEN
     DO i = 1, numProc-1
        CALL MPI_TYPE_FREE(my_array_type(i), ierr)
     END DO
  ELSE
     CALL MPI_TYPE_FREE(my_array_type(rank), ierr)
  END IF

  ! Broadcast the new alfa array to all processors
  subsizes(2) = npsi
  starts(2) = 0

  CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
       my_array_type2, ierr)
  CALL MPI_TYPE_COMMIT(my_array_type2, ierr)
  CALL MPI_BCAST(alfa, 1, my_array_type2, 0, iComm, ierr)
  CALL MPI_BARRIER(iComm, ierr)

  CALL MPI_TYPE_FREE(my_array_type2, ierr)

  !  boundary condition at themin and themax delta = 0 and phi = alphaVal
  DO k = 2, nzeta
     DO j = 1,npsi
        alfa(1,j,k) = alphaVal(k)
        alfa(nthe,j,k) = alphaVal(k)
     END DO
     !  extrapolate alfa to the j = 1 & npsi surfaces
     DO i = 2, nthem
        IF (iWantAlphaExtrapolation == 0) THEN
           alfa(i,npsi,k) = alfa(i,npsi-1,k) ! If the extrapolation is problematic
        ELSE
           CALL polint(radTot(i,npsi-3:npsi-1,k), alfa(i,npsi-3:npsi-1,k), radTot(i,npsi,k), alfa(i,npsi,k), dyDummy)
           !C CALL extap(alfa(i,npsi-3,k),alfa(i,npsi-2,k),alfa(i,npsi-1,k),alfa(i,npsi,k))
           !C alfa(i,npsi,k) = blendGlobal * alfa(i,npsi,k) + (1._dp - blendGlobal)*alfa(i,npsi-1,k)
           !C CALL polint(angle(i,npsi-4:npsi-1,k), alfa(i,npsi-4:npsi-1,k), angle(i,npsi,k), alfa(i,npsi,k), dyDummy)
           !C alfa(i,npsi,k) = (abs(i-nThetaEquator)/real(nthe-nThetaEquator,dp))**3 *alfa(i,npsi-1,k) + &
           !C     (1._dp - (abs(i-nThetaEquator)/real(nthe-nThetaEquator,dp))**3)*alfa(i,npsi,k)
           !C alfa(i,npsi,k) = alfa(i,npsi-1,k)
        END IF
        CALL extap(alfa(i,4,k),alfa(i,3,k),alfa(i,2,k),alfa(i,1,k)) ! This is never a problem - very close to Earth
     END DO
  END DO

  !...  set "blending" in alpha for outer iteration loop
  DO  j = 1, npsi
     DO  i = 1, nthe
        DO  k = 2, nzeta
           alfa(i,j,k) = alfa(i,j,k) * blendAlpha + alphaVal(k) * (1._dp - blendAlpha)
        END DO
        alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp * pi_d
        alfa(i,j,nzetp) = alfa(i,j,2) + 2._dp * pi_d
     END DO
  END DO

  IF (ALLOCATED(radTot)) DEALLOCATE(radTot, STAT = ierr)
  IF (ALLOCATED(angle)) DEALLOCATE(angle, STAT = ierr)


  RETURN

END SUBROUTINE iterateAlpha




