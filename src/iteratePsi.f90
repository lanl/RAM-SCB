!******************************************************************************
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

SUBROUTINE iteratePsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)
  !   iteration procedure to obtain matrix inversion

  USE nrtype
  USE Module1
  USE mpi
  USE ModRamMpi, ONLY: iComm

  IMPLICIT NONE

  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd,vec1,vec2,vec3,vec4,vec6, &
       vec7,vec8,vec9,vecr
  REAL(DP) :: om(na) 
  REAL(DP) :: omegaOpt
  REAL(DP) :: omc, anorm, anormf, resid, diff, ano, sumbtest, anormResid, anormError, &
       RESULT
  REAL(DP), ALLOCATABLE :: psiPrev(:,:,:), psiPrevTemp(:,:,:)
  INTEGER :: j, k, i, ierr, ni, ict, jz, jp, jm, iz, im, ip, nratio, myAlphaBegin, myAlphaEnd, &
       my_array_type_psi2, alphaRangeDiff, resultInt
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: alphaBeg(0:63), alphaEnd(0:63), param(0:63), my_array_type_psi(0:63)  ! Assumes no more than 64 processors
  INTEGER :: sizes(3), subsizes(3), starts(3)
  INTEGER :: request(3, 0:63)

  !     perform SOR iteration
  !..   choose for inner iteration loop

  ALLOCATE(psiPrevTemp(nthe,npsi,nzeta+1), STAT = ierr)
  ALLOCATE(psiPrev(nthe,npsi,nzeta+1), STAT = ierr)

  rjac = 1._dp - 2._dp*pi_d*pi_d / (REAL(nthe,dp)*REAL(nthe,dp) + REAL(npsi,dp)*REAL(npsi,dp)) ! Radius of conv. of Jacobi iteration,
  ! could be used to find omega optimal in SOR
  !rjac = (cos(pi_d/real(npsi,dp)) + (dr/dt)**2 * cos(pi_d/real(nthe,dp))) / &
  !     (1._dp +  (dr/dt)**2)

  om = 1._dp
  omegaOpt = 2._dp / (1._dp + SQRT(1._dp - rjac*rjac))

  ! Give each processor a chunk of the total number of alpha surfaces to compute (alfphas - k - are not connected here)

  nratio = CEILING(REAL(nzeta-1, dp)/REAL(numProc,dp))

  IF (numProc /= 1) THEN
     IF (rank == 0) THEN
        alphaBeg(rank) = 2
        alphaEnd(rank) = nratio + 1
     ELSE IF (rank /= numProc-1) THEN
        alphaBeg(rank) = nratio * rank + 2
        alphaEnd(rank) = MIN(nzeta, nratio * (rank+1) + 1)
     ELSEIF (rank == numProc-1) THEN
        alphaBeg(rank) = MIN(nzeta, nratio * (numProc-1) + 2)
        alphaEnd(rank) = nzeta
        ! print*, 'rank, alphaBeg, alphaEnd = ', rank, alphaBeg(rank), alphaEnd(rank); call flush(6)
     END IF
  ELSE
     alphaBeg(rank) = 2
     alphaEnd(rank) = nzeta
  END IF


  ! Make master know all alpha ranges

  CALL MPI_BARRIER(iComm, ierr)
  myAlphaBegin = alphaBeg(rank)
  CALL MPI_GATHER(myAlphaBegin, 1, MPI_INTEGER, param, 1, MPI_INTEGER, &
       0, iComm, ierr)
  IF (rank == 0)  alphaBeg(0:numProc-1) = param(0:numProc-1)
  myAlphaEnd = alphaEnd(rank)
  CALL MPI_GATHER(myAlphaEnd,1,MPI_INTEGER,param,1,MPI_INTEGER, &
       0, iComm, ierr)
  IF (rank == 0) alphaEnd(0:numProc-1) = param(0:numProc-1)

  psiPrev = psi
  nisave = 0

  alphaLoop: DO k = alphaBeg(rank), alphaEnd(rank)

     ni = 1

     anormf = 0._dp
     DO jz = 2, npsi-1
        DO iz = 2, nthem
           anormf = anormf + ABS(vecr(iz,jz,k))
        END DO
     END DO


     Iterations: DO WHILE (ni <= nimax)

        anormResid = 0._dp

        jLoop:       DO jz = 2, npsi-1
           jp = jz + 1
           jm = jz - 1
           iLoop:        DO  iz = 2, nthem
              im = iz - 1
              ip = iz + 1
              resid = -vecd(iz,jz,k)*psi(iz,jz,k) + vec1(iz,jz,k)*psi(im,jm,k) + & 
                   vec2(iz,jz,k)*psi(iz,jm,k) + vec3(iz,jz,k)*psi(ip,jm,k) + & 
                   vec4(iz,jz,k)*psi(im,jz,k) + vec6(iz,jz,k)*psi(ip,jz,k) + & 
                   vec7(iz,jz,k)*psi(im,jp,k) + vec8(iz,jz,k)*psi(iz,jp,k) + & 
                   vec9(iz,jz,k)*psi(ip,jp,k) - vecr(iz,jz,k)
              anormResid = anormResid + ABS(resid)
              psi(iz,jz,k) = psi(iz,jz,k) + om(k)*resid/vecd(iz,jz,k)
           END DO iLoop
        END DO jLoop

        om(k) = omegaOpt

        IF (isSorDetailNeeded == 1) THEN
           IF(MOD(ni,50) == 1 .AND. k == alphaBeg(rank)) THEN
              WRITE(*, '(A, 1X, I3, 1X, I4, A1, 1X, 3(E12.3))') 'RANK, ni, Residue for psi, RHS , residue/RHS at ni = ', RANK, &
                   ni, ':', anormResid, anormf, anormResid/anormf; CALL FLUSH(6)
           END IF
        END IF

        IF((ni > 1) .AND. (anormResid < 1e-3_dp*anormf)) THEN
           !           print*, 'exiting Iterations, ni = ', ni
           EXIT Iterations
        END IF

        ni = ni + 1

     END DO Iterations
     ! print*, 'k, ni = ', k, ni
     nisave = MAX(nisave, ni)
  END DO AlphaLoop

  sumdb = SUM(ABS(psi(2:nthem,2:npsim,alphaBeg(rank):alphaEnd(rank)) - psiprev(2:nthem,2:npsim,alphaBeg(rank):alphaEnd(rank))))
  sumb = SUM(ABS(psi(2:nthem,2:npsim,alphaBeg(rank):alphaEnd(rank))))
  diffmx = MAXVAL(ABS(psi(2:nthem,2:npsim,alphaBeg(rank):alphaEnd(rank))/&
       psiprev(2:nthem,2:npsim,alphaBeg(rank):alphaEnd(rank)) - 1)) * &
       REAL(npsi, DP)
  IF (ALLOCATED(psiPrevTemp))  DEALLOCATE(psiPrevTemp)
  IF (ALLOCATED(psiPrev)) DEALLOCATE(psiPrev)

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
 ! print*, 'nisave before reduce = ', nisave
  CALL MPI_REDUCE(nisave, resultInt, 1, MPI_INTEGER, MPI_MAX, &
       0, iComm, ierr)
 ! print*, 'resultInt = ', resultInt
  IF (rank == 0) nisave = resultInt

  ! Have to assemble psi on the master node from the calculations done by all processors

  CALL MPI_BARRIER(iComm, ierr)  ! necessary to make sure all data is computed by all processors

  sizes(1) = nthe
  sizes(2) = npsi
  sizes(3) = nzeta+1
  subsizes(1) = nthe
  subsizes(2) = npsi
  starts(1) = 0
  starts(2) = 0

  IF (rank /= 0) THEN
     subsizes(3) = alphaEnd(rank) - alphaBeg(rank)+1
     starts(3) = alphaBeg(rank) -1 ! because it starts from 0
     CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
          my_array_type_psi(rank), ierr)
     CALL MPI_TYPE_COMMIT(my_array_type_psi(rank), ierr)
     CALL MPI_SEND(psi, 1, my_array_type_psi(rank), 0, 27+3*rank, iComm, ierr)
  END IF

  IF (rank == 0 .AND. numProc > 1) THEN
     DO i = 1, numProc-1
        subsizes(3) = alphaEnd(i) - alphaBeg(i) + 1
        starts(3) = alphaBeg(i)-1
        CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
             my_array_type_psi(i), ierr)
        CALL MPI_TYPE_COMMIT(my_array_type_psi(i), ierr)
        CALL MPI_IRECV(psi, 1, my_array_type_psi(i), i, 27+3*i, iComm, request(1, i), ierr)
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
        CALL MPI_TYPE_FREE(my_array_type_psi(i), ierr)
     END DO
  ELSE
     CALL MPI_TYPE_FREE(my_array_type_psi(rank), ierr)
  END IF

  ! Broadcast the new psi array to all processors

  subsizes(3) = nzeta+1
  starts(3) = 0

  CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
       my_array_type_psi2, ierr)
  CALL MPI_TYPE_COMMIT(my_array_type_psi2, ierr)
  CALL MPI_BCAST(psi, 1, my_array_type_psi2, 0, iComm, ierr)
  CALL MPI_BARRIER(iComm, ierr)
  CALL MPI_TYPE_FREE(my_array_type_psi2, ierr)

 ! print*, 'nisave, ni, blendPsi = ', nisave, ni, blendPsi
  !...  set blend for outer iteration loop, and apply periodic boundary conditions
  DO  j = 1, npsi
     DO  i = 1, nthe
        DO  k = 2,nzeta
           psi(i,j,k) = psi(i,j,k) * blendPsi + psival(j) * (1._dp - blendPsi)
        END DO
        psi(i,j,1) = psi(i,j,nzeta)
        psi(i,j,nzetp) = psi(i,j,2)
     END DO
  END DO



  RETURN

END SUBROUTINE iteratePsi
