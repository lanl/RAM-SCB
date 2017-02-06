!******************************************************************************
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

SUBROUTINE computeBandJacob_initial

  USE nrtype
  USE Module1
  USE mpi

  IMPLICIT NONE

  INTERFACE EZ_Spline_3D_derivatives
     SUBROUTINE Spline_coord_derivs(x_1, x_2, x_3, f_input, derivsX1, derivsX2, derivsX3)
       USE Module1
       USE EZspline_obj ! import the modules
       USE EZspline  
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP
       !       REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_2
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_3 ! independent variables
       REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f_input
       REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: derivsX1, derivsX2, derivsX3 
     END SUBROUTINE Spline_coord_derivs
  END INTERFACE

  INTEGER :: i, j, k, ierr, idealerr
  REAL(DP) :: yyp, phi, deltaPhi
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: derivXTheta, derivXRho, derivXZeta, &
       & derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       & gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq
  ! gradRhoSq, gradRhoGradZeta are global

  !  print*, 'ComputeBandJacob_initial has been called.'

 ! Get my rank:
 !C CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  ! Get the total number of processors:
 !C CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProc, ierr)
  !C PRINT *, "computeBandJacob_initial: Processor", rank, "of", numProc, ": Hello World!" 

!C  end_time = MPI_WTIME( )    ! start timer, measured in seconds
!C  print*, 'FIRST in CBJ: RANK, MPI_WTIME_IS_GLOBAL, time = ', rank, MPI_WTIME_IS_GLOBAL, end_time

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), derivXTheta, &
       derivXRho, derivXZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), derivYTheta, &
       derivYRho, derivYZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), derivZTheta, &
       derivZRho, derivZZeta) 
  ! Now I have all the point derivatives

!C end_time = MPI_WTIME( )    ! start timer, measured in seconds
!C print*, 'SECOND in CBJ: RANK, MPI_WTIME_IS_GLOBAL, time = ', rank, MPI_WTIME_IS_GLOBAL, end_time

  ! Time to build the Jacobian

  jacobian = derivXRho * (derivYZeta * derivZTheta - derivYTheta * derivZZeta) + derivXZeta * &
       & (derivYTheta * derivZRho - derivYRho * derivZTheta) + derivXTheta * &
       & (derivYRho * derivZZeta - derivYZeta * derivZRho)

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

  gradThetaSq = gradThetaX**2 + gradThetaY**2 + gradThetaZ**2
  gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ

  gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2

  ! Compute magnetic field
  DO  k = 1,nzeta
     DO  j = 1,npsi
        DO  i = 1,nthe
           bsq(i,j,k) = (gradRhoSq(i,j,k) * gradZetaSq(i,j,k) - gradRhoGradZeta(i,j,k) **2) & 
                & * (f(j) * fzet(k)) **2 
           bfInitial(i,j,k) = SQRT(bsq(i,j,k))
           bf(i,j,k) = bfInitial(i,j,k)
           IF (ABS(bsq(i,j,k)) < 1e-30_dp) THEN
              PRINT*, i, j, k, bsq(i,j,k)
              STOP 'Problem with Bsq in computeBandJacob.' 
           END IF
        END DO
     END DO
  END DO


  RETURN

END SUBROUTINE computeBandJacob_initial


