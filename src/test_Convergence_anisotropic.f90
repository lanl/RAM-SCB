
SUBROUTINE test_Convergence_anisotropic

  USE nrtype
  USE Module1
  USE ModIoUnit, ONLY: UNITTMP_
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

  INTEGER :: i, j, k, id, ierr, idealerr

  REAL(DP) :: normDiffRel, volume
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: derivXTheta, derivXRho, derivXZeta, &
       & derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       & gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq, &
       derivBsqTheta, derivBsqRho, derivBsqZeta, derivNU1, derivNU2
  ! gradRhoSq, gradRhoGradZeta are global

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jGradRhoPartialTheta, derivjGradRhoPartialTheta, &
       & jGradRhoPartialZeta, derivjGradRhoPartialZeta, jGradRho, jGradRhoFactor, jGradZetaPartialRho, &
       & derivjGradZetaPartialRho, jGradZetaPartialTheta, derivjGradZetaPartialTheta, jGradZeta, &
       & jGradZetaFactor, jGradThetaPartialRho, derivjGradThetaPartialRho, jGradThetaPartialZeta, &
       & derivjGradThetaPartialZeta, jGradTheta, jGradThetaFactor, phiToroid, derivPhiRho, derivPhiZeta, &
       & derivPhiTheta, derivDiffPTheta

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jCrossBUpRho, jCrossBUpZeta, jCrossBUpTheta, &
       derivjCrossBUpRho, derivjCrossBUpZeta, derivjCrossBUpTheta, jCrossBMinusGradPSq, &
       jCrossBMinusGradPMod, jCrossBSq, jCrossB, gradPSq, gradP
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jrrInt, jrr, jzzInt, jzz, jrtInt, jrt, jztInt, jzt, &
       rhoCompSq, zetaCompSq, thetaCompSq, curlJCrossBSq, curlJCrossB
 
!  LOGICAL, EXTERNAL :: isnand ! Intrinsic for Portland Group Fortran

  !**********************************************************************************************************!


  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), derivXTheta, &
       & derivXRho, derivXZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), derivYTheta, &
       & derivYRho, derivYZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), derivZTheta, &
       & derivZRho, derivZZeta) 
  ! Now I have all the point derivatives

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

  ! B field squared is obtained now, then the B field bf; in the following, i, j, k can be 
  ! taken over the full domain because we have all the required quantities everywhere;
  ! thus, extrapolation for Bfield is not necessary anymore

  DO  k = 1,nzeta
     DO  j = 1,npsi
        DO  i = 1,nthe
           bsq(i,j,k) = (gradRhoSq(i,j,k) * gradZetaSq(i,j,k) - gradRhoGradZeta(i,j,k) **2) & 
                & * (f(j) * fzet(k)) **2 
           bf(i,j,k) = SQRT(bsq(i,j,k))
        END DO
     END DO
  END DO

  ! j dot gradRho

  DO j = 1, npsi
     DO k = 1, nzeta
        jGradRhoPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * &
             gradRhoGradZeta(:,j,k))
        jGradRhoPartialZeta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradZetaSq(:,j,k) - gradRhoGradZeta(:,j,k) **2)
     END DO
  END DO

  !  jGradRhoPartialZeta(:,:,nzetp) = jGradRhoPartialZeta(:,:,2) 
  !  DO i = 2, nthe-1
  !     derivjGradRhoPartialTheta(i,:,:) = 0.5 * rdt * (jGradRhoPartialTheta(i+1,:,:) - jGradRhoPartialTheta(i-1,:,:))
  !  END DO
  !  DO k = 2, nzeta
  !     derivjGradRhoPartialZeta(:,:,k) = 0.5 * rdp * (jGradRhoPartialZeta(:,:,k+1) - jGradRhoPartialZeta(:,:,k-1))
  !  END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialTheta, &
       & derivjGradRhoPartialTheta, derivNU1, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialZeta, &
       & derivNU1, derivNU2, derivjGradRhoPartialZeta)

  jGradRho = 1._dp / jacobian * (derivjGradRhoPartialTheta + derivjGradRhoPartialZeta)

  ! j dot gradZeta

  DO j = 1, npsi
     DO k = 1, nzeta
        jGradZetaPartialRho(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k) **2 - gradRhoSq(:,j,k) * &
             gradZetaSq(:,j,k))
        jGradZetaPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * &
             & gradZetaSq(:,j,k))
     END DO
  END DO

  !  DO j = 2, npsi-1
  !     derivjGradZetaPartialRho(:,j,:) = 0.5 * rdr * (jGradZetaPartialRho(:,j+1,:) - jGradZetaPartialRho(:,j-1,:))
  !  END DO
  !  DO i = 2, nthe-1
  !    derivjGradZetaPartialTheta(i,:,:) = 0.5 * rdt * (jGradZetaPartialTheta(i+1,:,:) - jGradZetaPartialTheta(i-1,:,:))
  ! END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialRho, &
       & derivNU1, derivjGradZetaPartialRho, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialTheta, &
       & derivjGradZetaPartialTheta, derivNU1, derivNU2)

  ! The one below is for calculating \partial Pper/\partial \theta and \partial(J(Pper-Ppar))/\partial \theta for the anisotropic case  
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jacobian(1:nthe,1:npsi,1:nzeta) * &
       (pper(1:nthe,1:npsi,1:nzeta)-ppar(1:nthe,1:npsi,1:nzeta)), derivDiffPTheta, derivNU1, derivNU2)

  jGradZeta = 1._dp / jacobian * (derivjGradZetaPartialRho + derivjGradZetaPartialTheta)

  !***************************************************************************************************

  ! Now we have jGradRho, jGradZeta

  ! Time to compute |j x B - div dot P|

     DO k = 1, nzeta
!        print*, 'tca: k = ', k
        DO j = 1, npsi
!           if (k == 101) print*, 'tca: j, sum(jac) = ', j, sum(jacobian(:,j,k))
        ! First one is the same as in the isotropic case
        jCrossBSq(:,j,k) = f(j)**2 * fzet(k)**2 * (gradRhoSq(:,j,k) * jGradZeta(:,j,k)**2 + &
             gradZetaSq(:,j,k) * jGradRho(:,j,k)**2 - 2._dp * &
             jGradZeta(:,j,k) * jGradRho(:,j,k) * gradRhoGradZeta(:,j,k))

    !    DO i = 1, nthe
           ! print*, 'test_Convergence_anisotropic: i, j, k = ', i, j, k, REAL(jCrossBSq(i,j,k),SP)
    !       IF (isnand(jCrossBSq(i,j,k))) PRINT*, 'NaN: test_Convergence_anisotropic: i, j, k = ', i, j, k
    !    END DO

        ! Second one is different from isotropic, but keep same notation; gradPSq equiv. (div P)**2
        gradPSq(:,j,k) = gradRhoSq(:,j,k)*dPperdRho(:,j,k)**2 + gradZetaSq(:,j,k)*dPperdZeta(:,j,k)**2 + &
             gradThetaSq(:,j,k)*dPperdTheta(:,j,k)**2 + 2.*dPperdRho(:,j,k)*dPperdZeta(:,j,k)*gradRhoGradZeta(:,j,k) + &
             2.*dPperdRho(:,j,k)*dPperdTheta(:,j,k)*gradRhoGradTheta(:,j,k) + &
             2.*dPperdZeta(:,j,k)*dPperdTheta(:,j,k)*gradThetaGradzeta(:,j,k) + (derivDiffPTheta(:,j,k)/jacobian(:,j,k))**2 - &
             2.*dPperdTheta(:,j,k) * derivDiffPTheta(:,j,k)/jacobian(:,j,k)

        ! Isotropic case formula
        ! fzet(k)**2 * dpdAlpha(1:nthe,j,k)**2 * gradZetaSq(:,j,k) + f(j)**2 * &
        !             dpdPsi(1:nthe,j,k)**2 * gradRhoSq(:,j,k) + 2._dp*dpdAlpha(1:nthe,j,k)*dpdPsi(1:nthe,j,k) * &
        !             f(j) * fzet(k) * gradRhoGradZeta(:,j,k)

        ! Similar here, keep the same name for the anisotropic case

        jCrossBMinusGradPSq(:,j,k) = gradRhoSq(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k))**2 + &
             gradZetaSq(:,j,k)*(f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k))**2 + gradThetaSq(:,j,k)*dPperdTheta(:,j,k)**2 + &
             (derivDiffPTheta(:,j,k)/jacobian(:,j,k))**2 - &
             2.*gradRhoGradZeta(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k)) * &
             (f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k)) - &
             2.*gradRhoGradTheta(:,j,k)*(f(j)*fzet(k)*jGradZeta(:,j,k)-dPperdRho(:,j,k))*dPperdTheta(:,j,k) + &
             2.*gradThetaGradzeta(:,j,k)*(f(j)*fzet(k)*jGradRho(:,j,k)+dPperdZeta(:,j,k))*dPperdTheta(:,j,k) - &
             2.*dPperdTheta(:,j,k)*derivDiffPTheta(:,j,k)/jacobian(:,j,k)

     !   DO i = 1, nthe
     !      ! print*, 'test_Convergence_anisotropic: i, j, k = ', i, j, k, REAL(jCrossBSq(i,j,k),SP)
     !      IF (isnand(jCrossBMinusGradPSq(i,j,k))) PRINT*, 'NaN2: test_Convergence_anisotropic: i, j, k = ', i, j, k
     !   END DO

        ! Isotropic case formula
        !        jCrossBMinusGradPSq(:,j,k) = f(j)**2 * fzet(k)**2 * (gradRhoSq(:,j,k) * (jGradZeta(:,j,k) - &
        !             & 1. / fzet(k) * dpdPsi(1:nthe,j,k))**2 + &
        !             gradZetaSq(:,j,k) * (jGradRho(:,j,k) + 1./f(j) * dpdAlpha(1:nthe,j,k))**2 - 2._dp * &
        !             (jGradZeta(:,j,k) - 1./fzet(k) * dpdPsi(1:nthe,j,k)) * (jGradRho(:,j,k) + 1./f(j)* &
        !             dpdAlpha(1:nthe,j,k)) * gradRhoGradZeta(:,j,k))

     END DO
  END DO

!C  print*, 'tca: min, max jCrossBSq = ', MINVAL(jCrossBSq), maxval(jCrossBSq)
!C  print*, 'tca: min, max jCrossBMinusGradPSq = ', MINVAL(jCrossBMinusGradPSq), maxval(jCrossBMinusGradPSq)
!C  print*, 'tca: min, max gradPSq = ', MINVAL(gradPSq), maxval(gradPSq)
  jCrossB = SQRT(jCrossBSq)
  jCrossBMinusGradPMod = SQRT(ABS(jCrossBMinusGradPSq))
  gradP = SQRT(ABS(gradPSq))

 !C PRINT*, 'tca: sum(jCrossB) = ', SUM(jCrossB)
 !C PRINT*, 'tca: sum(jCrossBMinusGradPMod) = ', SUM(jCrossBMinusGradPMod)
 !C PRINT*, 'tca: sum(gradP) = ', SUM(gradP)

  ! Force balance quantities
  OPEN(UNITTMP_, file = TRIM(ADJUSTL(prefixOut))//'Force_balance_equatorial', status = 'replace')  ! for Matlab plotting
  WRITE(UNITTMP_, '(3I7)') nthe, npsi, nzeta
  DO j = 1, npsi
     DO k = 1, nzeta
        WRITE(UNITTMP_,'(6G14.5)') x(nThetaEquator, j, k), y(nThetaEquator, j, k), &
             jCrossB(nThetaEquator,j,k), jCrossBMinusGradPMod(nThetaEquator,j,k), &
             gradP(nThetaEquator,j,k), jacobian(nThetaEquator,j,k)

     END DO
  END DO
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file =TRIM(ADJUSTL(prefixOut))//'Force_balance_midnight', status = 'replace') 
  WRITE(UNITTMP_, '(3I7)') nthe, npsi, nzeta
  DO j = 1, npsi
     DO i = 1, nthe
        WRITE(UNITTMP_,'(6G14.5)') x(i,j,nZetaMidnight), z(i,j,nZetaMidnight), &
             jCrossB(i,j,nZetaMidnight), &
             jCrossBMinusGradPMod(i,j,nZetaMidnight), gradP(i,j,nZetaMidnight), &
             jacobian(i,j,nZetaMidnight)
     END DO
  END DO
  CLOSE(UNITTMP_)

  !  OPEN(102, file = 'Force_balance_noon', status = 'replace') 
  !  WRITE(102, '(3I7)') nthe, npsi, nzeta
  !  DO j = 1, npsi
  !     DO i = 1, nthe
  !        WRITE(102,'(6G14.5)') x(i,j,2), z(i,j,2), &
  !             jCrossB(i,j,2), jCrossBMinusGradPMod(i,j,2), &
  !             gradP(i,j,2), jacobian(i,j,2)
  !     END DO
  !  END DO
  !  CLOSE(102)

  !  OPEN(20, file = 'jxBconv', action = 'write', status = 'replace')  
  !  DO i = 1, nthe
  !     DO j = 1, npsi
  !        WRITE(20, *) x(i,j, nZetaMidnight), z(i, j, nZetaMidnight), &
  !             & jacobian(i,j,nZetaMidnight)*jCrossB(i,j,nZetaMidnight)
  !     END DO
  !  END DO
  !  CLOSE(20)

  normDiff = 0.0_dp
  normDiffRel = 0.0_dp
  normJxB  = 0.0_dp
  normGradP = 0.0_dp
  volume   = 0.0_dp

  DO i = 2, nthe-1
     DO j = 2, npsi-1
        DO k = 2, nzeta
           IF (2.*pper(i,j,k) > 1.E-1_dp*bsq(i,j,k)) THEN  
              ! Only in regions with beta > 1E-2 
              ! (in regions of low plasma beta, the pressure does not change the magnetic field)
           normDiff = normDiff + jacobian(i,j,k) * dr * dpPrime * dt * &
                jCrossBMinusGradPMod(i,j,k)
           normDiffRel = normDiffRel + jacobian(i,j,k) * dr * dpPrime * dt * &
                jCrossBMinusGradPMod(i,j,k) / jCrossB(i,j,k) 
           normJxB = normJxB + jacobian(i,j,k) * dr * dpPrime * dt * &
                jCrossB(i,j,k)
           normGradP = normGradP + jacobian(i,j,k) * dr * dpPrime * dt * gradP(i,j,k)
           volume = volume + jacobian(i,j,k) * dr * dpPrime * dt
        END IF
        END DO
     END DO
  END DO

  ! Normalize to total computational volume
  normDiff = normDiff/volume
  normJxB = normJxB/volume
  normGradP = normGradP/volume

  !  Norms of |jxB-grad P|,      |jxB|,      |gradP| 
  WRITE(*, '(3 (3X, E14.5))') normDiff, normJxB, normGradP
  WRITE(iUnitLog, '(3 (3X, E14.5))') normDiff, normJxB, normGradP


  RETURN

END SUBROUTINE test_Convergence_anisotropic









