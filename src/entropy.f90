SUBROUTINE entropy(ent_local, vol_local, iteration_local)

  USE nrtype 
  USE Module1
  USE ezcdf

  IMPLICIT NONE

  INTERFACE Spline_2D_derivs
     SUBROUTINE Spline_2D_derivs(x1, x2, f, derivsX1, derivsX2)
       USE EZspline_obj ! import the modules
       USE EZspline  
       use nrtype, ONLY : DP
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP 
       REAL(r8), DIMENSION(:), INTENT(IN) :: x1, x2 ! independent variables
       REAL(r8), DIMENSION(:,:), INTENT(IN) :: f, derivsX1, derivsX2  ! function values
     END SUBROUTINE Spline_2D_derivs
  END INTERFACE


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

  real(DP) :: ent_local(:,:), vol_local(:,:)
  integer, intent(IN) :: iteration_local

  INTEGER :: i, j, k, ierr, idealerr, ncdfId
  REAL(DP) :: yyp, phi, deltaPhi
  REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: derivXTheta, derivXRho, derivXZeta, &
       & derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       & gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq, derivBsqRho, derivBsqZeta
  ! gradRhoSq, gradRhoGradZeta are global
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dVoldXEq, dVoldYEq, dVoldZeta, dVoldAlpha, dVoldRho, &
       dVoldPsi, dEntdXEq, dEntdYEq, dEntdZeta, dEntdAlpha, dEntdRho, &
       dEntdPsi, facVasGlobal, secondTermB 
  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  REAL(DP) :: delS
  REAL(DP) :: rr1, rr2, zangle, thangle, thangleOnEarth, rr, dza, dya 
  REAL(DP) ::  dipoleFactor, dipoleFactor4RE, factorIncrease

  ALLOCATE(dVoldZeta(npsi,nzeta), dVoldAlpha(npsi,nzeta), &
       dVoldRho(npsi,nzeta), dVoldPsi(npsi,nzeta), dEntdRho(npsi,nzeta), dEntdPsi(npsi,nzeta), &
       dEntdZeta(npsi,nzeta), dEntdAlpha(npsi,nzeta))

  vol_local = 0.0_dp
  if (iteration_local == 1) ent_local = 0.0_dp

  DO j = 1, npsi
     DO k = 1, nzeta
        DO i = 1, nThetaEquator-1  ! only (nThetaEquator-1) integration intervals
           vol_local(j,k) = vol_local(j,k) + jacobian(i,j,k) * 1./(f(j)*fzet(k)) * dt
        END DO
        if (iteration_local == 1) ent_local(j,k) = pressure3D(nThetaEquator,j,k) * (vol_local(j,k))**(5./3.)
     END DO
  END DO

  ! Compute grad(fluxVolume) components

  CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), vol_local(:,1:nzeta), &
       dVoldRho(:,1:nzeta), dVoldZeta(:,1:nzeta))
  CALL Spline_2D_derivs(rhoVal, zetaVal(1:nzeta), ent_local(:,1:nzeta), &
       dEntdRho(:,1:nzeta), dEntdZeta(:,1:nzeta))

  DO j = 1, npsi
     dVoldPsi(j,:) = 1._dp / f(j) * dVoldRho(j,:)
     dEntdPsi(j,:) = 1._dp / f(j) * dEntdRho(j,:)
  END DO

  DO k = 1, nzeta
     dVoldAlpha(:,k) = dVoldZeta(:,k) / fzet(k)
     dEntdAlpha(:,k) = dEntdZeta(:,k) / fzet(k)
  END DO

  ! Allocate derivXTheta etc.

  ALLOCATE(derivXTheta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivXRho(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivXZeta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivYTheta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivYRho(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivYZeta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivZTheta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivZRho(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(derivZZeta(nthe, npsi, nzeta), STAT = ierr)

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), derivXTheta, &
       derivXRho, derivXZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), derivYTheta, &
       derivYRho, derivYZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), derivZTheta, &
       derivZRho, derivZZeta) 
  ! Now I have all the point derivatives

  ! Time to build the Jacobian

  jacobian = derivXRho * (derivYZeta * derivZTheta - derivYTheta * derivZZeta) + derivXZeta * &
       & (derivYTheta * derivZRho - derivYRho * derivZTheta) + derivXTheta * &
       & (derivYRho * derivZZeta - derivYZeta * derivZRho)

  ! allocate gradRhoX, etc.

  ALLOCATE(gradRhoX(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradRhoY(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradRhoZ(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradZetaX(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradZetaY(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradZetaZ(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradThetaX(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradThetaY(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradThetaZ(nthe, npsi, nzeta), STAT = ierr)

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

  gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ

  gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2

  DEALLOCATE(derivXTheta, STAT = idealerr)
  DEALLOCATE(derivXRho, STAT = idealerr)
  DEALLOCATE(derivXZeta, STAT = idealerr)
  DEALLOCATE(derivYTheta, STAT = idealerr)
  DEALLOCATE(derivYRho, STAT = idealerr)
  DEALLOCATE(derivYZeta, STAT = idealerr)
  DEALLOCATE(derivZTheta, STAT = idealerr)
  DEALLOCATE(derivZRho, STAT = idealerr)
  DEALLOCATE(derivZZeta, STAT = idealerr)

  ALLOCATE(dVoldXEq(npsi,nzeta), stat = ierr)
  ALLOCATE(dVoldYEq(npsi,nzeta), stat = ierr) 
  ALLOCATE(dEntdXEq(npsi,nzeta), stat = ierr)
  ALLOCATE(dEntdYEq(npsi,nzeta), stat = ierr) 

  DO j = 1, npsi
     DO k = 1, nzeta
        dVoldXEq(j,k) = dVoldAlpha(j,k) * fzet(k) * gradZetaX(nThetaEquator,j,k) + &
             dVoldPsi(j,k) * f(j) * gradRhoX(nThetaEquator,j,k)
        dVoldYEq(j,k) = dVoldAlpha(j,k) * fzet(k) * gradZetaY(nThetaEquator,j,k) + &
             dVoldPsi(j,k) * f(j) * gradRhoY(nThetaEquator,j,k)
        dEntdXEq(j,k) = dEntdAlpha(j,k) * fzet(k) * gradZetaX(nThetaEquator,j,k) + &
             dEntdPsi(j,k) * f(j) * gradRhoX(nThetaEquator,j,k)
        dEntdYEq(j,k) = dEntdAlpha(j,k) * fzet(k) * gradZetaY(nThetaEquator,j,k) + &
             dEntdPsi(j,k) * f(j) * gradRhoY(nThetaEquator,j,k)
     END DO
  END DO


  DEALLOCATE(gradRhoX, STAT = ierr)
  DEALLOCATE(gradRhoY, STAT = ierr)
  DEALLOCATE(gradRhoZ, STAT = ierr)
  DEALLOCATE(gradZetaX, STAT = ierr)
  DEALLOCATE(gradZetaY, STAT = ierr)
  DEALLOCATE(gradZetaZ, STAT = ierr)
  DEALLOCATE(gradThetaX, STAT = ierr)
  DEALLOCATE(gradThetaY, STAT = ierr)
  DEALLOCATE(gradThetaZ, STAT = ierr)

  ALLOCATE(facVasGlobal(npsi,nzeta), stat = ierr)
  ALLOCATE(secondTermB(npsi,nzeta), stat = ierr)

  DO j = 1, npsi
     DO k = 1, nzeta
        facVasGlobal(j,k) = (dVoldPsi(j,k)*dPdAlpha(nThetaEquator,j,k) - &
             dVoldAlpha(j,k)*dPdPsi(nThetaEquator,j,k)) * bf(1,j,k)
        secondTermB(j,k) = jacobian(1,j,k)/bf(1,j,k)* &
             (f(j)*fzet(k)**2 * gradRhoGradTheta(1,j,k)*gradZetaSq(1,j,k)*dPdAlpha(1,j,k) + &
             f(j)**2 * fzet(k) * gradRhoGradTheta(1,j,k) * gradRhoGradZeta(1,j,k)*dPdPsi(1,j,k) - &
             f(j)*fzet(k)**2 * gradRhoGradZeta(1,j,k) * gradThetaGradZeta(1,j,k)*dPdAlpha(1,j,k) - &
             f(j)**2 * fzet(k) * gradRhoSq(1,j,k) * gradThetaGradZeta(1,j,k) * dPdPsi(1,j,k))
     END DO
  END DO

  DO k = 1, nzeta
     DO j = 1, npsi
        xx(1,j,k) = SQRT(x(1,j,k)**2 + y(1,j,k)**2)
        rr2 = x(1,j,k)**2 + y(1,j,k)**2 + z(1,j,k)**2
        rr1 = SQRT(rr2)
        ! thangle is the polar angle on the inner sphere delimiting the computational domain
        thangle = ASIN(z(1,j,k) / rr1)
        ! thangleOnEarth is the polar angle on Earth's surface
        thangleOnEarth = ACOS(SQRT((COS(thangle))**2/r0Start))

        ! Dipole field at the Earth's surface
        dipoleFactor = SQRT(1. + 3. * (SIN(thangleOnEarth))**2) 
        dipoleFactor4RE = SQRT(1. + 3. * (SIN(thangle))**2)
        factorIncrease = dipoleFactor * r0Start**3 / dipoleFactor4RE

        IF (INT(r0Start) /= 1) facVasGlobal(j,k) = facVasGlobal(j,k) * factorIncrease
     END DO
  END DO


  CALL cdf_open(ncdfId, 'entropy.cdf', 'w')

  dimlens(1) = npsi
  dimlens(2) = 1
  CALL cdf_define(ncdfId, 'psival', dimlens, 'R8')

  dimlens(1) = npsi
  dimlens(2) = nzeta
  CALL cdf_define(ncdfId, 'xEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'yEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'fluxVolume', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'dVoldXEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'dVoldYEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'entropy', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'dEntdXEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'dEntdYEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'facVasG', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'secondTermB', dimlens, 'R8')

  CALL cdf_write(ncdfId, 'psival', psival(1:npsi))
  CALL cdf_write(ncdfId, 'xEq', x(nThetaEquator,1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'yEq', y(nThetaEquator,1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'fluxVolume', vol_local(1:npsi, 1:nzeta))
  CALL cdf_write(ncdfId, 'dVoldXEq', dVoldXEq(1:npsi, 1:nzeta))
  CALL cdf_write(ncdfId, 'dVoldYEq', dVoldYEq(1:npsi, 1:nzeta))
  CALL cdf_write(ncdfId, 'entropy', ent_local(1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'dEntdXEq', dEntdXEq(1:npsi, 1:nzeta))
  CALL cdf_write(ncdfId, 'dEntdYEq', dEntdYEq(1:npsi, 1:nzeta))
  CALL cdf_write(ncdfId, 'facVasG', facVasGlobal(1:npsi, 1:nzeta)*pjconst)
  CALL cdf_write(ncdfId, 'secondTermB', secondTermB(1:npsi, 1:nzeta)*pjconst)

  CALL cdf_close(ncdfId)

  DEALLOCATE(facVasGlobal, stat = idealerr)
  DEALLOCATE(secondTermB, stat = idealerr)

  DEALLOCATE(dVoldXEq, stat = idealerr)
  DEALLOCATE(dVoldYEq, stat = idealerr)
  DEALLOCATE(dEntdXEq, stat = idealerr)
  DEALLOCATE(dEntdYEq, stat = idealerr)

  ! Can de-allocate derivXRho etc.



1000 RETURN

END SUBROUTINE entropy



