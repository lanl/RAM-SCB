SUBROUTINE Write_pressure_isotropic

  USE nrtype
  USE Module1
  USE ezcdf

  IMPLICIT NONE

  INTERFACE EZ_Spline_3D_derivatives
     SUBROUTINE Spline_coord_derivs(x_1, x_2, x_3, f_input, derivsX1, derivsX2, derivsX3)
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind; same as DP
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_2
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_3 ! independent variables
       REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f_input
       REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: derivsX1, derivsX2, derivsX3 
     END SUBROUTINE Spline_coord_derivs
  END INTERFACE

  !  INCLUDE 'netcdf.inc'

  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  INTEGER :: ncdfId

  INTEGER :: i, j, k, ierr, idealerr
  REAL(DP) :: yyp, phi, deltaPhi
  REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: derivXTheta, derivXRho, derivXZeta, &
       & derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       & gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq, derivBsqRho, derivBsqZeta
  ! gradRhoSq, gradRhoGradZeta are global
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dPdXEq, dPdYEq, dBsqdXEq, dBsqdYEq, &
       dPsidXEq, dPsidYEq


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

  ALLOCATE(dPdXEq(npsi,nzeta), stat = ierr)
  ALLOCATE(dPdYEq(npsi,nzeta), stat = ierr) 
  ALLOCATE(dPsidXEq(npsi,nzeta), stat = ierr)
  ALLOCATE(dPsidYEq(npsi,nzeta), stat = ierr) 

  DO j = 1, npsi
     DO k = 1, nzeta
        dPdXEq(j,k) = dPdAlpha(nThetaEquator,j,k) * fzet(k) * gradZetaX(nThetaEquator,j,k) + &
             dPdPsi(nThetaEquator,j,k) * f(j) * gradRhoX(nThetaEquator,j,k)
        dPdYEq(j,k) = dPdAlpha(nThetaEquator,j,k) * fzet(k) * gradZetaY(nThetaEquator,j,k) + &
             dPdPsi(nThetaEquator,j,k) * f(j) * gradRhoY(nThetaEquator,j,k)
        dPsidXEq(j,k) = gradRhoX(nThetaEquator,j,k)*f(j)
        dPsidYEq(j,k) = gradRhoY(nThetaEquator,j,k)*f(j)
     END DO
  END DO

  ! Allocate gradRhoSq. etc.

  !  ALLOCATE(gradRhoSq(nthe, npsi, nzeta), STAT = ierr)
  !  ALLOCATE(gradRhoGradZeta(nthe, npsi, nzeta), STAT = ierr)
  ! ALLOCATE(gradRhoGradTheta(nthe, npsi, nzeta), STAT = ierr)
  ALLOCATE(gradThetaSq(nthe, npsi, nzeta), STAT = ierr)
  ! ALLOCATE(gradThetaGradZeta(nthe, npsi, nzeta), STAT = ierr)
  ! ALLOCATE(gradZetaSq(nthe, npsi, nzeta), STAT = ierr)


  gradRhoSq = gradRhoX**2 + gradRhoY**2 + gradRhoZ**2
  gradRhoGradZeta = gradRhoX * gradZetaX + gradRhoY * gradZetaY + gradRhoZ * gradZetaZ
  gradRhoGradTheta = gradRhoX * gradThetaX + gradRhoY * gradThetaY + gradRhoZ * gradThetaZ

  gradThetaSq = gradThetaX**2 + gradThetaY**2 + gradThetaZ**2
  gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ

  gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2

  DO  k = 1,nzeta
     DO  j = 1,npsi
        DO  i = 1,nthe
           bsq(i,j,k) = (gradRhoSq(i,j,k) * gradZetaSq(i,j,k) - gradRhoGradZeta(i,j,k) **2) & 
                & * (f(j) * fzet(k)) **2 
           bf(i,j,k) = SQRT(bsq(i,j,k))
        END DO
     END DO
  END DO

  ALLOCATE(derivBsqRho(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(derivBsqZeta(nthe, npsi, nzeta), stat = ierr)

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, bsq(1:nthe, 1:npsi, 1:nzeta), &
       & dBsqdTheta, derivBsqRho, derivBsqZeta)

  ALLOCATE(dBSqdXEq(npsi,nzeta), stat = ierr)
  ALLOCATE(dBsqdYEq(npsi,nzeta), stat = ierr) 

  DO j = 1, npsi
     DO k = 1, nzeta
        dBSqdXEq(j,k) = derivBsqZeta(nThetaEquator,j,k) * gradZetaX(nThetaEquator,j,k) + &
             derivBsqRho(nThetaEquator,j,k) * gradRhoX(nThetaEquator,j,k)
        dBSqdYEq(j,k) = derivBsqZeta(nThetaEquator,j,k) * gradZetaY(nThetaEquator,j,k) + &
             derivBsqRho(nThetaEquator,j,k) * gradRhoY(nThetaEquator,j,k)
     END DO
  END DO

  DEALLOCATE(derivBsqRho, stat = idealerr)
  DEALLOCATE(derivBsqZeta, stat = idealerr)

  ! Opens file
  CALL     cdf_open(ncdfId, TRIM(ADJUSTL(prefixOut))//'pressure_isotropic.cdf', 'w')

  dimlens(1) = npsi
  dimlens(2) = 1
  CALL cdf_define(ncdfId, 'psival', dimlens, 'R4')

  dimlens(1) = npsi
  dimlens(2) = nzeta
  CALL      cdf_define(ncdfId, 'xEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'yEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'presEq', dimlens, 'R4')
  CALL      cdf_setatt(ncdfId, 'presEq', 'Pressure in the equatorial plane', 'nPa')
  CALL cdf_define(ncdfId, 'xIono', dimlens, 'R4')
  CALL cdf_define(ncdfId, 'yIono', dimlens, 'R4')
  CALL cdf_define(ncdfId, 'presIono', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'betaEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdAlphaEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdPsiEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdXEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdYEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dBSqdXEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dBSqdYEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPsidXEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPsidYEq', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'BOverBdipoleEq', dimlens, 'R4')


  dimlens(1) = nthe
  dimlens(2) = npsi
  CALL      cdf_define(ncdfId, 'xMid', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'zMid', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'presMid', dimlens, 'R4')
  CALL      cdf_setatt(ncdfId, 'presMid', 'Pressure in midnight merid. plane', 'nPa')
  CALL      cdf_define(ncdfId, 'betaMid', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdPsiMid', dimlens, 'R4')
  !     CALL      cdf_define(ncdfId, 'BOverBdipoleMid', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'presTotalMid', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'xNoo', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'zNoo', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'presNoo', dimlens, 'R4')
  CALL      cdf_setatt(ncdfId, 'presNoo', 'Pressure in noon merid. plane', 'nPa')
  CALL      cdf_define(ncdfId, 'betaNoo', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'dPdPsiNoo', dimlens, 'R4')
  !     CALL      cdf_define(ncdfId, 'BOverBdipoleNoo', dimlens, 'R4')
  CALL      cdf_define(ncdfId, 'presTotalNoo', dimlens, 'R4')

  CALL      cdf_write(ncdfId, 'psival', REAL(psival(1:npsi),sp))
  CALL      cdf_write(ncdfId, 'xEq', REAL(x(nThetaEquator,:,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'yEq', REAL(y(nThetaEquator,:,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'presEq', REAL(pressure3D(nThetaEquator,:,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'xIono', REAL(x(1,:,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'yIono', REAL(y(1,:,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'presIono', REAL(pressure3D(1,:,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'betaEq', REAL(2.*pressure3D(nThetaEquator,:,1:nzeta)/bsq(nThetaEquator,1:npsi,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'dPdAlphaEq', REAL(dpdAlpha(nThetaEquator,1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'dPdPsiEq', REAL(dpdPsi(nThetaEquator,1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'dPdXEq', REAL(dPdXEq(1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'dPdYEq', REAL(dPdYEq(1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'dPsidXEq', REAL(dPsidXEq(1:npsi,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'dPsidYEq', REAL(dPsidYEq(1:npsi,1:nzeta),sp))
  CALL      cdf_write(ncdfId, 'dBSqdXEq', REAL(dBSqdXEq(1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'dBSqdYEq', REAL(dBSqdYEq(1:npsi,1:nzeta)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'BOverBdipoleEq', REAL(bf(nThetaEquator,1:npsi,1:nzeta),sp) / & 
       ABS(REAL((xzero/xx(nThetaEquator,1:npsi,1:nzeta))**3,sp)))
  CALL      cdf_write(ncdfId, 'xMid', REAL(x(:,:,nZetaMidnight),sp))
  CALL      cdf_write(ncdfId, 'zMid', REAL(z(:,:,nZetaMidnight),sp))
  CALL      cdf_write(ncdfId, 'presMid', REAL(pressure3D(:,:,nZetaMidnight)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'betaMid', REAL(2.*pressure3D(:,:,nZetaMidnight)/bsq(:,:,nZetaMidnight),sp))
  CALL      cdf_write(ncdfId, 'dPdPsiMid', REAL(dpdPsi(:,:,nZetaMidnight)*pnormal,sp))
  !     CALL      cdf_write(ncdfId, 'BOverBdipoleMid', bf(:,:,nZetaMidnight) / (dipolefactorMid(:,:)* & 
  !     ABS((xzero/xx(nThetaEquator,:,nZetaMidnight))**3)))
  CALL      cdf_write(ncdfId, 'presTotalMid', REAL(pressure3D(:,:,nZetaMidnight) + 0.5_dp*bsq(:,:,nZetaMidnight),sp))
  CALL      cdf_write(ncdfId, 'xNoo', REAL(x(:,:,2),sp))
  CALL      cdf_write(ncdfId, 'zNoo', REAL(z(:,:,2),sp))
  CALL      cdf_write(ncdfId, 'presNoo', REAL(pressure3D(:,:,2)*pnormal,sp))
  CALL      cdf_write(ncdfId, 'betaNoo', REAL(2.*pressure3D(:,:,2)/bsq(:,:,2),sp))
  CALL      cdf_write(ncdfId, 'dPdPsiNoo', REAL(dpdPsi(:,:,2)*pnormal,sp))
  !     CALL      cdf_write(ncdfId, 'BOverBdipoleMid', bf(:,:,2) / (dipolefactorMid(:,:)* & 
  !     ABS((xzero/xx(nThetaEquator,:,2))**3)))
  CALL      cdf_write(ncdfId, 'presTotalNoo', REAL(pressure3D(:,:,2) + 0.5_dp*bsq(:,:,2),sp))

  CALL      cdf_close(ncdfId)

  DEALLOCATE(dPdXEq, stat = idealerr)
  DEALLOCATE(dPdYEq, stat = idealerr)
  DEALLOCATE(dPsidXEq, stat = idealerr)
  DEALLOCATE(dPsidYEq, stat = idealerr)
  DEALLOCATE(dBSqdXEq, stat = idealerr)
  DEALLOCATE(dBSqdYEq, stat = idealerr)


  ! Can de-allocate derivXRho etc.

  DEALLOCATE(derivXTheta, STAT = idealerr)
  DEALLOCATE(derivXRho, STAT = idealerr)
  DEALLOCATE(derivXZeta, STAT = idealerr)
  DEALLOCATE(derivYTheta, STAT = idealerr)
  DEALLOCATE(derivYRho, STAT = idealerr)
  DEALLOCATE(derivYZeta, STAT = idealerr)
  DEALLOCATE(derivZTheta, STAT = idealerr)
  DEALLOCATE(derivZRho, STAT = idealerr)
  DEALLOCATE(derivZZeta, STAT = idealerr)
  DEALLOCATE(gradRhoX, STAT = ierr)
  DEALLOCATE(gradRhoY, STAT = ierr)
  DEALLOCATE(gradRhoZ, STAT = ierr)
  DEALLOCATE(gradZetaX, STAT = ierr)
  DEALLOCATE(gradZetaY, STAT = ierr)
  DEALLOCATE(gradZetaZ, STAT = ierr)
  DEALLOCATE(gradThetaX, STAT = ierr)
  DEALLOCATE(gradThetaY, STAT = ierr)
  DEALLOCATE(gradThetaZ, STAT = ierr)

  RETURN

END SUBROUTINE Write_pressure_isotropic
