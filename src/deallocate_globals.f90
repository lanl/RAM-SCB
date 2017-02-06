SUBROUTINE deallocate_globals

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTEGER :: ierr

  IF(ALLOCATED(psi)) DEALLOCATE(psi, stat = ierr)
  IF(ALLOCATED(alfa)) DEALLOCATE(alfa, stat = ierr)
  IF(ALLOCATED(x)) DEALLOCATE(x, stat = ierr)
  IF(ALLOCATED(y)) DEALLOCATE(y, stat = ierr)
  IF(ALLOCATED(z)) DEALLOCATE(z, stat = ierr)
  IF(ALLOCATED(xx)) DEALLOCATE(xx, stat = ierr)
  IF(ALLOCATED(yy)) DEALLOCATE(yy, stat = ierr)
  IF(ALLOCATED(bsq)) DEALLOCATE(bsq, stat = ierr)
  IF(ALLOCATED(bf)) DEALLOCATE(bf, stat = ierr)
  IF(ALLOCATED(bfInitial)) DEALLOCATE(bfInitial, stat = ierr)
  IF(ALLOCATED(phij)) DEALLOCATE(phij, stat = ierr)
  IF (ALLOCATED(bj)) DEALLOCATE(bj, stat = ierr)
  IF (ALLOCATED(jParDirect)) DEALLOCATE(jParDirect, stat = ierr)

  IF (isotropy == 1) THEN 
     IF(ALLOCATED(dpdAlpha)) DEALLOCATE(dpdAlpha, stat = ierr)
     IF(ALLOCATED(dpdPsi)) DEALLOCATE(dpdPsi, stat = ierr)
  ELSE   ! if anisotropic pressure case
     IF(ALLOCATED(dPperdAlpha)) DEALLOCATE(dPperdAlpha, stat = ierr)
     IF(ALLOCATED(dPperdTheta)) DEALLOCATE(dPperdTheta, stat = ierr)
     IF(ALLOCATED(dPperdRho)) DEALLOCATE(dPperdRho, stat = ierr)
     IF(ALLOCATED(dPperdZeta)) DEALLOCATE(dPperdZeta, stat = ierr)
     IF(ALLOCATED(dPpardAlpha)) DEALLOCATE(dPpardAlpha, stat = ierr)
     IF(ALLOCATED(dPperdPsi)) DEALLOCATE(dPperdPsi, stat = ierr)
     IF(ALLOCATED(dPpardPsi)) DEALLOCATE(dPpardPsi, stat = ierr)
     IF (ALLOCATED(dBBdAlpha)) DEALLOCATE(dBBdAlpha, stat = ierr)
     IF (ALLOCATED(dBBdPsi)) DEALLOCATE(dBBdPsi, stat = ierr)
     IF(ALLOCATED(ppar)) DEALLOCATE(ppar, stat = ierr)
     IF(ALLOCATED(pper)) DEALLOCATE(pper, stat = ierr)
     IF(ALLOCATED(ratioEq)) DEALLOCATE(ratioEq, stat = ierr)
     IF(ALLOCATED(sigma)) DEALLOCATE(sigma, stat = ierr)
     IF(ALLOCATED(tau)) DEALLOCATE(tau, stat = ierr)
  END IF

  IF(ALLOCATED(dBsqdAlpha)) DEALLOCATE(dBsqdAlpha, stat = ierr)
  IF(ALLOCATED(dBsqdPsi)) DEALLOCATE(dBsqdPsi, stat = ierr)
  IF(ALLOCATED(dBsqdTheta)) DEALLOCATE(dBsqdTheta, stat = ierr)
  IF(ALLOCATED(bpar)) DEALLOCATE(bpar, stat = ierr)
  IF(ALLOCATED(bper)) DEALLOCATE(bper, stat = ierr)
  IF(ALLOCATED(xEqMidNew)) DEALLOCATE(xEqMidNew, stat = ierr)
  IF (ALLOCATED(radEqMidNew)) DEALLOCATE(radEqMidNew, stat = ierr)
  IF(ALLOCATED(pressure3D)) DEALLOCATE(pressure3D, stat = ierr)

  IF(ALLOCATED(jacobian)) DEALLOCATE(jacobian, stat = ierr)
  IF(ALLOCATED(gradRhoSq)) DEALLOCATE(gradRhoSq, stat = ierr)
  IF(ALLOCATED(gradZetaSq)) DEALLOCATE(gradZetaSq, stat = ierr)
  IF(ALLOCATED(gradRhoGradZeta)) DEALLOCATE(gradRhoGradZeta, stat = ierr)
  IF(ALLOCATED(gradRhoGradTheta)) DEALLOCATE(gradRhoGradTheta, stat = ierr)
  IF(ALLOCATED(gradThetaGradZeta)) DEALLOCATE(gradThetaGradZeta, stat = ierr)
  IF(ALLOCATED(gradPsiGradAlpha)) DEALLOCATE(gradPsiGradAlpha, stat = ierr)
  IF(ALLOCATED(curvaturePsi)) DEALLOCATE(curvaturePsi, stat = ierr)
  IF(ALLOCATED(curvatureGradPsi)) DEALLOCATE(curvatureGradPsi, stat = ierr)
  IF(ALLOCATED(gradAlphaSq)) DEALLOCATE(gradAlphaSq, stat = ierr)
  IF(ALLOCATED(radialCurvature)) DEALLOCATE(radialCurvature, stat = ierr)

  IF(ALLOCATED(g)) DEALLOCATE(g, stat = ierr)
  IF(ALLOCATED(f)) DEALLOCATE(f, stat = ierr)
  IF(ALLOCATED(fp)) DEALLOCATE(fp, stat = ierr)
  IF(ALLOCATED(fzet)) DEALLOCATE(fzet, stat = ierr)
  IF(ALLOCATED(fzetp)) DEALLOCATE(fzetp, stat = ierr)  
  IF(ALLOCATED(dela)) DEALLOCATE(dela, stat = ierr)
  IF(ALLOCATED(psival)) DEALLOCATE(psival, stat = ierr)
  IF(ALLOCATED(alphaVal)) DEALLOCATE(alphaVal, stat = ierr)
  IF(ALLOCATED(alphaValInitial)) DEALLOCATE(alphaValInitial, stat = ierr)
  IF(ALLOCATED(chiVal)) DEALLOCATE(chiVal, stat = ierr)
  IF(ALLOCATED(rhoVal)) DEALLOCATE(rhoVal, stat = ierr)
  IF(ALLOCATED(zetaVal)) DEALLOCATE(zetaVal, stat = ierr)
  IF(ALLOCATED(thetaVal)) DEALLOCATE(thetaVal, stat = ierr)
  IF(ALLOCATED(factor)) DEALLOCATE(factor, stat = ierr)

  IF(ALLOCATED(bX)) DEALLOCATE(bX, stat = ierr)
  IF(ALLOCATED(bY)) DEALLOCATE(bY, stat = ierr)
  IF(ALLOCATED(bZ)) DEALLOCATE(bZ, stat = ierr)


  RETURN

END SUBROUTINE deallocate_globals
