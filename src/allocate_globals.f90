SUBROUTINE allocate_globals

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTEGER :: ierr

  ALLOCATE(psi(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(alfa(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(x(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(y(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(z(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(xx(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(yy(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(bsq(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(bf(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(bfInitial(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(phij(nthe, npsi, nzeta+1), stat = ierr)

  ! print*, 'rank, first block done.', rank; call flush_unit(6)

  IF (isotropy == 1) THEN 
     ALLOCATE(dpdAlpha(nthe, npsi, nzeta+1), stat = ierr)
     ALLOCATE(dpdPsi(nthe, npsi, nzeta+1), stat = ierr)
     !     IF (iOuterMethod == 2) THEN ! Newton method
     ALLOCATE(dSqPdAlphaSq(nthe, npsi, nzeta+1), stat = ierr)
     ALLOCATE(dSqPdPsiSq(nthe, npsi, nzeta+1), stat = ierr)
     !     END IF
  ELSE   ! if anisotropic pressure case
     ALLOCATE(dPperdAlpha(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPperdTheta(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPperdRho(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPperdZeta(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPpardAlpha(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPperdPsi(nthe, npsi, nzeta), stat = ierr)
     ALLOCATE(dPpardPsi(nthe, npsi, nzeta), stat = ierr)
     IF (iOuterMethod == 2) THEN
        ALLOCATE(dBBdAlpha(nthe,npsi,nzeta), stat = ierr)
        ALLOCATE(dBBdPsi(nthe,npsi,nzeta), stat = ierr)
     END IF
     ALLOCATE(ppar(nthe, npsi, nzeta+1), stat = ierr)
     ALLOCATE(pper(nthe, npsi, nzeta+1), stat = ierr)
     ALLOCATE(ratioEq(npsi,nzeta), stat = ierr)
     ALLOCATE(sigma(nthe, npsi, nzeta+1), stat = ierr)
     ALLOCATE(tau(nthe, npsi, nzeta+1), stat = ierr)
  END IF

  ! print*, 'rank, second block done.', rank; call flush_unit(6)

  ALLOCATE(dBsqdAlpha(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(dBsqdPsi(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(dBsqdTheta(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(bpar(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(bper(nthe, npsi, nzeta+1), stat = ierr)
  ALLOCATE(pressure3D(nthe, npsi, nzeta+1), stat = ierr)

  ALLOCATE(jacobian(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradRhoSq(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradZetaSq(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradRhoGradZeta(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradRhoGradTheta(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradThetaGradZeta(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradPsiGradAlpha(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(curvaturePsi(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(curvatureGradPsi(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(gradAlphaSq(nthe, npsi, nzeta), stat = ierr)
  ALLOCATE(radialCurvature(nthe, npsi, nzeta), stat = ierr)

  ! print*, 'rank, third block done.', rank; call flush_unit(6)

  ALLOCATE(xEqMidNew(npsi), stat = ierr)
  ALLOCATE(radEqMidNew(npsi), stat = ierr)
  ALLOCATE(g(npsi), stat = ierr)
  ALLOCATE(f(npsi), stat = ierr)
  ALLOCATE(fp(npsi), stat = ierr)
  ALLOCATE(fzet(nzeta+1), stat = ierr)
  ALLOCATE(fzetp(nzeta+1), stat = ierr)  
  ALLOCATE(dela(nzeta+1), stat = ierr)
  ALLOCATE(psival(npsi), stat = ierr)
  ALLOCATE(alphaval(nzeta+1), stat = ierr)
  ALLOCATE(alphavalInitial(nzeta+1), stat = ierr)
  ALLOCATE(chival(nthe), stat = ierr)
  ALLOCATE(rhoVal(npsi), stat = ierr)
  ALLOCATE(zetaVal(nzeta), stat = ierr)
  ALLOCATE(thetaVal(nthe), stat = ierr)
  ALLOCATE(factor(nzeta+1), stat = ierr)

  ! cartesian components of the magnetic field now; allocate them (they are global in Module1)  
  ALLOCATE(bX(nthe,npsi,nzeta), STAT = ierr)
  ALLOCATE(bY(nthe,npsi,nzeta), STAT = ierr)
  ALLOCATE(bZ(nthe,npsi,nzeta), STAT = ierr)

  ! print*, 'rank, final block done.', rank; call flush_unit(6)

  RETURN

END SUBROUTINE allocate_globals
