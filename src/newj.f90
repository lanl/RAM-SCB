SUBROUTINE newj(vec_r, vec_dj)
  !     define right hand side of J dot grad(alfa) equilibrium equation

  USE Module1
  USE nrtype

  IMPLICIT NONE

  INTEGER         :: i, j, k 
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: vec_r
  REAL(DP), DIMENSION(:,:,:), INTENT(IN OUT) :: vec_dj


  Isotropic_case: IF (isotropy == 1) THEN 
     DO k = 1, nzeta
        Newton_method:   IF (iOuterMethod == 2) THEN
           DO j = 1, npsi
              vec_r(1:nthe, j, k) = jacobian(1:nthe, j, k)/fzet(k)**2 * (dpdPsi(1:nthe, j, k) - &
                   dSqPdPsiSq(1:nthe,j,k)*psival(j))
              vec_dj(1:nthe, j, k) = vec_dj(1:nthe, j, k) + jacobian(1:nthe,j,k)/fzet(k)**2 *dSqPdPsiSq(1:nthe,j,k)
           END DO
        ELSE            ! Picard method
           vec_r(1:nthe, 1:npsi, k) = jacobian(1:nthe, 1:npsi, k) * dpdPsi(1:nthe, 1:npsi, k) / &
                fzet(k)**2
        END IF Newton_method
     END DO
  ELSE       ! Anisotropic case
     DO k = 1, nzeta
        DO j = 1, npsi
           DO i = 1, nthe 
              IF (iOuterMethod == 1) THEN ! Picard method
                 vec_r(i,j,k) = jacobian(i,j,k) / fzet(k)**2 * (1./sigma(i,j,k) * dPperdPsi(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j) * fzet(k)**2 * (gradRhoGradZeta(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradZetaSq(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k)) * 0.5_dp * dBsqdTheta(i,j,k)) + &
                      (1.-sigma(i,j,k)) / sigma(i,j,k) * 0.5_dp * dBsqdPsi(i,j,k))
              ELSE ! Newton method
                 vec_r(i,j,k) = jacobian(i,j,k) / fzet(k)**2 * (1./sigma(i,j,k) * dPperdPsi(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j) * fzet(k)**2 * (gradRhoGradZeta(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradZetaSq(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k)) * 0.5_dp * dBsqdTheta(i,j,k)) + &
                      (1.-sigma(i,j,k)) / sigma(i,j,k) * 0.5_dp * dBsqdPsi(i,j,k))
                 vec_r(i,j,k) = vec_r(i,j,k) - jacobian(i,j,k)/(fzet(k)**2) * dBBdPsi(i,j,k)*psival(j)
                 vec_dj(i,j,k) = vec_dj(i,j,k) + jacobian(i,j,k)/(fzet(k)**2) * dBBdPsi(i,j,k)
              END IF
           END DO
        END DO
     END DO

  END IF Isotropic_case


  RETURN

END SUBROUTINE newj
