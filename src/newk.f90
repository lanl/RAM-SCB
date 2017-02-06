SUBROUTINE newk(vec_x, vec_dk)
  !     define right hand side of alpha equation  

  USE nrtype
  USE Module1
  IMPLICIT NONE

  INTEGER         :: i, j, k 
  REAL(DP), DIMENSION(:, :, :), INTENT(OUT) :: vec_x
  REAL(DP), DIMENSION(:, :, :), INTENT(IN OUT) ::  vec_dk
!  LOGICAL :: isnand

  IF (isotropy == 1) THEN      ! Isotropic case
     DO j = 2, npsi-1
        IF (iOuterMethod == 2) THEN   ! Newton method
           DO k = 2, nzeta
              vec_x(1:nthe, j, k) = - jacobian(1:nthe,j,k)/f(j)**2 * (dpdAlpha(1:nthe,j,k) - dSqpdAlphaSq(1:nthe,j,k) * &
                   alphaVal(k))
              vec_dk(1:nthe, j, k) = vec_dk(1:nthe, j, k) - jacobian(1:nthe,j,k)/f(j)**2 * dSqpdAlphaSq(1:nthe,j,k)
           END DO
        ELSE                          ! Picard method
           vec_x(1:nthe, j, 1:nzeta) = - dpdAlpha(1:nthe, j, 1:nzeta) * &
                & jacobian(1:nthe, j, 1:nzeta) * 1._dp / f(j)**2 
        END IF
     END DO

  ELSE                    ! Anisotropic case
     DO i = 2, nthe-1
        DO j = 2, npsi-1
           DO k = 2, nzeta
              IF (iOuterMethod == 1) THEN ! Picard method
                 vec_x(i,j,k) = jacobian(i,j,k) / f(j)**2 * (-1./sigma(i,j,k) * dPperdAlpha(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j)**2 * fzet(k) * (gradRhoSq(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradRhoGradZeta(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k))*0.5*dBsqdTheta(i,j,k)) - &
                      (1. - sigma(i,j,k)) / sigma(i,j,k) * 0.5 * dBsqdAlpha(i,j,k))
               !           print*, 'newk: i, j, k, f(j), jac, dpdAlpha, gradRhoSq, gradRhoGradZeta, sigma, bsq, vec_x: ', &
               !                i, j, k, real(f(j),sp), real(jacobian(i,j,k),sp), real(dPperdAlpha(i,j,k),sp), &
               !                real(gradRhoSq(i,j,k),sp), real(gradRhoGradZeta(i,j,k),sp), &
               !                real(sigma(i,j,k),sp), real(bsq(i,j,k),sp), real(vec_x(i,j,k),sp)
 !                IF (isnand(vec_x(i,j,k))) STOP 'newk NaN.'
              ELSE ! Newton method
                 vec_x(i,j,k) = jacobian(i,j,k) / f(j)**2 * (-1./sigma(i,j,k) * dPperdAlpha(i,j,k) &
                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j)**2 * fzet(k) * (gradRhoSq(i,j,k)* &
                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradRhoGradZeta(i,j,k)) * &
                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k))*0.5*dBsqdTheta(i,j,k)) - &
                      (1. - sigma(i,j,k)) / sigma(i,j,k) * 0.5 * dBsqdAlpha(i,j,k)) 
                 vec_x(i,j,k) = vec_x(i,j,k) - jacobian(i,j,k)/((f(j))**2) * dBBdAlpha(i,j,k) *alphaVal(k)
                 vec_dk(i,j,k) = vec_dk(i,j,k) + jacobian(i,j,k)/(f(j)**2) * dBBdAlpha(i,j,k)
              END IF
           END DO
        END DO
     END DO

  END IF


  RETURN
END SUBROUTINE newk
