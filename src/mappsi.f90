
SUBROUTINE mapPsi(iSmoothMove)
  ! new psiMap, without computing linear distance
  ! first and last flux surface remain unchanged

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTERFACE splint_interface
     FUNCTION splint(xa,ya,y2a,x)
       USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
       USE nr, ONLY: locate
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: splint
     END FUNCTION splint
  END INTERFACE

  INTERFACE spline_interface
     SUBROUTINE spline(x,y,yp1,ypn,y2)
       USE nrtype; USE nrutil, ONLY : assert_eq
       USE nr, ONLY : tridag
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(DP), INTENT(IN) :: yp1,ypn
       REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
     END SUBROUTINE spline
  END INTERFACE

  INTEGER, INTENT(IN) :: iSmoothMove
  INTEGER :: ierr, k, i, j
  REAL(DP) :: distance1deriv1, distance1derivn
  REAL(DP), DIMENSION(npsi) :: xOld, yOld, zOld, psiOld, distanceOld, x2derivs, y2derivs, &
       & z2derivs 
  REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
  REAL(DP) :: blend

  blend = 0.1_dp**iPsiMove

  IF (iSmoothMove /= 0 .AND. iPsiMove > 1) GOTO 10

  xPrev = x
  yPrev = y
  zPrev = z

  ierr = 0

  kloop: DO k = 2, nzeta
     iloop: DO i = 1, nthe
        xOld(1:npsi) = x(i,1:npsi,k)
        yOld(1:npsi) = y(i,1:npsi,k)
        zOld(1:npsi) = z(i,1:npsi,k)
        psiOld(1:npsi) = psi(i,1:npsi,k)

        CALL spline(psiold, xOld, 1.E31_dp, 1.E31_dp, x2derivs)
        CALL spline(psiold, yOld, 1.E31_dp, 1.E31_dp, y2derivs)
        CALL spline(psiold, zOld, 1.E31_dp, 1.E31_dp, z2derivs)

        DO j = 2, npsim
           x(i,j,k) = splint(psiold, xOld, x2derivs, psival(j))
           y(i,j,k)  = splint(psiold, yOld, y2derivs, psival(j))
           z(i,j,k)  = splint(psiold, zOld, z2derivs, psival(j))
           psi(i,j,k) = psival(j)  
        END DO
     END DO iloop
  END DO kloop

  ! periodic boundary condition in zeta

  x(:,:,1) = x(:,:,nzeta)
  y(:,:,1) = y(:,:,nzeta)
  z(:,:,1) = z(:,:,nzeta)
  x(:,:,nzetp) = x(:,:,2)
  y(:,:,nzetp) = y(:,:,2)
  z(:,:,nzetp) = z(:,:,2)

  DO j = 1, npsi
     psi(:,j, 1) = psival(j)
     psi(:,j,nzetp) = psival(j)
  END DO

  !C diffmx = MAXVAL(sqrt((x-xPrev)**2 + (y-yPrev)**2 + (z-zPrev)**2))

10 IF (iSmoothMove /= 0 .AND. iPsiMove > 1) THEN
     ! Add these in difficult equilibria
     PRINT*, 'mappsi: blend = ', blend
     x = 1*blend*x + (1.-1*blend)*xPrev 
     y = 1*blend*y + (1.-1*blend)*yPrev 
     z = 1*blend*z + (1.-1*blend)*zPrev 
     !C z(nThetaEquator,:,:) = 0._dp ! Symmetry
  END IF


  RETURN

END SUBROUTINE mapPsi


