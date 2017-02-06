SUBROUTINE mapAlpha(iSmoothMove)
  ! new cubic spline interpolation, without involving linear distance calculation

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

  !  ignore i = 1 & nthe lines since the boundary condition is delta = 0 there (on Earth surface)

  INTEGER, INTENT(IN) :: iSmoothMove
  REAL(DP), DIMENSION(nzetp) :: xOld, yOld, zOld, alfaOld, x2derivs, y2derivs, &
       & z2derivs   
  REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: xPrev, yPrev, zPrev
  REAL(DP) :: blend 
  INTEGER :: i, j, k, ierr

  blend = 0.1_dp**iAlphaMove

  IF (iSmoothMove /= 0 .AND. iAlphaMove > 1) GOTO 10

  xPrev = x
  yPrev = y
  zPrev = z

  ierr = 0

  jloop : DO j = 1, npsi
     iloop: DO i = 2, nthem

        xOld(1:nzetp) = x(i,j,1:nzetp)
        yOld(1:nzetp) = y(i,j,1:nzetp)
        zOld(1:nzetp) = z(i,j,1:nzetp)
        alfaOld(1:nzetp) = alfa(i,j,1:nzetp) 

        CALL spline(alfaOld, xOld, 1.E31_dp, 1.E31_dp, x2derivs)
        CALL spline(alfaOld, yOld, 1.E31_dp, 1.E31_dp, y2derivs)
        CALL spline(alfaOld, zOld, 1.E31_dp, 1.E31_dp, z2derivs)

        DO k = 2,nzeta
           x(i,j,k) = splint(alfaOld, xOld, x2derivs, alphaval(k))
           y(i,j,k) = splint(alfaOld, yOld, y2derivs, alphaval(k))
           z(i,j,k) = splint(alfaOld, zOld, z2derivs, alphaval(k))
           alfa(i,j,k) = alphaval(k)
        END DO

        ! Periodic boundary conditions
        x(i,j,1) = x(i,j,nzeta)
        y(i,j,1) = y(i,j,nzeta)
        z(i,j,1) = z(i,j,nzeta)
        alfa(i,j,1) = alfa(i,j,nzeta) - 2._dp*pi_d 
        x(i,j,nzetp) = x(i,j,2)
        y(i,j,nzetp) = y(i,j,2)
        z(i,j,nzetp) = z(i,j,2)
        alfa(i,j,nzetp) = alfa(i,j,2) + 2._dp*pi_d 

     END DO iloop
  END DO jloop

10 IF (iSmoothMove /= 0 .AND. iAlphaMove > 1) THEN
     ! Add these in difficult equilibria
     PRINT*, 'mapalpha: blend = ', blend
     x = 1*blend*x + (1.-1*blend)*xPrev 
     y = 1*blend*y + (1.-1*blend)*yPrev 
     z = 1*blend*z + (1.-1*blend)*zPrev 
  END IF


  RETURN

END SUBROUTINE mapAlpha








