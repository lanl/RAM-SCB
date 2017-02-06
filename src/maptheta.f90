SUBROUTINE mapTheta

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
       INTEGER(I4B) :: khi,klo,n
       REAL(DP) :: a,b,h
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
       INTEGER(I4B) :: n
       REAL(DP), DIMENSION(SIZE(x)) :: a,b,c,r
     END SUBROUTINE spline
  END INTERFACE

  REAL(DP), DIMENSION(nthe) :: xOld, yOld, zOld, distance, chiValOld, chi2derivsX, &
       chi2derivsY, chi2derivsZ
  INTEGER :: i, j, k

  !     now move theta coordinates along each surface
  !     equal arc length along the i grids

  zetaloop: DO k = 2, nzeta
     fluxloop: DO j = 1, npsi
        distance(1) = 0._dp
        xOld(:) = x(1:nthe,j,k)
        yOld(:) = y(1:nthe,j,k)
        zOld(:) = z(1:nthe,j,k)
        chiValOld(1) = 0._dp

        DO i = 2, nthe
           distance(i) = distance(i-1) + SQRT((x(i,j,k)-x(i-1,j,k))**2 &
                & +(y(i,j,k)-y(i-1,j,k))**2 +(z(i,j,k)-z(i-1,j,k))**2)
        END DO

        chiValOld = distance / distance(nthe) * pi_d

        ! "Natural" splines

        CALL spline(chiValOld, xOld, 1.E31_dp, 1.E31_dp, chi2derivsX)
        CALL spline(chiValOld, yOld, 1.E31_dp, 1.E31_dp, chi2derivsY)
        CALL spline(chiValOld, zOld, 1.E31_dp, 1.E31_dp, chi2derivsZ)

        DO i = 2, nthem
           x(i,j,k) = splint(chiValOld, xOld, chi2derivsX, chiVal(i))
           y(i,j,k) = splint(chiValOld, yOld, chi2derivsY, chiVal(i))
           z(i,j,k) = splint(chiValOld, zOld, chi2derivsZ, chiVal(i))
        END DO
     END DO fluxloop
  END DO zetaloop
  !
  !  periodic boundary conditions
  !

  x(:,:,1) = x(:,:,nzeta)
  y(:,:,1) = y(:,:,nzeta)
  z(:,:,1) = z(:,:,nzeta)
  x(:,:,nzetp) = x(:,:,2)
  y(:,:,nzetp) = y(:,:,2)
  z(:,:,nzetp) = z(:,:,2)

!C z(nThetaEquator,:,:) = 0.0_dp ! Symmetry


  RETURN

END SUBROUTINE mapTheta
                

