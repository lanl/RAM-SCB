SUBROUTINE newlin2(xp,yp,zp,xpn,ypn,zpn,s_p,spn,np,np2)
  !     interpolate to new position (xpn,ypn,zpn)
  !     linear version
  ! np2(old) greater or equal to np

  USE nrtype
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

  INTEGER, INTENT(IN)                          :: np
  INTEGER, INTENT(IN)                          :: np2
  REAL(DP), INTENT(IN)                         :: xp(np2)
  REAL(DP), INTENT(IN)                         :: yp(np2)
  REAL(DP), INTENT(IN)                         :: zp(np2)
  REAL(DP), INTENT(OUT)                        :: xpn(np)
  REAL(DP), INTENT(OUT)                        :: ypn(np)
  REAL(DP), INTENT(OUT)                        :: zpn(np)
  REAL(DP), INTENT(IN)                         :: s_p(np2)
  REAL(DP), INTENT(IN)                         :: spn(np)

  REAL(DP), DIMENSION(np2) :: distance2derivsX, distance2derivsY, distance2derivsZ
  INTEGER :: j

  xpn(1) = xp(1)
  ypn(1) = yp(1)
  zpn(1) = zp(1)

  CALL spline(s_p, xp, 1.E31_dp, 1.E31_dp, distance2derivsX)
  CALL spline(s_p, yp, 1.E31_dp, 1.E31_dp, distance2derivsY)
  CALL spline(s_p, zp, 1.E31_dp, 1.E31_dp, distance2derivsZ)

  DO j = 2, np
     xpn(j) = splint(s_p, xp, distance2derivsX, spn(j))
     ypn(j) = splint(s_p, yp, distance2derivsY, spn(j))
     zpn(j) = splint(s_p, zp, distance2derivsZ, spn(j))
  END DO

  RETURN
END SUBROUTINE newlin2
