SUBROUTINE InterpolatePsi   

  ! Interpolates the new values of psi at the new locations xnew on midnight 
  ! equator

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

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: xEqmid, psival1D, psi2Deriv
  INTEGER :: ialloc, ierr, j

  ALLOCATE(xEqmid(npsi), stat = ialloc)
  ALLOCATE(psival1D(npsi), stat = ialloc)
  ALLOCATE(psi2Deriv(npsi), stat = ialloc)

  xEqmid(1:npsi) = x(nThetaEquator, 1:npsi, nZetaMidnight)
  psival1D(1:npsi) = psival(1:npsi)
 
  CALL spline(xEqmid, psival1D, 1.E31_dp, 1.E31_dp, psi2Deriv)

  DO j = 2, npsi-1
     psival(j) = splint(xEqmid, psival1D, psi2Deriv, xEqmidNew(j))
  END DO

  IF (ALLOCATED(xEqmid)) DEALLOCATE(xEqmid, stat = ierr)
  IF (ALLOCATED(psival1D)) DEALLOCATE(psival1D, stat = ierr)
  IF (ALLOCATED(psi2Deriv)) DEALLOCATE (psi2Deriv, stat = ierr)


  RETURN
END SUBROUTINE InterpolatePsi


