SUBROUTINE secant(f,x0,x1,f1,xtol,ftol,ntol,iflag)

  USE nrtype

  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT)                     :: x0
  REAL(DP), INTENT(IN OUT)                     :: x1
  REAL(DP), INTENT(OUT)                        :: f1
  REAL(DP), INTENT(IN)                         :: xtol
  REAL(DP), INTENT(IN)                         :: ftol
  INTEGER, INTENT(IN)                      :: ntol
  INTEGER, INTENT(OUT)                     :: iflag
  REAL(DP), EXTERNAL :: f

  REAL(DP) :: f0, deltax, deltaf
  INTEGER :: n

  iflag=0 
  f0=f(x0) 
  deltax=x1 -x0
  CYCLE:  DO  n=1,ntol
     f1=f(x1)
     IF ( ABS(f1) < ftol) THEN
        iflag = 1
        RETURN
     END IF
     deltaf=f0-f1
     IF (deltaf == 0.0) EXIT CYCLE
     deltax=f1/deltaf*deltax
     x0=x1
     x1=x1+deltax
     IF ( ABS(deltax) < xtol) RETURN
     f0=f1
  END DO CYCLE
  !c  iflag=1 or 2, the convergence is bad
999 iflag=2
  RETURN
END SUBROUTINE secant
!*****************************
