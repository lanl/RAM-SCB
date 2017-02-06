SUBROUTINE polint(xa,ya,x,y,dy)
  USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(OUT) :: y, dy
  INTEGER(I4B) :: m,n,ns
  REAL(DP), DIMENSION(SIZE(xa)) :: c,d,den,ho
  n = assert_eq(SIZE(xa),SIZE(ya),'polint')
  c=ya
  d=ya
  ho=xa-x
  ns=iminloc(ABS(x-xa))
  y=ya(ns)
  ns=ns-1
  DO m=1,n-1
     den(1:n-m)=ho(1:n-m)-ho(1+m:n)
     IF (ANY(den(1:n-m) == 0.0)) THEN
        PRINT*, 'xa = ', xa
        PRINT*, 'ya = ', ya
        PRINT*, 'x = ', x
        PRINT*, 'y = ', y
        PAUSE 'Problem'
        CALL nrerror('polint: calculation failure')
     END IF
     den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
     d(1:n-m)=ho(1+m:n)*den(1:n-m)
     c(1:n-m)=ho(1:n-m)*den(1:n-m)
     IF (2*ns < n-m) THEN
        dy=c(ns+1)
     ELSE
        dy=d(ns)
        ns=ns-1
     END IF
     y=y+dy
  END DO

return
END SUBROUTINE polint
