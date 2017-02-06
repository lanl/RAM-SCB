SUBROUTINE psiges
  !   initial guess of poloidal flux
 
  USE nrtype
  USE Module1

  INTEGER :: i, j, k

  !     initial guess for psi
  DO  k=1,nzetp
     DO  j=1,npsi
        DO  i=1,nthe
           psi(i,j,k)=psival(j)
        END DO
     END DO
  END DO

  RETURN

END SUBROUTINE psiges
