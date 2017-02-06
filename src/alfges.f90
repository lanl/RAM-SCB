SUBROUTINE alfges
  !   initial guess of alpha 

  USE nrtype
  USE Module1

  INTEGER :: k

  DO k = 1, nzetp
     alfa(1:nthe,1:npsi,k) = alphaval(k)
  END DO


  RETURN

END SUBROUTINE alfges
