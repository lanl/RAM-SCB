SUBROUTINE dmat_solve (n, nrhs, a, info )

  !*******************************************************************************
  !
  !! DMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
  !
  !  Modified:
  !
  !    29 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer N, the order of the matrix.
  !
  !    Input, integer NRHS, the number of right hand sides.  NRHS
  !    must be at least 0.
  !
  !    Input/output, real ( kind = 8 ) A(N,N+NRHS), contains in rows and columns 1
  !    to N the coefficient matrix, and in columns N+1 through
  !    N+NRHS, the right hand sides.  On output, the coefficient matrix
  !    area has been destroyed, while the right hand sides have
  !    been overwritten with the corresponding solutions.
  !
  !    Output, integer INFO, singularity flag.
  !    0, the matrix was not singular, the solutions were computed;
  !    J, factorization failed on step J, and the solutions could not
  !    be computed.
  !
 
  use nrtype, ONLY : DP

 IMPLICIT NONE

  INTEGER, intent(IN) :: n
  INTEGER, intent(IN) :: nrhs

  REAL(DP), intent(IN OUT) :: a(n,n+nrhs)
  REAL(DP) :: apivot
  REAL(DP) :: factor
  INTEGER i
  INTEGER info
  INTEGER ipivot
  INTEGER j

  info = 0

  DO j = 1, n
     !
     !  Choose a pivot row.
     !
     ipivot = j
     apivot = a(j,j)

     DO i = j+1, n
        IF ( ABS ( apivot ) < ABS ( a(i,j) ) ) THEN
           apivot = a(i,j)
           ipivot = i
        END IF
     END DO

     IF ( apivot == 0.0D+00 ) THEN
        info = j
        RETURN
     END IF
     !
     !  Interchange.
     !
     DO i = 1, n + nrhs
        CALL d_swap ( a(ipivot,i), a(j,i) )
     END DO
     !
     !  A(J,J) becomes 1.
     !
     a(j,j) = 1.0D+00
     a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
     !
     !  A(I,J) becomes 0.
     !
     DO i = 1, n

        IF ( i /= j ) THEN

           factor = a(i,j)
           a(i,j) = 0.0D+00
           a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)

        END IF

     END DO

  END DO

  RETURN
END SUBROUTINE dmat_solve
