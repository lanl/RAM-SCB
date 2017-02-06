MODULE Module_inversion

CONTAINS

  SUBROUTINE MIGS (A,N,X,INDX)
    !
    ! Subroutine for matrix inversion. A(N,N) is the matrix to be 
    ! inverted, and the inverse is stored in X(N,N) in the output.
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: NMAX=100
    INTEGER, INTENT (IN) :: N
    INTEGER :: I,J,K
    INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
    REAL, INTENT (INOUT), DIMENSION (N,N):: A
    REAL, INTENT (OUT), DIMENSION (N,N):: X
    REAL, DIMENSION (NMAX,NMAX) :: B
    LOGICAL :: IR
    !
    IF(N.GT.NMAX) STOP 'The matrix dimension is too large.'
    !

    DO I = 1, N
       DO J = 1, N
          B(I,J) = 0.0
       END DO
    END DO
    DO I = 1, N
       B(I,I) = 1.0
    END DO
    !
    IR = .FALSE.
    CALL ELGS (A,N,IR,INDX)
    !
    DO I = 1, N-1
       DO J = I+1, N
          DO K = 1, N
             B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
          END DO
       END DO
    END DO
    !
    DO I = 1, N
       X(N,I) = B(INDX(N),I)/A(INDX(N),N)
       DO J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
          DO K = J+1, N
             X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
          END DO
          X(J,I) =  X(J,I)/A(INDX(J),J)
       END DO
    END DO
  END SUBROUTINE MIGS
!

SUBROUTINE ELGS (A,N,IR,INDX)
!
! Subroutine for the partial pivoting Gaussian elimination.
! A is the coefficient matrix in the input and the transformed
! matrix in the output.  INDX records the pivoting order.  IR
! is the indicator for the rescaling.  Copyright (c) Tao Pang 1997.
!
   IMPLICIT NONE
   INTEGER, PARAMETER :: NMAX=100
   INTEGER, INTENT (IN) :: N
   INTEGER :: I,J,K,ITMP
   INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
   REAL :: C1,PI,PI1,PJ
   REAL, INTENT (INOUT), DIMENSION (N,N) :: A
   REAL, DIMENSION (NMAX) :: C
   LOGICAL, INTENT (IN) :: IR
!
  IF (N.GT.NMAX) STOP 'The matrix dimension is too large.'
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Rescale the coefficients
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      C1 = AMAX1(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
  DO I = 1, N
    DO J = 1, N
      A(I,J) = A(I,J)/C(I)
    END DO
  END DO
!
! Search the pivoting elements
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Eliminate the elements below the diagonal
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
      A(INDX(I),J) = PJ
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
! Rescale the coefficients back to the original magnitudes
!
  IF (IR) THEN
    DO I = 1, N
      DO J = 1, N
        A(I,J) = C(I)*A(I,J)
      END DO
    END DO
  END IF
END SUBROUTINE ELGS

END MODULE Module_inversion
