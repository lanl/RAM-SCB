MODULE module_Gauss
  IMPLICIT NONE
CONTAINS
  FUNCTION matrixInverse(matrix)
    ! This function finds solution to set of linear simultaneous equations
    ! using Gauss-Jordan elimination.  It uses all array handling capabilities
    ! of Fortran 90 and effectively does the back substitution at the same
    ! time as the elimination.  It uses the intrinsic function SPREAD which,
    ! as it is used here produces a vector containing n or n+1 elements
    ! each with the value i.

    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(IN) :: matrix
    ! Find number of unknowns/equations and set shape of function, working
    ! array and logical mask array
    REAL, DIMENSION(SIZE(matrix,1)) :: gauss_elim_3
    REAL, DIMENSION(SIZE(matrix,1),SIZE(matrix,2)) :: work
    REAL, DIMENSION(SIZE(matrix,1),SIZE(matrix,2)-1) :: matrixInverse
    LOGICAL, DIMENSION(SIZE(matrix,1),SIZE(matrix,2)) :: not_pivot
    INTEGER :: n, i
    ! Check shape of matrix and if incorrect exit
    IF(SIZE(matrix,2) /= (SIZE(matrix,1)+1)) THEN
       STOP 'Bad matrix shape'
    END IF
    n = SIZE(matrix,1)
    ! Copy matrix to working array

    work = matrix
    DO i = 1,n
       ! Do elimination checking for near zero pivot

       IF (ABS(work(i,i)) <= 1.0e-9) THEN
          STOP 'Zero encountered in pivot'
       END IF
       ! Normalize i'th equation

       work(i,:) = work(i,:)/work(i,i)
       not_pivot = .TRUE.
       ! Set mask array false along row i

       not_pivot(i,:) = .FALSE.
       ! Operate on rows below i'th

       WHERE (not_pivot)
          ! This is where everything happens!  You'll need to think a bit about
          ! what is going on.
          work = work - work(:,SPREAD(i,1,n+1))*work(SPREAD(i,1,n),:)
       END WHERE
    END DO
    gauss_elim_3 = work(:,n+1)
    matrixInverse = work(:, 1:n)

  END FUNCTION matrixInverse
END MODULE module_Gauss

