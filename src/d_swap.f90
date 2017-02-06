subroutine d_swap ( x, y )

!*******************************************************************************
!
!! D_SWAP switches two real values.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
 
  use nrtype, ONLY : DP
  
 implicit none

  real (DP) :: x
  real (DP) :: y
  real (DP) :: z

  z = x
  x = y
  y = z

  return
end
