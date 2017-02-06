SUBROUTINE Spline_2D_array(x_1, x_2, f, x, y, func)

  USE EZspline_obj ! import the modules
  USE EZspline  
  use nrtype, ONLY : DP

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8 = DP 
  REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8

  REAL(r8), DIMENSION(:), INTENT(IN) :: x_1, x_2 ! independent variable
  REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
  REAL(r8), DIMENSION(:,:), INTENT(OUT) :: func   ! interpolated values and derivs.
  INTEGER           :: n1, n2, n_x, n_y, ier, BCS1(2), BCS2(2)
  TYPE(EZspline2_r8) :: spline_o ! 2-D EZspline object

  INTEGER           :: j, k
  REAL(r8), DIMENSION(:), INTENT(IN)    :: x, y  ! arrays of points for output

  n1 = SIZE(x_1)
  n2 = SIZE(x_2)

  ! Boundary conditions for interpolation
  BCS1 = (/ 0, 0 /) ! "not-a-knot" spline in first coordinate
  ! BCS1 = (/ 2, 2 /)    ! "natural" spline 

  BCS2 = (/ 0, 0 /) ! "not-a-knot" spline in second coordinate
  ! BCS2 = (/ 2, 2 /) !  "natural" spline

  ! initialize/allocate memory
  !  PRINT *,'initializing...'

  CALL EZspline_init(spline_o, n1, n2, BCS1, BCS2, ier)
  CALL EZspline_error(ier)


  spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
  ! and is aliased to x_1 now
  spline_o%x2 = x_2


  ! need to set explicitly the following if boundary conditions 
  ! are anything but not-a-knot or periodic, i.e. BCS(n) has
  ! element /= -1, 0.
  spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
  spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs


  !  PRINT *,'setting up 2D spline ...'

  CALL EZspline_setup(spline_o, f, ier)
  CALL EZspline_error(ier)

  ! save object

  CALL EZspline_save(spline_o, "spline.nc", ier)
  CALL EZspline_error(ier)

  ! Cloud interpolation

  n_x = SIZE(x)
  n_y = SIZE(y)   ! x and y have the same shape and size

  CALL EZspline_interp(spline_o, n_x, n_y, x, y, func, ier)
  CALL EZspline_error(ier)

  !  PRINT *,'cleaning up'

  CALL Ezspline_free(spline_o, ier)
  CALL EZspline_error(ier)


  RETURN
END SUBROUTINE Spline_2D_array
