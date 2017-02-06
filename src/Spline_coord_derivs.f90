SUBROUTINE Spline_coord_derivs(x_1, x_2, x_3, f_input, derivsX1, derivsX2, derivsX3)

  USE EZspline_obj ! import the modules
  USE EZspline  
  USE Module1 ! for nthe, npsi, nzeta

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8 = DP

  REAL(r8), DIMENSION(:), INTENT(IN) :: x_1  
  REAL(r8), DIMENSION(:), INTENT(IN) :: x_2
  REAL(r8), DIMENSION(:), INTENT(IN) :: x_3 ! independent variables
  REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f_input
  REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: derivsX1, derivsX2, derivsX3 
  ! derivative values
  INTEGER       :: n1, n2, n3, ier, ierrgrad, BCS1(2), BCS2(2), BCS3(2)
  TYPE(EZspline3_r8) :: spline_o ! 3-D EZspline object

  REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: grad

  n1 = SIZE(x_1)
  n2 = SIZE(x_2)
  n3 = SIZE(x_3)


  BCS1 = (/ 0, 0 /) ! "not-a-knot" spline in theta
  !   BCS1 = (/ 2, 2 /)    ! "natural" spline in theta
  BCS2 = (/ 0, 0 /) ! "not-a-knot" spline in rho
  !   BCS2 = (/ 2, 2 /)    ! "natural" spline in rho
  ! BCS2 = (/ 1, 1 /) ! first derivative imposed in rho
  BCS3 = (/ -1, -1 /)  ! periodic in zeta
  ! BCS3 = (/ 0, 0 /)

  ! initialize/allocate memory
  ! PRINT *,'initializing...'

  CALL EZspline_init(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
  ! CALL EZLinear_init(spline_o, n1, n2, n3, ier)
  CALL EZspline_error(ier)

  spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]
  spline_o%x2 = x_2    ! necessary if spline_o%x2 not in [0, 2 pi]

  ! set explicitely spline_o%x3

  spline_o%x3 = x_3

  spline_o%isHermite = 1 ! Akima spline


  ! need to set explicitly the following if boundary conditions 
  ! are anything but not-a-knot or periodic, i.e. BCS(n) has
  ! element /= -1, 0.
  ! spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
  ! spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs
  ! spline_o%bcval3min = ... ; spline_o%bcval3max = ...


  CALL EZspline_setup(spline_o, f_input, ier)
  CALL EZspline_error(ier)

  ! save object

  ! The save to disk for all processors is extremely slow with MPI - solved by removing the save (is not necessary since it is defined for each call)
  !C CALL EZspline_save(spline_o, "spline.nc", ier)
  CALL EZspline_error(ier)


  ! Array interpolation

  if (.not. allocated(grad)) ALLOCATE(grad(n1, n2, n3, 3), stat = ierrgrad)

  !  PRINT *,'array interpolation...'

  !  CALL EZspline_interp(spline_o, k1, k2, k3, z1, z2, z3, fz, ier)
  !  CALL EZspline_error(ier)

  !  PRINT*, 'gradient is being called ...'

  CALL EZspline_gradient(spline_o, n1, n2, n3, x_1, x_2, x_3, grad, ier)

  derivsX1 = grad(:, :, :, 1)
  derivsX2 = grad(:, :, :, 2)
  derivsX3 = grad(:, :, :, 3)

  DEALLOCATE(grad)

  CALL Ezspline_free(spline_o, ier)
  CALL EZspline_error(ier)

  RETURN

END SUBROUTINE Spline_coord_derivs
