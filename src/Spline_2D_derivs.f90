SUBROUTINE Spline_2D_derivs(x1, x2, f, derivsX1, derivsX2)

  ! example of ezspline calls
  ! -------------------------
  ! On PPPL cluster machine use the following command at compile:
  ! On Alpha-Linux:
  ! f90 -assume no2underscores -I/usr/ntcc/mod -o spline_test spline_test.f90 -L/usr/ntcc/lib -lpspline -lezcdf -L/usr/local/lib -lnetcdf
  ! On Alpha-OSF1:
  ! f90 -I/usr/ntcc/mod -o spline_test spline_test.f90 -L/usr/ntcc/lib -lpspline -lezcdf -L/usr/local/lib -lnetcdf
  ! On Solaris:
  ! f90 -M/usr/ntcc/mod -dalign -o spline_test spline_test.f90 -L/usr/ntcc/lib -lpspline -lezcdf -L/usr/local/lib -lnetcdf
  ! On Linux (Fujitsu):
  ! f90 -Am -g -I /usr/ntcc/ffc/mod -o spline_test spline_test.f90 -L/usr/ntcc/ffc/lib -lpspline -lezcdf -L/usr/local/ffc/lib -lnetcdf

  USE EZspline_obj ! import the modules
  USE EZspline  
  use nrtype, ONLY : DP

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8 = DP 
  REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8

  REAL(r8), DIMENSION(:), INTENT(IN) :: x1, x2 ! independent variable
  REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
  REAL(r8), DIMENSION(:,:) :: derivsX1, derivsX2  ! derivative values

  REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: grad

  INTEGER           :: n1, n2, ier, BCS1(2), BCS2(2), j
  TYPE(EZspline2_r8) :: spline_o ! 1-D EZspline object
  INTEGER :: ierralloc, ierrdealloc


  n1 = SIZE(x1)
  n2 = SIZE(x2)

  IF (.NOT. ALLOCATED(grad)) ALLOCATE(grad(n1,n2,2), stat = ierralloc)

  ! Boundary conditions for interpolation
  BCS1 = (/ 0, 0 /) ! "not-a-knot" spline in first coordinate
  !  BCS1 = (/ 2, 2 /)    ! "natural" spline in first coordinate
  BCS2 = (/ -1, -1 /) ! periodic spline in second coordinate (zeta)

  ! initialize/allocate memory
  !  PRINT *,'initializing...'

  ! CALL EZspline_init(spline_o, n1, n2, BCS1, BCS2, ier)
  CALL EZlinear_init(spline_o, n1, n2, ier)
  CALL EZspline_error(ier)

  spline_o%x1 = x1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
  ! and is aliased to x1 now
  spline_o%x2 = x2

  ! need to set explicitly the following if boundary conditions 
  ! are anything but not-a-knot or periodic, i.e. BCS(n) has
  ! element /= -1, 0.
  ! spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
  ! spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs

  !  PRINT *,'setting up 2D spline ...'

  CALL EZspline_setup(spline_o, f, ier)
  CALL EZspline_error(ier)

  ! save object

  !C  CALL EZspline_save(spline_o, "spline.nc", ier)  ! Not necessary, it is define at each call
  CALL EZspline_error(ier)

  ! Cloud interpolation

  !  PRINT *,'cloud interpolation ...'
  !  CALL EZspline_interp(spline_o, k, x, y, func, ier)
  !  CALL EZspline_error(ier)

  !  PRINT*, 'derivatives 2D spline ...'
  CALL EZspline_gradient(spline_o, n1, n2, x1, x2, grad, ier)
  ! partial wrt the coordinate for all points
  derivsX1 = grad(:,:,1)
  derivsX2 = grad(:,:,2)
 
 ! clean up and free up memory
  DEALLOCATE(grad, stat = ierrdealloc)

  CALL Ezspline_free(spline_o, ier)
  CALL EZspline_error(ier)


  RETURN

END SUBROUTINE Spline_2D_derivs
