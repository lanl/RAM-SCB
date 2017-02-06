SUBROUTINE Spline_derivs_1D_zeta(x1, f, func, derivsX1)

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
  USE nrtype, ONLY : DP

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8 = DP
  REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8

  REAL(r8), DIMENSION(:), INTENT(IN) :: x1 ! independent variable
  REAL(r8), DIMENSION(:), INTENT(IN) :: f
  REAL(r8), DIMENSION(:), INTENT(OUT) :: func, derivsX1 ! interpolated values
  INTEGER n1, ier, BCS1(2), j
  TYPE(EZspline1_r8) :: spline_o ! 1-D EZspline object

  INTEGER k1
  REAL(r8), DIMENSION(:), ALLOCATABLE :: z1

  n1 = SIZE(x1)

  BCS1 = (/ -1, -1 /) ! periodic spline for zeta
  ! BCS1 = (/ 0, 0 /) ! "not-a-knot" spline 
  ! BCS1 = (/ 2, 2 /)    ! "natural" spline 

  ! initialize/allocate memory
  !  PRINT *,'initializing...'

  CALL EZspline_init(spline_o, n1, BCS1, ier)
  CALL EZspline_error(ier)

  spline_o%x1 = x1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
  ! and is aliased to x1 now

  ! need to set explicitly the following if boundary conditions 
  ! are anything but not-a-knot or periodic, i.e. BCS(n) has
  ! element /= -1, 0.
  ! spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs

  !  PRINT *,'setting up 1D spline ...'

  CALL EZspline_setup(spline_o, f, ier)
  CALL EZspline_error(ier)

  ! save object

  CALL EZspline_save(spline_o, "spline.nc", ier)
  CALL EZspline_error(ier)

  ! Array interpolation


  k1 = SIZE(derivsX1)

  ALLOCATE(z1(k1))

  z1 = x1

  !  PRINT *,'array interpolation 1D spline ...'

  CALL EZspline_interp(spline_o, k1, z1, func, ier)
  CALL EZspline_error(ier)

  !  PRINT*, 'derivatives 1D spline ...'

  CALL EZspline_derivative(spline_o, 1, k1, z1, derivsX1, ier)
  ! partial wrt the coordinate for all points

  ! clean up and free up memory

  !  PRINT *,'cleaning up'

  CALL Ezspline_free(spline_o, ier)
  CALL EZspline_error(ier)

  DEALLOCATE(z1)

  RETURN
END SUBROUTINE Spline_derivs_1D_zeta
