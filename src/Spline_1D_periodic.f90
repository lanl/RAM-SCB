SUBROUTINE Spline_1D_periodic(x_1, f, x, func, ierrDomain)

  USE nrtype, ONLY : DP
  USE EZspline_obj ! import the modules
  USE EZspline  

  IMPLICIT NONE

  INTEGER, PARAMETER :: r8 = DP 
  REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8

  REAL(r8), DIMENSION(:), INTENT(IN) :: x_1 ! independent variable
  REAL(r8), DIMENSION(:), INTENT(IN) :: f
  REAL(r8), DIMENSION(:), INTENT(OUT) :: func   ! interpolated values 
  INTEGER           :: n1, n_x, ier, BCS1(2)
  TYPE(EZspline1_r8) :: spline_o ! 1-D EZspline object

  INTEGER           :: j
  INTEGER, INTENT(OUT) :: ierrDomain
  REAL(r8), DIMENSION(:), INTENT(IN)    :: x  ! grid of points for output

  n1 = SIZE(x_1)

  ! Boundary conditions for interpolation

  BCS1 = (/-1, -1/)  ! periodic spline 
  ! BCS2 = (/ 0, 0 /) ! "not-a-knot" spline 
  ! BCS2 = (/ 2, 2 /) !  "natural" spline

  ! initialize/allocate memory
  !  PRINT *,'initializing...'

  CALL EZspline_init(spline_o, n1, BCS1, ier)
  !C CALL EZLinear_init(spline_o, n1, n2, ier)
  CALL EZspline_error(ier)

  spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
  ! and is aliased to x_1 now
  spline_o%isHermite = 1 ! Akima spline; smoother, more "natural"  interpolation (see Akima's paper)

    ! need to set explicitly the following if boundary conditions 
  ! are anything but not-a-knot or periodic, i.e. BCS(n) has
  ! element /= -1, 0.
!   spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
!  spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs


  !  PRINT *,'setting up 1D spline ...'

  CALL EZspline_setup(spline_o, f, ier)
  CALL EZspline_error(ier)

  ! save object

  !C CALL EZspline_save(spline_o, "spline.nc", ier)  ! Not necessary to save it if not loading again; slow by NFS on hyrax nodes
  CALL EZspline_error(ier)

  ! Cloud interpolation

  n_x = SIZE(x,1)

  DO j = 1, n_x
     CALL EZspline_isInDomain(spline_o, x(j), ierrDomain)
     IF (ierrDomain <= 0) THEN
        CALL EZspline_interp(spline_o, x(j), func(j), ier)
        CALL EZspline_error(ier)
     END IF
  END DO

  !C PRINT *,'Spline_1D_point: cleaning up'
  
  CALL Ezspline_free(spline_o, ier)
  CALL EZspline_error(ier)
  
  
  RETURN
  
END SUBROUTINE Spline_1D_periodic
