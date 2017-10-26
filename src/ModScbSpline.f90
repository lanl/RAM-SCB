MODULE ModScbSpline
  ! Contains subroutines related to calculating splines for the scb code
  implicit none
  
  contains
!==============================================================================
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

!==============================================================================
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
  
!==============================================================================
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
  
!==============================================================================
  SUBROUTINE Spline_2D_periodic(x_1, x_2, f, x, y, func, ierrDomain)
  
    USE EZspline_obj ! import the modules
    USE EZspline  
    use nrtype, ONLY : DP
  
    IMPLICIT NONE
  
    INTEGER, PARAMETER :: r8 = DP 
    REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8
  
    REAL(r8), DIMENSION(:), INTENT(IN) :: x_1, x_2 ! independent variable
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(r8), DIMENSION(:,:), INTENT(OUT) :: func   ! interpolated values 
    INTEGER           :: n1, n2, n_x, n_y, ier, BCS1(2), BCS2(2)
    TYPE(EZspline2_r8) :: spline_o ! 2-D EZspline object
  
    INTEGER           :: j, k
    INTEGER, INTENT(OUT) :: ierrDomain
    REAL(r8), DIMENSION(:,:), INTENT(IN)    :: x, y  ! grid of points for output
  
    n1 = SIZE(x_1)
    n2 = SIZE(x_2)
  
    ! Boundary conditions for interpolation
    BCS1 = (/ 0, 0 /) ! "not-a-knot" spline in first coordinate
    !C BCS1 = (/ 2, 2 /)    ! "natural" spline 
  
    BCS2 = (/-1, -1/)  ! periodic spline in second coordinate
    ! BCS2 = (/ 0, 0 /) ! "not-a-knot" spline in second coordinate
    ! BCS2 = (/ 2, 2 /) !  "natural" spline
  
    ! initialize/allocate memory
    !  PRINT *,'initializing...'
  
    CALL EZspline_init(spline_o, n1, n2, BCS1, BCS2, ier)
    !C CALL EZLinear_init(spline_o, n1, n2, ier)
    CALL EZspline_error(ier)
  
    spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
    ! and is aliased to x_1 now
    spline_o%x2 = x_2
    spline_o%isHermite = 1 ! Akima spline; smoother, more "natural"  interpolation (see Akima's paper)
  
      ! need to set explicitly the following if boundary conditions 
    ! are anything but not-a-knot or periodic, i.e. BCS(n) has
    ! element /= -1, 0.
  !   spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
  !  spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs
  
  
    !  PRINT *,'setting up 2D spline ...'
  
    CALL EZspline_setup(spline_o, f, ier)
    CALL EZspline_error(ier)
  
    ! save object
  
    !C CALL EZspline_save(spline_o, "spline.nc", ier)  ! Not necessary to save it if not loading again; slow by NFS on hyrax nodes
    CALL EZspline_error(ier)
  
    ! Cloud interpolation
  
    n_x = SIZE(x,1)
    n_y = SIZE(x,2)   ! x and y have the same shape and size
  
    DO j = 1, n_x
       DO k = 1, n_y
          CALL EZspline_isInDomain(spline_o, x(j,k), y(j,k), ierrDomain)
          IF (ierrDomain <= 0) THEN
             CALL EZspline_interp(spline_o, x(j,k), y(j,k), func(j,k), ier)
             CALL EZspline_error(ier)
          END IF
       END DO
    END DO
  
    !C PRINT *,'Spline_2D_point: cleaning up'
  
    CALL Ezspline_free(spline_o, ier)
    CALL EZspline_error(ier)
  
  
    RETURN
  
  END SUBROUTINE Spline_2D_periodic
  
  
!==============================================================================
  SUBROUTINE Spline_2D_point(x_1, x_2, f, x, y, func, ierrDomain)
  
    USE EZspline_obj ! import the modules
    USE EZspline  
    use nrtype, ONLY : DP
  
    IMPLICIT NONE
  
    INTEGER, PARAMETER :: r8 = DP 
    REAL(r8), PARAMETER :: twopi = 6.2831853071795862320_r8
  
    REAL(r8), DIMENSION(:), INTENT(IN) :: x_1, x_2 ! independent variable
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(r8), DIMENSION(:,:), INTENT(OUT) :: func   ! interpolated values 
    INTEGER           :: n1, n2, n_x, n_y, ier, BCS1(2), BCS2(2)
    TYPE(EZspline2_r8) :: spline_o ! 2-D EZspline object
  
    INTEGER           :: j, k
    INTEGER, INTENT(OUT) :: ierrDomain
    REAL(r8), DIMENSION(:,:), INTENT(IN)    :: x, y  ! grid of points for output
  
    n1 = SIZE(x_1)
    n2 = SIZE(x_2)
  
    ! Boundary conditions for interpolation
    BCS1 = (/ 0, 0 /) ! "not-a-knot" spline in first coordinate
    !C BCS1 = (/ 2, 2 /)    ! "natural" spline 
  
    BCS2 = (/-1, -1/)  ! periodic spline in second coordinate
    ! BCS2 = (/ 0, 0 /) ! "not-a-knot" spline in second coordinate
    ! BCS2 = (/ 2, 2 /) !  "natural" spline
  
    ! initialize/allocate memory
    !  PRINT *,'initializing...'
  
    ! CALL EZspline_init(spline_o, n1, n2, BCS1, BCS2, ier)
    CALL EZLinear_init(spline_o, n1, n2, ier)
    CALL EZspline_error(ier)
  
    spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, 
    ! and is aliased to x_1 now
    spline_o%x2 = x_2
    spline_o%isHermite = 1 ! Akima spline; smoother, more "natural"  interpolation (see Akima's paper)
  
      ! need to set explicitly the following if boundary conditions 
    ! are anything but not-a-knot or periodic, i.e. BCS(n) has
    ! element /= -1, 0.
  !   spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
  !  spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs
  
  
    !  PRINT *,'setting up 2D spline ...'
  
    CALL EZspline_setup(spline_o, f, ier)
    CALL EZspline_error(ier)
  
    ! save object
  
    !C CALL EZspline_save(spline_o, "spline.nc", ier)  ! Not necessary to save it if not loading again; slow by NFS on hyrax nodes
    CALL EZspline_error(ier)
  
    ! Cloud interpolation
  
    n_x = SIZE(x,1)
    n_y = SIZE(x,2)   ! x and y have the same shape and size
  
    DO j = 1, n_x
       DO k = 1, n_y
          CALL EZspline_isInDomain(spline_o, x(j,k), y(j,k), ierrDomain)
          IF (ierrDomain <= 0) THEN
             CALL EZspline_interp(spline_o, x(j,k), y(j,k), func(j,k), ier)
             CALL EZspline_error(ier)
          END IF
       END DO
    END DO
  
    !C PRINT *,'Spline_2D_point: cleaning up'
  
    CALL Ezspline_free(spline_o, ier)
    CALL EZspline_error(ier)

    RETURN
  
  END SUBROUTINE Spline_2D_point
  
!==========================================
  SUBROUTINE Spline_coord_2nd_derivs(x1, x2, x3, f, derivsX1, derivsX2, derivsX3)
  
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
  
    REAL(r8), DIMENSION(:), INTENT(IN) :: x1, x2, x3 ! independent variables
    REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f
    REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: derivsX1, derivsX2, derivsX3 ! function values
    INTEGER n1, n2, n3, ier, BCS1(2), BCS2(2), BCS3(2)
    TYPE(EZspline3_r8) :: spline_o ! 3-D EZspline object
  
    INTEGER k1, k2, k3
    REAL(r8), DIMENSION(:), ALLOCATABLE :: z1, z2, z3
    REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fz
  
    n1 = SIZE(x1)
    n2 = SIZE(x2)
    n3 = SIZE(x3)
  
    !   ALLOCATE(x1(n1), x2(n2), x3(n3), f(n1,n2,n3))
  
    BCS1 = (/ 2, 2 /)    ! "natural" spline in theta
    BCS2 = (/ 2, 2 /)    ! "natural" spline in psi
    BCS3 = (/ -1, -1 /)  ! periodic in zeta
  
    ! initialize/allocate memory
    PRINT *,'initializing...'
  
    CALL EZspline_init(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
    CALL EZspline_error(ier)
  
    spline_o%x1 = x1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, and is aliased to x1 now
    spline_o%x2 = x2    ! necessary if spline_o%x2 not in [0, 2 pi]
    ! set explicitely spline_o%x3
    ! spline_o%x3 = -1._r8 + 2._r8*(/ ( (REAL(i-1,r8)/REAL(n3-1,r8))**2,  i=1, n3 ) /)
    spline_o%x3 = x3
  
    ! need to set explicitly the following if boundary conditions 
    ! are anything but not-a-knot or periodic, i.e. BCS(n) has
    ! element /= -1, 0.
    spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
    spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs
    ! spline_o%bcval3min = ... ; spline_o%bcval3max = ...
  
  
    ! compute cubic spline coefficients
  
    !  DO i3=1, n3
    !     DO i2 = 1, n2
    !        DO i1 = 1, n1
    !           f(i1,i2,i3) = (COS(twopi*i1/(n1-1))* &
    !                & SIN(twopi*i2/(n2-1))**2)*EXP(spline_o%x3(i3))
    !        ENDDO
    !     ENDDO
    !  ENDDO
  
    PRINT *,'setting up...'
  
    CALL EZspline_setup(spline_o, f, ier)
    CALL EZspline_error(ier)
  
    ! save object
  
    CALL EZspline_save(spline_o, "spline.nc", ier)
    CALL EZspline_error(ier)
    
    ! point interpolation
  
    ! p1 = twopi/2._r8
    ! p2 = 0._r8
    ! p3 = -0.5_r8
    ! PRINT *,'point interpolation...'
  
    ! CALL EZspline_interp(spline_o, p1, p2, p3, fp, ier)
    ! CALL EZspline_error(ier)
  
    ! cloud interpolation
  
    ! m  = 21
    ! ALLOCATE(y1(m), y2(m), y3(m), fy(m))
  
    ! y1 = spline_o%x1min + (spline_o%x1max-spline_o%x1min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
    ! y2 = spline_o%x2min + (spline_o%x2max-spline_o%x2min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
    ! y3 = spline_o%x3min + (spline_o%x3max-spline_o%x3min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
  
    ! PRINT *,'cloud interpolation...'
  
    ! CALL EZspline_interp(spline_o, m, y1, y2, y3, fy, ier)
    ! CALL EZspline_error(ier)

    ! Array interpolation
  
    k1 = n1
    k2 = n2
    k3 = n3
  
    ALLOCATE(z1(k1), z2(k2), z3(k3), fz(k1,k2,k3))
  
    z1 = x1 ! spline_o%x1min + (spline_o%x1max-spline_o%x1min)* &
    ! & (/ ( REAL(i-1,r8)/REAL(k1-1,r8),  i=1, k1 ) /)
    z2 = x2 !spline_o%x2min + (spline_o%x2max-spline_o%x2min)* &
    ! & (/ ( REAL(i-1,r8)/REAL(k2-1,r8),  i=1, k2 ) /)
    z3 = x3 ! spline_o%x3min + (spline_o%x3max-spline_o%x3min)* &
    ! & (/ ( REAL(i-1,r8)/REAL(k3-1,r8),  i=1, k3 ) /)
  
  !  PRINT *,'array interpolation...'
  !  CALL EZspline_interp(spline_o, k1, k2, k3, z1, z2, z3, fz, ier)
  !  CALL EZspline_error(ier)
  
    PRINT*, 'derivatives...'
  
    CALL EZspline_derivative(spline_o, 2, 0, 0, k1, k2, k3, z1, z2, z3, derivsX1, ier)
    ! partial wrt theta for all points
  
    CALL EZspline_derivative(spline_o, 0, 2, 0, k1, k2, k3, z1, z2, z3, derivsX2, ier)
    ! partial wrt psi for all points
  
    CALL EZspline_derivative(spline_o, 0, 0, 2, k1, k2, k3, z1, z2, z3, derivsX3, ier)
    ! partial wrt zeta for all points
  
    ! clean up and free up memory
  
    PRINT *,'cleaning up'
  
    CALL Ezspline_free(spline_o, ier)
    CALL EZspline_error(ier)
  
    !  DEALLOCATE(y1, y2, y3, fy)
    DEALLOCATE(z1, z2, z3, fz)
  
    RETURN
  END SUBROUTINE Spline_coord_2nd_derivs
  
  
!==============================================================================
  SUBROUTINE Spline_coord_derivs(x_1, x_2, x_3, f_input, derivsX1, derivsX2, derivsX3)
  
    USE EZspline_obj ! import the modules
    USE EZspline  
    use nrtype, ONLY: DP
  
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
   !BCS1 = (/ 2, 2 /) ! "natural" spline in theta
    BCS2 = (/ 0, 0 /) ! "not-a-knot" spline in rho
   !BCS2 = (/ 2, 2 /) ! "natural" spline in rho
   !BCS2 = (/ 1, 1 /) ! first derivative imposed in rho
    BCS3 = (/ -1,-1 /)  ! periodic in zeta
   !BCS3 = (/ 0, 0 /)
  
    ! initialize/allocate memory
    ! PRINT *,'initializing...'
  
    !CALL EZspline_init(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
    CALL EZLinear_init(spline_o, n1, n2, n3, ier)
    CALL EZspline_error(ier)
  
    spline_o%x1 = x_1    ! necessary if spline_o%x1 not in [0, 2 pi]
    spline_o%x2 = x_2    ! necessary if spline_o%x2 not in [0, 2 pi]
  
    ! set explicitely spline_o%x3
  
    spline_o%x3 = x_3
  
    spline_o%isHermite = 1 ! Akima spline
  
  
    ! need to set explicitly the following if boundary conditions 
    ! are anything but not-a-knot or periodic, i.e. BCS(n) has element /= -1, 0.
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
  
!==========================================
  SUBROUTINE Spline_coord_derivs_large(x1, x2, x3, f, func, deriv1, deriv2)
  
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
  
    REAL(r8), DIMENSION(:), INTENT(IN) :: x1, x2, x3 ! independent variables
    REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: f
    REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: func, deriv1, deriv2 ! function values
    INTEGER n1, n2, n3, ier, BCS1(2), BCS2(2), BCS3(2), j
    TYPE(EZspline3_r8) :: spline_o ! 3-D EZspline object
  
    INTEGER k1, k2, k3
    REAL(r8), DIMENSION(:), ALLOCATABLE :: z1, z2, z3
  
    n1 = SIZE(x1)
    n2 = SIZE(x2)
    n3 = SIZE(x3)
  
    !   ALLOCATE(x1(n1), x2(n2), x3(n3), f(n1,n2,n3))
  
    BCS1 = (/ 2, 2 /)    ! "natural" spline in theta
      BCS2 = (/ 0, 0 /) ! "not-a-knot spline in rho
   ! BCS2 = (/ 2, 2 /)    ! "natural" spline in rho
    BCS3 = (/ -1, -1 /)  ! periodic in zeta
  
  print*, 'In spline large, size(f)', size(f,2)
  
    ! initialize/allocate memory
    PRINT *,'initializing...'
  
    CALL EZspline_init(spline_o, n1, n2, n3, BCS1, BCS2, BCS3, ier)
    CALL EZspline_error(ier)
  
    spline_o%x1 = x1    ! necessary if spline_o%x1 not in [0, 2 pi]; spline_o%x1 is a pointer, and is aliased to x1 now
    spline_o%x2 = x2    ! necessary if spline_o%x2 not in [0, 2 pi]
    ! set explicitely spline_o%x3
    ! spline_o%x3 = -1._r8 + 2._r8*(/ ( (REAL(i-1,r8)/REAL(n3-1,r8))**2,  i=1, n3 ) /)
    spline_o%x3 = x3
  
    ! need to set explicitly the following if boundary conditions 
    ! are anything but not-a-knot or periodic, i.e. BCS(n) has
    ! element /= -1, 0.
    spline_o%bcval1min = 0._r8 ; spline_o%bcval1max = 0._r8  ! values for 2nd derivatives at BCs
   ! spline_o%bcval2min = 0._r8 ; spline_o%bcval2max = 0._r8  ! values for 2nd derivatives at BCs
    ! spline_o%bcval3min = ... ; spline_o%bcval3max = ...
  
  
    ! compute cubic spline coefficients
  
    !  DO i3=1, n3
    !     DO i2 = 1, n2
    !        DO i1 = 1, n1
    !           f(i1,i2,i3) = (COS(twopi*i1/(n1-1))* &
    !                & SIN(twopi*i2/(n2-1))**2)*EXP(spline_o%x3(i3))
    !        ENDDO
    !     ENDDO
    !  ENDDO
  
    PRINT *,'setting up...'
  
    CALL EZspline_setup(spline_o, f, ier)
    CALL EZspline_error(ier)
  
    ! save object
  
    CALL EZspline_save(spline_o, "spline.nc", ier)
    CALL EZspline_error(ier)
  
  
    ! point interpolation
  
    ! p1 = twopi/2._r8
    ! p2 = 0._r8
    ! p3 = -0.5_r8
    ! PRINT *,'point interpolation...'
  
    !   CALL EZspline_interp(spline_o, p1, p2, p3, fp, ier)
    !   CALL EZspline_error(ier)
  
    ! cloud interpolation
  
    ! m  = 21
    ! ALLOCATE(y1(m), y2(m), y3(m), fy(m))
  
    ! y1 = spline_o%x1min + (spline_o%x1max-spline_o%x1min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
    ! y2 = spline_o%x2min + (spline_o%x2max-spline_o%x2min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
    ! y3 = spline_o%x3min + (spline_o%x3max-spline_o%x3min)* &
    !      & (/ ( REAL(i-1,r8)/REAL(m-1,r8),  i=1, m ) /)
  
  
    ! PRINT *,'cloud interpolation...'
  
    ! CALL EZspline_interp(spline_o, m, y1, y2, y3, fy, ier)
    ! CALL EZspline_error(ier)
  
    ! Array interpolation
  
    k1 = 1  ! On Eq. plane
    k2 = SIZE(func, 2)
    print*, 'k2 =', k2
    k3 = 1  ! At midnight
  
    ALLOCATE(z1(k1), z2(k2), z3(k3))
  
    z1(1) = twopi / 4._r8 ! spline_o%x1min + (spline_o%x1max-spline_o%x1min)* &
    ! & (/ ( REAL(i-1,r8)/REAL(k1-1,r8),  i=1, k1 ) /)
  
    DO j = 1, k2
       z2(j) = (float(j) - 1._r8) / (float(k2) - 1._r8)
    END DO
  
    ! need the midnight value here for z3
    z3(1) = 19_r8 * twopi / 35_r8 ! spline_o%x3min + (spline_o%x3max-spline_o%x3min)* &
    ! & (/ ( REAL(i-1,r8)/REAL(k3-1,r8),  i=1, k3 ) /)
  
    PRINT *,'array interpolation...'
    CALL EZspline_interp(spline_o, k1, k2, k3, z1, z2, z3, func, ier)
    CALL EZspline_error(ier)
  
    PRINT*, 'derivatives...'
  
  !  CALL EZspline_derivative(spline_o, 1, 0, 0, k1, k2, k3, z1, z2, z3, derivsX1, ier)
    ! partial wrt theta for all points
  
    CALL EZspline_derivative(spline_o, 0, 1, 0, k1, k2, k3, z1, z2, z3, deriv1, ier)
    ! partial wrt psi for all points
  
    CALL EZspline_derivative(spline_o, 0, 2, 0, k1, k2, k3, z1, z2, z3, deriv2, ier)
  
  !  CALL EZspline_derivative(spline_o, 0, 0, 1, k1, k2, k3, z1, z2, z3, derivsX3, ier)
    ! partial wrt zeta for all points
  
    ! clean up and free up memory
  
    PRINT *,'cleaning up'
  
    CALL Ezspline_free(spline_o, ier)
    CALL EZspline_error(ier)
  
    !  DEALLOCATE(y1, y2, y3, fy)
    DEALLOCATE(z1, z2, z3)
  
    RETURN
  END SUBROUTINE Spline_coord_derivs_large
  
  
!==========================================
  SUBROUTINE Spline_derivs_1D(x1, f, func, derivsX1)
  
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
  
    BCS1 = (/ 0, 0 /) ! "not-a-knot" spline 
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
  
    !C CALL EZspline_save(spline_o, "spline.nc", ier) ! Not necessary, it is defined at each call
    CALL EZspline_error(ier)
  
    ! Array interpolation  
    k1 = SIZE(derivsX1) 
    ALLOCATE(z1(k1))
  
    DO j = 1, k1
       z1(j) = (float(j) - 1._r8) / (float(k1) - 1._r8)
    END DO
  
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
  END SUBROUTINE Spline_derivs_1D

!==============================================================================
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
  
!==============================================================================
  SUBROUTINE spline(x,y,yp1,ypn,y2)   ! MOD everything is DP
    USE nrtype
    USE nrutil, ONLY: assert_eq
    USE nrmod,  ONLY: tridag
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN)          ::    x,y
    REAL(DP), INTENT(IN)                        ::    yp1,ypn
    REAL(DP), DIMENSION(:), INTENT(OUT)         ::    y2
    INTEGER(I4B)                                ::    n
    REAL(DP), DIMENSION(size(x))                ::    a,b,c,r
  
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_dp) then
       r(1)=0.0
       c(1)=0.0
    else
       r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5
    end if
    if (ypn > 0.99e30_dp) then
       r(n)=0.0
       a(n)=0.0
    else
       r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
  END SUBROUTINE spline
  
!==============================================================================
  FUNCTION splint(xa,ya,y2a,x)
    USE nrtype
    USE nrutil, ONLY: assert_eq,nrerror
    USE nrmod,  ONLY: locate
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: splint
    INTEGER(I4B) :: khi,klo,n
    REAL(DP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
  END FUNCTION splint
  
END MODULE ModScbSpline
