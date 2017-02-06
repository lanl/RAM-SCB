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
