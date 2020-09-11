!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamGSL

  implicit none

  ! C Interface Bindings
  interface
     subroutine RamGSL_Initialization_c(a) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: a
     end subroutine RamGSL_Initialization_c

     subroutine Interpolation_Smooth_c(i1,i2,xa,fa,xb,fb,err) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: i1, i2
        integer(c_int) :: err(*)
        real(c_double) :: xa(*), fa(*)
        real(c_double) :: xb(*), fb(*)
     end subroutine Interpolation_Smooth_c

     subroutine Interpolation_1D_c(n,i1,i2,xa,fa,xb,fb,err) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n, i1, i2
        integer(c_int) :: err(*)
        real(c_double) :: xa(*), fa(*)
        real(c_double) :: xb(*), fb(*)
     end subroutine Interpolation_1D_c

     subroutine Interpolation_2D_c(i1,j1,i2,j2,xa,ya,fa,xb,yb,fb,err) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: i1, j1, i2, j2
        integer(c_int) :: err(*)
        real(c_double) :: xa(*), ya(*), fa(*)
        real(c_double) :: xb(*), yb(*), fb(*)
     end subroutine Interpolation_2D_c

     subroutine Interpolation_3D_c(i1,j1,k1,i2,j2,k2,xa,ya,za,fa,xb,yb,zb,fb,err) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: i1, j1, k1, i2, j2, k2
        integer(c_int) :: err(*)
        real(c_double) :: xa(*), ya(*), za(*), fa(*)
        real(c_double) :: xb(*), yb(*), zb(*), fb(*)
     end subroutine Interpolation_3D_c

     subroutine Interpolation_Derivs_c(n,i1,xa,fa,dx,err) BIND(C)
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n, i1
        integer(c_int) :: err(*)
        real(c_double) :: xa(*), fa(*)
        real(c_double) :: dx(*)
     end subroutine Interpolation_Derivs_c

     subroutine Integrator_c(nT,nPa,mirror,theta,field,yI,yH) BIND(C)
         use, intrinsic :: iso_c_binding
        integer(c_int), value :: nT, nPa
        real(c_double) :: theta(*), field(*) ! These have dimension nT
        real(c_double) :: mirror(*), yI(*), yH(*) ! These have dimension nPa
     end subroutine Integrator_c

     subroutine BounceAverage_c(nT,nPa,mirror,theta,field,variable,yV) BIND(C)
         use, intrinsic :: iso_c_binding
        integer(c_int), value :: nT, nPa
        real(c_double) :: theta(*), field(*), variable(*) ! These have dimension nT
        real(c_double) :: mirror(*), yV(*) ! These have dimension nPa
     end subroutine BounceAverage_c

  end interface

  ! Spline Interpolation interfaces
  !! 1D Interpolations
  interface GSL_Interpolation_1D
     module procedure Interpolation_1D_array, Interpolation_1D_point
  end interface
  !! 2D Interpolations
  interface GSL_Interpolation_2D
     module procedure Interpolation_2D_array, Interpolation_2D_cloud, &
                      Interpolation_2D_NN_point, Interpolation_2D_point, &
                      Interpolation_2D_NN_array
  end interface

  ! Spline Derivative interfaces
  interface GSL_Derivs
     module procedure Interpolation_3D_Derivs, Interpolation_2D_Derivs, &
                      Interpolation_1D_Derivs
  end interface

  ! Nearest Neigbor Interfaces
  interface GSL_NN
    module procedure NN_Interpolation_2D, NN_Interpolation_3D
  end interface

  contains
!==================================================================================================
!  subroutine GSL_error(modName, subName, message, errno)
!
!    implicit none
!
!    character(len=*) :: modName, subName, message
!    integer :: errno
!
!    return
!  endsubroutine GSL_error
!==================================================================================================
  subroutine GSL_Initialize
    ! Initialized GSL error checking and other c constructs
    use, intrinsic :: iso_c_binding


    implicit none

    integer(c_int) :: a

    a = 1
    call RamGSL_Initialization_c(a)

    return

  end subroutine GSL_Initialize

!==================================================================================================
  subroutine GSL_BounceAverage(variable,bM,theta,field,y_v)
    ! Bounce average a variable

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP

    implicit none

    real(DP), INTENT(IN) :: theta(:), field(:), variable(:)
    real(DP), INTENT(INOUT) :: y_v(:), bM(:)

    integer(c_int), SAVE :: nT_c, nPa_c
    real(c_double), ALLOCATABLE, SAVE :: bM_c(:), y_v_c(:)
    real(c_double), ALLOCATABLE, SAVE :: theta_c(:), field_c(:), variable_c(:)
    !$OMP THREADPRIVATE(theta_C,field_C,variable_C,bM_c,y_v_c,nT_c,nPa_c)

    nT_c = size(theta,1)
    nPa_c = size(bM,1)

    allocate(theta_c(nT_c),field_c(nT_c),variable_c(nT_c))
    allocate(bM_c(nPa_c),y_v_c(nPa_c))

    bM_c = REAL(bM,c_double)
    theta_c = REAL(theta,c_double)
    field_c = REAL(field,c_double)
    variable_c = REAL(variable,c_double)

    call BounceAverage_c(nT_c,nPa_c,bM_c,theta_c,field_c,variable_c,y_v_c)

    y_v = REAL(y_v_c,DP)
    bM = REAL(bM_c,DP)

    deallocate(theta_c, field_c, variable_c, bM_c, y_v_c)

    return

  end subroutine GSL_BounceAverage

!==================================================================================================
  subroutine GSL_Integration_hI(bM, theta, field, y_i, y_h)
    ! Compute h and I integrals

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP

    implicit none

    real(DP), INTENT(IN) :: theta(:), field(:)
    real(DP), INTENT(INOUT) :: y_i(:), y_h(:), bM(:)

    integer(c_int), SAVE :: nT_c, nPa_c
    real(c_double), ALLOCATABLE, SAVE :: bM_c(:), y_i_c(:), y_h_c(:)
    real(c_double), ALLOCATABLE, SAVE :: theta_c(:), field_c(:)
    !$OMP THREADPRIVATE(theta_C,field_C,bM_c,y_i_c,y_h_c,nT_c,nPa_c)

    nT_c = size(theta,1)
    nPa_c = size(bM,1)

    allocate(theta_c(nT_c),field_c(nT_c))
    allocate(bM_c(nPa_c),y_i_c(nPa_c),y_h_c(nPa_c))

    bM_c = REAL(bM,c_double)
    theta_c = REAL(theta,c_double)
    field_c = REAL(field,c_double)

    call Integrator_c(nT_c,nPa_c,bM_c,theta_c,field_c,y_i_c,y_h_c)

    y_i = REAL(y_i_c,DP)
    y_h = REAL(y_h_c,DP)
    bM = REAL(bM_c,DP)

    deallocate(theta_c, field_c, bM_c, y_i_c, y_h_c)

    return

  end subroutine GSL_Integration_hI

!==================================================================================================
  subroutine GSL_Smooth_1D(x1,f1,x2,f2,err)
    ! 1D interpolation using B-Splines

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, intent(out) :: err
    real(DP), DIMENSION(:), INTENT(IN)  :: x1, f1
    real(DP), DIMENSION(:), INTENT(INOUT) :: x2, f2

    INTEGER(c_int) :: n1, n2, err_c(1)
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, fa, xb, fb

    n1 = size(x1,1)
    n2 = size(x2,1)

    ALLOCATE(xa(n1),fa(n1),xb(n2),fb(n2))

    xa = REAL(x1,c_double)
    fa = REAL(f1,c_double)
    xb = REAL(x2,c_double)

    call Interpolation_Smooth_c(n1,n2,xa,fa,xb,fb,err_c)

    x2 = REAL(xb,DP)
    f2 = REAL(fb,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,fa,xb,fb)

    RETURN

  end subroutine GSL_Smooth_1D

!==================================================================================================
  subroutine Interpolation_1D_array(x1,f1,x2,f2,err,ctype)
    ! Subroutine for performing 1D interpolations using the GSL package.
    ! The default interpolation uses the Steffen spline. If ctype is included in
    ! the call you can specify Cubic, Akima, or Linear splines be used instead.

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP

    implicit none

    integer :: i,j,n0
    integer, intent(out) :: err
    real(DP), DIMENSION(:), INTENT(IN)  :: x1, f1, x2
    real(DP), DIMENSION(:), INTENT(INOUT) :: f2
    real(DP), DIMENSION(:), ALLOCATABLE :: y1, g1
    character(len=*), optional, INTENT(IN) :: ctype

    INTEGER(c_int) :: n1, n2, err_c(1), ntype
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, fa, xb, fb

    ! Check for monotonicity in x1
    n0 = size(x1,1)
    allocate(y1(n0), g1(n0))
    y1(1) = x1(1)
    g1(1) = f1(1)
    j = 1
    do i = 2, n0
       if (x1(i) > y1(j)) then
          j = j + 1
          y1(j) = x1(i)
          g1(j) = f1(i)
       endif
    enddo

    n1 = j
    n2 = size(x2,1)

    ALLOCATE(xa(n1),fa(n1),xb(n2),fb(n2))

    xa = REAL(y1(1:j),c_double)
    fa = REAL(g1(1:j),c_double)
    xb = REAL(x2,c_double)

    if (present(ctype)) then
       if (ctype.eq.'Cubic') then
          ntype = 1
       elseif (ctype.eq.'Akima') then
          ntype = 2
       elseif (ctype.eq.'Steffen') then
          ntype = 3
       elseif (ctype.eq.'Linear') then
          ntype = 0
       elseif (ctype.eq.'Smooth') then
          ntype = 4
       else
          ntype = 3
       endif
    else
       ntype = 3
    endif

    call Interpolation_1D_c(ntype,n1,n2,xa,fa,xb,fb,err_c)

    f2 = REAL(fb,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,fa,xb,fb)

    RETURN

  end subroutine Interpolation_1D_array

!==================================================================================================
  subroutine Interpolation_1D_point(x1,f1,x2,f2,err,ctype)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP

    implicit none

    integer, INTENT(OUT)  :: err
    real(DP), INTENT(IN)  :: x1(:), f1(:), x2
    real(DP), INTENT(INOUT) :: f2
    character(len=*), optional, INTENT(IN) :: ctype

    INTEGER(c_int) :: n1, n2, err_c(1), ntype
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, fa, xb, fb

    n1 = size(x1,1)
    n2 = 1

    ALLOCATE(xa(n1),fa(n1),xb(n2),fb(n2))

    xa = REAL(x1,c_double)
    fa = REAL(f1,c_double)
    xb(1) = REAL(x2,c_double)

    if (present(ctype)) then
       if (ctype.eq.'Cubic') then
          ntype = 1
       elseif (ctype.eq.'Akima') then
          ntype = 2
       elseif (ctype.eq.'Steffen') then
          ntype = 3
       elseif (ctype.eq.'Linear') then
          ntype = 0
       elseif (ctype.eq.'Smooth') then
          ntype = 4
       else
          ntype = 3
       endif
    else
       ntype = 3
    endif

    call Interpolation_1D_c(ntype,n1,n2,xa,fa,xb,fb,err_c)

    f2 = REAL(fb(1),DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,fa,xb,fb)

    RETURN

  end subroutine Interpolation_1D_point

!==================================================================================================
  subroutine Interpolation_2D_NN_point(x1,y1,f1,x2,y2,f2,err)

    use nrtype, ONLY: DP


    implicit none

    integer, parameter :: NN = 9
    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:,:), INTENT(IN) :: x1, y1, f1
    real(DP), INTENT(IN)  :: x2, y2
    real(DP), INTENT(INOUT) :: f2

    integer :: i,j,k, n1, m1, nTotal, iTotal, iTemp(1), stat
    real(DP) :: xNear(NN), yNear(NN), fNear(NN)
    real(DP), ALLOCATABLE :: distance(:), xScatter(:), yScatter(:), fScatter(:)

    n1 = SIZE(x1,1)
    m1 = SIZE(x1,2)

    nTotal = n1*m1
    ALLOCATE(distance(nTotal))
    ALLOCATE(xScatter(nTotal))
    ALLOCATE(yScatter(nTotal))
    ALLOCATE(fScatter(nTotal))

    iTotal = 0
    DO i = 1,n1
       DO j = 1,m1
          iTotal = iTotal + 1
          xScatter(iTotal) = x1(i,j)
          yScatter(iTotal) = y1(i,j)
          fScatter(iTotal) = f1(i,j)
       END DO
    END DO

    ! Use a nearest neighbor search to interpolate
    distance = (xScatter - x2)**2 &
              +(yScatter - y2)**2
    do k = 1,NN
       iTemp = minloc(distance)
       xNear(k) = xScatter(iTemp(1))
       yNear(k) = yScatter(iTemp(1))
       fNear(k) = fScatter(iTemp(1))
       distance(iTemp(1)) = 999999.9
    end do
    CALL NN_Interpolation_2D(xNear,yNear,fNear,x2,y2,f2,stat)

    err = stat

    DEALLOCATE(distance,xScatter,yScatter,fScatter)

    RETURN

  end subroutine Interpolation_2D_NN_point

!==================================================================================================
  subroutine Interpolation_2D_NN_array(x1,y1,f1,x2,y2,f2,err)

    use nrtype, ONLY: DP


    implicit none

    integer, parameter :: NN = 9
    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:,:), INTENT(IN)  :: x1, y1, f1, x2, y2
    real(DP), DIMENSION(:,:), INTENT(INOUT) :: f2

    integer :: i,j,k, n1, m1, n2, m2, nTotal, iTotal, iTemp(1), stat
    real(DP) :: xo, yo, xNear(NN), yNear(NN), fNear(NN)
    real(DP), ALLOCATABLE :: distance(:), xScatter(:), yScatter(:), fScatter(:)

    n1 = SIZE(x1,1)
    m1 = SIZE(x1,2)
    n2 = SIZE(x2,1)
    m2 = SIZE(x2,2)

    nTotal = n1*m1
    ALLOCATE(distance(nTotal))
    ALLOCATE(xScatter(nTotal))
    ALLOCATE(yScatter(nTotal))
    ALLOCATE(fScatter(nTotal))

    iTotal = 0
    DO i = 1,n1
       DO j = 1,m1
          iTotal = iTotal + 1
          xScatter(iTotal) = x1(i,j)
          yScatter(iTotal) = y1(i,j)
          fScatter(iTotal) = f1(i,j)
       END DO
    END DO

    ! Use a nearest neighbor search to interpolate
    DO i = 1,n2
       DO j = 1,m2
          xo = x2(i,j)
          yo = y2(i,j)
          distance = (xScatter - xo)**2 &
                    +(yScatter - yo)**2
          do k = 1,NN
             iTemp = minloc(distance)
             xNear(k) = xScatter(iTemp(1))
             yNear(k) = yScatter(iTemp(1))
             fNear(k) = fScatter(iTemp(1))
             distance(iTemp(1)) = 999999.9
          end do
          CALL NN_Interpolation_2D(xNear,yNear,fNear,xo,yo,f2(i,j),stat)
       ENDDO
    ENDDO

    err = stat

    DEALLOCATE(distance,xScatter,yScatter,fScatter)

    RETURN

  end subroutine Interpolation_2D_NN_array

!==================================================================================================
  subroutine Interpolation_2D_point(x1,y1,f1,x2,y2,f2,err)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), INTENT(IN) :: x1(:), y1(:), x2, y2
    real(DP), INTENT(IN) :: f1(:,:)
    real(DP), INTENT(INOUT) :: f2
    integer :: i, j

    INTEGER(c_int) :: n1, m1, n2, m2, err_c(1)
    real(c_double), DIMENSION(:), ALLOCATABLE   :: xa, ya
    real(c_double), DIMENSION(:,:), ALLOCATABLE :: xb, yb, fa, fb

    n1 = size(x1,1)
    m1 = size(y1,1)
    n2 = 1
    m2 = 1

    ALLOCATE(xa(n1),ya(m1),fa(n1,m1),xb(n2,m2),yb(n2,m2),fb(n2,m2))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    fa = REAL(f1,c_double)
    do j = 1,m2
       xb(:,j) = REAL(x2,c_double)
    enddo
    do i = 1,n2
       yb(i,:) = REAL(y2,c_double)
    enddo

    call Interpolation_2D_c(n1,m1,n2,m2,xa,ya,fa,xb,yb,fb,err_c)

    f2 = REAL(fb(1,1),DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,ya,fa,xb,yb,fb)

    RETURN

  end subroutine Interpolation_2D_point

!==================================================================================================
  subroutine Interpolation_2D_array(x1,y1,f1,x2,y2,f2,err)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN)  :: x1, y1, x2, y2
    real(DP), DIMENSION(:,:), INTENT(IN) :: f1
    real(DP), DIMENSION(:,:), INTENT(INOUT) :: f2
    integer :: i, j

    INTEGER(c_int) :: n1, m1, n2, m2, err_c(1)
    real(c_double), DIMENSION(:), ALLOCATABLE   :: xa, ya
    real(c_double), DIMENSION(:,:), ALLOCATABLE :: xb, yb, fa, fb

    n1 = size(x1,1)
    m1 = size(y1,1)
    n2 = size(x2,1)
    m2 = size(y2,1)

    ALLOCATE(xa(n1),ya(m1),fa(n1,m1),xb(n2,m2),yb(n2,m2),fb(n2,m2))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    fa = REAL(f1,c_double)
    do j = 1,m2
       xb(:,j) = REAL(x2,c_double)
    enddo
    do i = 1,n2
       yb(i,:) = REAL(y2,c_double)
    enddo

    call Interpolation_2D_c(n1,m1,n2,m2,xa,ya,fa,xb,yb,fb,err_c)

    f2 = REAL(fb,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,ya,fa,xb,yb,fb)

    RETURN

  end subroutine Interpolation_2D_array

!==================================================================================================
  subroutine Interpolation_2D_cloud(x1,y1,f1,x2,y2,f2,err)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN)  :: x1, y1
    real(DP), DIMENSION(:,:), INTENT(IN) :: f1, x2, y2
    real(DP), DIMENSION(:,:), INTENT(INOUT) :: f2

    INTEGER(c_int) :: n1, m1, n2, m2, err_c(1)
    real(c_double), DIMENSION(:), ALLOCATABLE   :: xa, ya
    real(c_double), DIMENSION(:,:), ALLOCATABLE :: xb, yb, fa, fb

    n1 = size(x1,1)
    m1 = size(y1,1)
    n2 = size(x2,1)
    m2 = size(x2,2)

    ALLOCATE(xa(n1),ya(m1),fa(n1,m1),xb(n2,m2),yb(n2,m2),fb(n2,m2))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    fa = REAL(f1,c_double)
    xb = REAL(x2,c_double)
    yb = REAL(y2,c_double)

    call Interpolation_2D_c(n1,m1,n2,m2,xa,ya,fa,xb,yb,fb,err_c)

    f2 = REAL(fb,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,ya,fa,xb,yb,fb)

    RETURN

  end subroutine Interpolation_2D_cloud

!==================================================================================================
  subroutine Interpolation_3D(x1,y1,z1,f1,x2,y2,z2,f2,err)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN) :: x1, y1, z1, x2, y2, z2
    real(DP), DIMENSION(:,:,:), INTENT(IN)  :: f1
    real(DP), DIMENSION(:,:,:), INTENT(INOUT) :: f2

    INTEGER(c_int) :: n1, m1, l1, n2, m2, l2, err_c(1)
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, ya, za, xb, yb, zb
    real(c_double), DIMENSION(:,:,:), ALLOCATABLE :: fa, fb

    write(*,*) '3D Interpolation is not currently working, do you really need it?'
    return

    n1 = size(x1,1)
    m1 = size(y1,1)
    l1 = size(z1,1)
    n2 = size(x2,1)
    m2 = size(y2,1)
    l2 = size(z2,1)

    ALLOCATE(xa(n1),ya(m1),za(l1),fa(n1,m1,l1), &
             xb(n2),yb(m2),zb(l2),fb(n2,m2,l1))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    za = REAL(z1,c_double)
    fa = REAL(f1,c_double)
    xb = REAL(x2,c_double)
    yb = REAL(y2,c_double)
    zb = REAL(z2,c_double)

    call Interpolation_3D_c(n1,m1,l1,n2,m2,l2,xa,ya,za,fa,xb,yb,zb,fb,err_c)

    f2 = REAL(fb,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,ya,za,fa,xb,yb,zb,fb)

    RETURN

  end subroutine Interpolation_3D

!==================================================================================================
  subroutine Interpolation_1D_Derivs(x1,f1,dfdx,err,ctype)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN)  :: x1
    real(DP), DIMENSION(:), INTENT(IN)  :: f1
    real(DP), DIMENSION(:), INTENT(INOUT) :: dfdx
    character(len=*), optional, INTENT(IN) :: ctype

    INTEGER(c_int) :: n1, err_c(1), ntype
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa
    real(c_double), DIMENSION(:), ALLOCATABLE :: fa, dx

    n1 = size(x1,1)

    ALLOCATE(xa(n1),fa(n1),dx(n1))

    xa = REAL(x1,c_double)
    fa = REAL(f1,c_double)

    if (present(ctype)) then
       if (ctype.eq.'Cubic') then
          ntype = 1
       elseif (ctype.eq.'Akima') then
          ntype = 2
       elseif (ctype.eq.'Steffen') then
          ntype = 3
       elseif (ctype.eq.'Linear') then
          ntype = 0
       elseif (ctype.eq.'Smooth') then
          ntype = 4
       else
          ntype = 3
       endif
    else
       ntype = 3
    endif

    call Interpolation_Derivs_c(ntype,n1,xa,fa,dx,err_c)

    dfdx = REAL(dx,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,fa,dx)

    RETURN

  end subroutine Interpolation_1D_Derivs

!==================================================================================================
  subroutine Interpolation_2D_Derivs(x1,y1,f1,dfdx,dfdy,err,ctype)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN)    :: x1, y1
    real(DP), DIMENSION(:,:), INTENT(IN)  :: f1
    real(DP), DIMENSION(:,:), INTENT(INOUT) :: dfdx, dfdy
    character(len=*), optional, INTENT(IN) :: ctype
    integer :: i, j

    INTEGER(c_int) :: n1, m1, err_c(1), ntype
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, ya
    real(c_double), DIMENSION(:,:), ALLOCATABLE :: fa, dx, dy

    n1 = size(x1,1)
    m1 = size(y1,1)

    ALLOCATE(xa(n1),ya(m1),fa(n1,m1),dx(n1,m1),dy(n1,m1))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    fa = REAL(f1,c_double)

    if (present(ctype)) then
       if (ctype.eq.'Cubic') then
          ntype = 1
       elseif (ctype.eq.'Akima') then
          ntype = 2
       elseif (ctype.eq.'Steffen') then
          ntype = 3
       elseif (ctype.eq.'Linear') then
          ntype = 0
       elseif (ctype.eq.'Smooth') then
          ntype = 4
       else
          ntype = 3
       endif
    else
       ntype = 3
    endif

    do i = 1,n1
       call Interpolation_Derivs_c(ntype,m1,ya,fa(i,:),dy(i,:),err_c)
    end do
    do j = 1,m1
       call Interpolation_Derivs_c(ntype,n1,xa,fa(:,j),dx(:,j),err_c)
    end do

    dfdx = REAL(dx,DP)
    dfdy = REAL(dy,DP)
    err = int(err_c(1),kind=4)

    DEALLOCATE(xa,ya,fa,dx,dy)

    RETURN

  end subroutine Interpolation_2D_Derivs

!==================================================================================================
  subroutine Interpolation_3D_Derivs(x1,y1,z1,f1,dfdx,dfdy,dfdz,err,ctype)

    use, intrinsic :: iso_c_binding
    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), DIMENSION(:), INTENT(IN)      :: x1, y1, z1
    real(DP), DIMENSION(:,:,:), INTENT(IN)  :: f1
    real(DP), DIMENSION(:,:,:), INTENT(INOUT) :: dfdx, dfdy, dfdz
    character(len=*), optional, INTENT(IN) :: ctype
    integer :: i, j, k, e1, e2, e3

    INTEGER(c_int) :: n1, m1, l1, e1_c(1), e2_c(1), e3_c(1), ntype
    real(c_double), DIMENSION(:), ALLOCATABLE :: xa, ya, za
    real(c_double), DIMENSION(:,:,:), ALLOCATABLE :: fa, dx, dy, dz

    n1 = size(x1,1)
    m1 = size(y1,1)
    l1 = size(z1,1)

    ALLOCATE(xa(n1),ya(m1),za(l1),fa(n1,m1,l1), &
             dx(n1,m1,l1),dy(n1,m1,l1),dz(n1,m1,l1))

    xa = REAL(x1,c_double)
    ya = REAL(y1,c_double)
    za = REAL(z1,c_double)
    fa = REAL(f1,c_double)

    if (present(ctype)) then
       if (ctype.eq.'Cubic') then
          ntype = 1
       elseif (ctype.eq.'Akima') then
          ntype = 2
       elseif (ctype.eq.'Steffen') then
          ntype = 3
       elseif (ctype.eq.'Linear') then
          ntype = 0
       elseif (ctype.eq.'Smooth') then
          ntype = 4
       else
          ntype = 3
       endif
    else
       ntype = 3
    endif

    do i = 1,n1
       do j = 1,m1
          call Interpolation_Derivs_c(ntype,l1,za,fa(i,j,:),dz(i,j,:),e1_c)
       end do
       do k = 1,l1
          call Interpolation_Derivs_c(ntype,m1,ya,fa(i,:,k),dy(i,:,k),e2_c)
       end do
    end do
    do j = 1,m1
       do k = 1,l1
          call Interpolation_Derivs_c(ntype,n1,xa,fa(:,j,k),dx(:,j,k),e3_c)
       end do
    end do

    dfdx = REAL(dx,DP)
    dfdy = REAL(dy,DP)
    dfdz = REAL(dz,DP)
    e1 = int(e1_c(1),kind=4)
    e2 = int(e2_c(1),kind=4)
    e3 = int(e3_c(1),kind=4)
    err = abs(e1) + abs(e2) + abs(e3)

    DEALLOCATE(xa,ya,za,fa,dx,dy,dz)

    RETURN

  end subroutine Interpolation_3D_Derivs

!==================================================================================================
  subroutine NN_Interpolation_2D(x1,y1,f1,x2,y2,f2,err)

    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT) :: err
    real(DP), INTENT(IN)  :: x1(:), y1(:), f1(:), x2, y2
    real(DP), INTENT(INOUT) :: f2

    integer :: i,ii
    real(DP) :: d,wsum
    real(DP), ALLOCATABLE :: w(:)

    ii = size(x1,1)

    ALLOCATE(w(ii))

    ! Calculate NN weights
    wsum = 0._dp
    do i = 1,ii
       d = sqrt((x1(i)-x2)**2 + (y1(i)-y2)**2)
       if (abs(d).le.1e-9) then
          w = 0._dp
          w(i) = 1._dp
          wsum = 1._dp
          exit
       endif
       w(i) = 1/d**2
       wsum = wsum + w(i)
    enddo

    ! Calculated interpolated value
    f2 = 0._dp
    do i = 1,ii
       f2 = f2 + f1(i)*w(i)/wsum
    enddo

    err = 0

    DEALLOCATE(w)

    return

  end subroutine NN_Interpolation_2D

!==================================================================================================
  subroutine NN_Interpolation_3D(x1,y1,z1,f1,x2,y2,z2,f2,err)

    use nrtype, ONLY: DP


    implicit none

    integer, INTENT(OUT)  :: err
    real(DP), INTENT(IN)  :: x1(:), y1(:), z1(:), f1(:), x2, y2, z2
    real(DP), INTENT(INOUT) :: f2

    integer :: i,ii
    real(DP) :: d,wsum
    real(DP), ALLOCATABLE :: w(:)

    ii = size(x1,1)

    ALLOCATE(w(ii))

    ! Calculate NN weights
    wsum = 0._dp
    do i = 1,ii
       d = sqrt((x1(i)-x2)**2 + (y1(i)-y2)**2 + (z1(i)-z2)**2)
       if (abs(d).le.1e-9) then
          w = 0._dp
          w(i) = 1._dp
          wsum = 1._dp
          exit
       endif
       w(i) = 1/d**2
       wsum = wsum + w(i)
    enddo

    ! Calculated interpolated value
    f2 = 0._dp
    do i = 1,ii
       f2 = f2 + f1(i)*w(i)/wsum
    enddo

    err = 0

    DEALLOCATE(w)

    return

  end subroutine NN_Interpolation_3D

END MODULE ModRamGSL
