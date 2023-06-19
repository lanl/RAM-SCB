!=============================================================================
module ModRamFunctions
!    Copyright (c) 2016, Los Alamos National Laboratory
!    All rights reserved.
!=============================================================================

  implicit none
  
  contains
!==============================================================================
  function RamFileName(PrefixIn, SuffixIn, TimeIn)
    ! Create a file name using the new RAM-SCB I/O filename standards:
    ! FileNameOut = PrefixIn_dYYYYMMDD_tHHMMSS.SuffixIn
    ! PrefixIn and SuffixIn are short, descriptive strings (e.g. 
    ! 'pressure', 'ram_o' for prefixes and 'dat', 'cdf' for suffixes.

    use ModTimeConvert


    implicit none

    character(len=200)           :: RamFileName
    type(TimeType),   intent(in) :: TimeIn
    character(len=*), intent(in) :: PrefixIn, SuffixIn
    !------------------------------------------------------------------------
    write(RamFileName, "(a,'_d',i4.4,2i2.2,'_t',3i2.2,'.',a)") &
         trim(PrefixIn), &
         TimeIn % iYear, &
         TimeIn % iMonth, &
         TimeIn % iDay, &
         TimeIn % iHour, &
         TimeIn % iMinute, &
         TimeIn % iSecond, &
         trim(SuffixIn)

  end function RamFileName
!===========================================================================
  subroutine get_ramdst(dstOut)
    ! Use simple DPS relationship to calculate Dst from RAM domain.
    ! Sums over all species, corrects for internal currents (factor of 1.3).

    use ModRamMain, ONLY: Real8_
    use ModRamGrids, ONLY: NR, NE, NT, NPA, nS
    use ModRamVariables, ONLY: f2, rfactor, upa, we, wmu, ekev, outsideMGNP

    implicit none

    real(kind=Real8_), intent(out) :: dstOut

    ! Factor2 includes conversions and factor of 1.3.
    real(kind=Real8_)           :: sumEnergy, factor2
    integer                     :: i, j, k, l, s
    !------------------------------------------------------------------------
    sumEnergy = 0.0
    dstOut = 0.0
    factor2 =-5.174E-30
    ! Sum energy over whole domain and all species.
    do s=1,nS; do i=2,nR; do k=2,nE; do l=1,nPa
       if(l.ge.uPa(i))cycle
       do j=1, nT-1
          if (outsideMGNP(i,j) == 0) then
             sumEnergy=sumEnergy+f2(s,i,j,k,l)*wE(S,k)*wMu(L)*eKeV(s,k)
          endif
       end do
    end do; end do; end do; end do

    dstout = factor2 * rfactor * sumEnergy

  end subroutine get_ramdst

!=============================================================================
  function Gcoul(x)
    
    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI 

    implicit none

    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_) :: G1, Gcoul
  !-----------------------------------------------------------------------
    ! G1 is erf(x) - x erf'(x)
    ! The derivative of erf(x) is (2/sqrt(pi))*exp(-x^2)
    G1=erf(X)-2.*X/sqrt(PI)*exp(-X*X)
    Gcoul=G1/2./X/X
    return
  end function Gcoul

!=============================================================================
  function funt(x)
    ! function f(y) from Ejiri, JGR, 978

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI

    implicit none
    
    real(kind=Real8_), intent(in)  :: x
    real(kind=Real8_) :: y, alpha, beta, a1, a2, a3, a4
    real(kind=Real8_) :: funt
    !-----------------------------------------------------------------------

    Y=sqrt(1-X*X)
    ALPHA=1.+LOG(2.+SQRT(3.))/2./SQRT(3.)
    BETA=ALPHA/2.-PI*SQRT(2.)/12.
    A1=0.055
    A2=-0.037
    A3=-0.074
    A4=0.056
    FUNT=ALPHA-BETA*(Y+SQRT(Y))+A1*Y**(1./3.)+ &
         A2*Y**(2./3.)+A3*Y+A4*Y**(4./3.)

    return
  end function funt

!=============================================================================
  function funi(x)
    ! Function I(y) from Ejiri, JGR, 1978

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI

    implicit none

    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_) :: y, alpha, beta, A1, A2, A3, A4
    real(kind=Real8_) :: funi, ylog
    !-----------------------------------------------------------------------
    ylog = 0.0
    Y=SQRT(1-X*X)
    if (y>0) ylog=log(y)
    ALPHA=1.+LOG(2.+SQRT(3.))/2./SQRT(3.)
    BETA=ALPHA/2.-PI*SQRT(2.)/12.
    A1=0.055
    A2=-0.037
    A3=-0.074
    A4=0.056
    FUNI=2.*ALPHA*(1.-Y)+2.*BETA*Y*ylog+4.*BETA*(Y-SQRT(Y))+ &
         3.*A1*(Y**(1./3.)-Y)+6.*A2*(Y**(2./3.)-Y)+6.*A4*(Y-Y**(4./3.)) &
         -2.*A3*Y*ylog
    
    return
  end function funi

!=============================================================================
   function atan2d(y, x)
    
    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI


    implicit none
    real(kind=Real8_), intent(in) :: x, y
    real(kind=Real8_) :: atan2d

    !-----------------------------------------------------------------------
    atan2d=180.0/pi*atan2(y,x)
    return
  end function atan2d
 
!=============================================================================
  function acosd(x)

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI
    

    implicit none
    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_) :: acosd
    !-----------------------------------------------------------------------
    acosd=180.0/pi*acos(x)
    return
  end function acosd
  
!=============================================================================
  function asind(x)

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI
    

    implicit none
    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_):: asind
    !-----------------------------------------------------------------------
    asind=180.0/pi*asin(x)
    return
  end function asind

!=============================================================================
  function cosd(x)

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI


    implicit none
    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_) :: cosd
    !-----------------------------------------------------------------------
    cosd=cos(pi/180.0 * x)
    return
  end function cosd
    
!=============================================================================
  function sind(x)

    use ModRamMain, ONLY: Real8_
    use ModRamConst, ONLY: PI


    implicit none
    real(kind=Real8_), intent(in) :: x
    real(kind=Real8_) :: sind
    !-----------------------------------------------------------------------
    sind=sin(pi/180.0 * x)
    return
  end function sind
    
!=============================================================================
  subroutine get_dipole_trace(xIn, nPoints, xOut, yOut, zOut, bxOut, byOut, bzOut)
    ! Create nPoints cartesian points along dipole from xIn to Earth's surface.
    use ModRamMain,  ONLY: Real8_
    use ModRamConst, ONLY: b0dip

    implicit none
    
    ! Argument declarations
    real(kind=Real8_), intent(in) :: xIn(3)
    integer,           intent(in) :: nPoints
    real(kind=Real8_), intent(out):: xOut(nPoints), yOut(nPoints), zOut(npoints)
    real(kind=Real8_), intent(out):: bxOut(nPoints), byOut(nPoints), bzOut(nPoints)

    ! Internal declarations
    integer       :: i
    real(kind=Real8_) :: lat, lon, r, rXyStart, dS
    real(kind=Real8_), dimension(nPoints) :: distance, rads, lats, rXy, rNew, bfact
    !-----------------------------------------------------------------------
    ! Obtain starting latitude, longitude, and radius:
    r   = sqrt(  xIn(1)**2+xIn(2)**2+xIn(3)**2 )
    lat = asin(  xIn(3)/r )
    lon = atan2( xIn(2), xIn(1) )
    
    ! Create spacing between points.
    dS = 1.0/real(nPoints)
    do i=1, nPoints
       distance(i) = real(i)*dS
    end do
    
    ! Create arrays of new locations.
    rads = r - (r-1.)*distance
    lats = acos( sqrt(rads/r*cos(lat)**2.0 ))
    
    ! Convert to new cartesian points.
    rXyStart = r*cos(lat)
    rXy      = rads*cos(lats)
    xOut = rXy*xIn(1)/rXyStart
    yOut = rXy*xIn(2)/rXyStart
    zOut = rads*sin(lats)

    ! Get dipole magnetic field values
    rNew  = sqrt(xOut**2+yOut**2+zOut**2)
    bfact = 3.0*b0dip/10**9/rNew**5
    bxOut = -bfact*xOut*zOut
    byOut = -bfact*yOut*zOut
    bzOut = -bfact*zOut*zOut + b0dip/10**9/rNew**3

  end subroutine get_dipole_trace

  !=============================================================================
  subroutine ram_sum_pressure

    use ModRamVariables, ONLY: PAllSum, PParT, PPerT, PParSum
    use ModRamGrids,     ONLY: NR, NT, nS
    use nrtype, ONLY: DP

    implicit none
    
    real(DP), parameter :: onethird=1.0/3.0, twothird=2.0/3.0
    integer :: i, j, iS
    !------------------------------------------------------------------------
   
    do i=1, nR
       do j=1, nT
          PAllSum(i,j) = 0.0
          PParSum(i,j) = 0.0
          do iS = 1, nS
             PAllSum(i,j) = PAllSum(i,j) &
                          + twothird*PPerT(iS,i,j) + onethird*PParT(iS,i,j)
             PParSum(i,j) = PParSum(i,j) + PParT(iS,i,j)
          enddo
       end do
    end do
    

  end subroutine ram_sum_pressure

!=============================================================================
  function NewtFitLarge(n, x,y,z, u, xi,yi,zi, I4)
    ! Newton interpolation with 3D vectors.
    ! X, Y, Z are length-n vectors containing the original grid.
    ! U is an n-length vector of the values to be interpolated.
    ! xi, yi, and zi are the coordinates at which to interpolate.
    ! I4 is a 4-element vector of indices specifying the 4 nearest points
    ! to our location, X,Y, and Z.
    ! The return value, NewtFit, is the value of U interpolated to xi,yi,zi.
    use ModRamMain, ONLY: Real8_

    implicit none
    ! Output value:
    real(kind=Real8_) :: NewtfitLarge
    ! Input Values:
    integer, intent(in) :: n, I4(4)
    real(kind=Real8_), intent(in) :: x(n), y(n), z(n), u(n)
    ! Local variables:
    real(kind=Real8_) :: xi, yi, zi
    real(kind=Real8_) :: f0, f1, f2, f3, D10, D12, D20, D30, D31, D32, &
         A1, A2, A3, dv0, dv1, dv2, dv3
    !-----------------------------------------------------------------------
    f0 = u(I4(1))
    f1 = u(I4(2))
    f2 = u(I4(3))
    f3 = u(I4(4))
    D10=SQRT((x(I4(2))-x(I4(1)))**2+(y(I4(2))-y(I4(1)))**2+ &
         (z(I4(2))-z(I4(1)))**2);
    D12=SQRT((x(I4(2))-x(I4(3)))**2+(y(I4(2))-y(I4(3)))**2+ &
         (z(I4(2))-z(I4(3)))**2);
    D20=SQRT((x(I4(3))-x(I4(1)))**2+(y(I4(3))-y(I4(1)))**2+ &
         (z(I4(3))-z(I4(1)))**2);
    D30=SQRT((x(I4(4))-x(I4(1)))**2+(y(I4(4))-y(I4(1)))**2+ &
         (z(I4(4))-z(I4(1)))**2);
    D31=SQRT((x(I4(4))-x(I4(2)))**2+(y(I4(4))-y(I4(2)))**2+ &
         (z(I4(4))-z(I4(2)))**2);
    D32=SQRT((x(I4(4))-x(I4(3)))**2+(y(I4(4))-y(I4(3)))**2+&
         (z(I4(4))-z(I4(3)))**2);
    
    A1=(f1-f0) / D10;
    A2=((f2-f0) / D20 + (f0-f1)/D10)/D12;
    A3=(D12*((f3-f0)/D30+(f0-f1)/D10) + &
         D31*((f1-f0)/D10+(f0-f2)/D20))/D12/D31/D32;
    dv0=SQRT((xi-x(I4(1)))**2+(yi-y(I4(1)))**2+(zi-z(I4(1)))**2)
    dv1=SQRT((xi-x(I4(2)))**2+(yi-y(I4(2)))**2+(zi-z(I4(2)))**2)
    dv2=SQRT((xi-x(I4(3)))**2+(yi-y(I4(3)))**2+(zi-z(I4(3)))**2)
    dv3=SQRT((xi-x(I4(4)))**2+(yi-y(I4(4)))**2+(zi-z(I4(4)))**2)   
    
    NewtfitLarge = f0+A1*dv0+A2*dv0*dv1+A3*dv0*dv1*dv2
    
  end function NewtFitLarge

!=============================================================================
  function NewtFit(x,y,z, u, xnew)
    ! Newton interpolation with 3D vectors.
    ! This version passes only the nearest-neighbor information and not
    ! the entire domain, as is the case with NewtFitLarge.
    ! X, Y, Z are 4 element vectors containing the nearest-neighbor coords.
    ! U is an 4-element vector of the values to be interpolated.
    ! xi, yi, and zi are the coordinates at which to interpolate.
    ! The return value, NewtFit, is the value of U interpolated to xi,yi,zi.
    
    ! This function is sensitive to the order of the variables in each
    ! vector.  For example, x(1) should correspond to the x coordinate
    ! of the first nearest neighbor, while x(4) should correspond to the
    ! fourth nearest neighbor, etc.

    use ModRamMain, ONLY: DP

    implicit none

    ! Output value:
    real(DP) :: Newtfit
    ! Input Values:
    real(DP), intent(in) :: x(4), y(4), z(4), u(4), xnew(3)
    ! Local variables:
    real(DP) :: xi, yi, zi
    real(DP) :: D10, D12, D20, D30, D31, D32, A1, A2, A3, dv0, dv1, dv2, dv3
    !-----------------------------------------------------------------------
    xi = xnew(1); yi = xnew(2); zi = xnew(3)
    D10=SQRT( (x(2)-x(1))**2 + (y(2)-y(1))**2 + (z(2)-z(1))**2 )
    D12=SQRT( (x(2)-x(3))**2 + (y(2)-y(3))**2 + (z(2)-z(3))**2 )
    D20=SQRT( (x(3)-x(1))**2 + (y(3)-y(1))**2 + (z(3)-z(1))**2 )
    D30=SQRT( (x(4)-x(1))**2 + (y(4)-y(1))**2 + (z(4)-z(1))**2 )
    D31=SQRT( (x(4)-x(2))**2 + (y(4)-y(2))**2 + (z(4)-z(2))**2 )
    D32=SQRT( (x(4)-x(3))**2 + (y(4)-y(3))**2 + (z(4)-z(3))**2 )

    A1=(u(2)-u(1)) / D10

    A2=( (u(3)-u(1))/D20 + (u(1)-u(2))/D10 ) / D12

    A3=( D12 * ((u(4)-u(1))/D30+(u(1)-u(2))/D10) + &
         D31 * ((u(2)-u(1))/D10+(u(1)-u(3))/D20) ) /D12/D31/D32

    dv0=SQRT( (xi-x(1))**2 + (yi-y(1))**2 + (zi-z(1))**2 )
    dv1=SQRT( (xi-x(2))**2 + (yi-y(2))**2 + (zi-z(2))**2 )
    dv2=SQRT( (xi-x(3))**2 + (yi-y(3))**2 + (zi-z(3))**2 )
    dv3=SQRT( (xi-x(4))**2 + (yi-y(4))**2 + (zi-z(4))**2 )   
    
    NewtFit = u(1) + A1*dv0 + A2*dv0*dv1 + A3*dv0*dv1*dv2

    ! If the result is not within the bounds of surrounding points,
    ! just use nearest neighbor.  
    if( (NewtFit > maxval(u)) .or. (NewtFit < minval(u)) ) &
         NewtFit = u(1)
 
  end function NewtFit
  !=============================================================================
  function WeightFit(n, x,y,z, u, xnew)
    ! Weighted linear interpolation (weights are distance inverses) for n 
    ! nearest neighbor points.  
    ! Inputs: 
    ! x, y, z are the locations of the n nearest neighbors.  Must be in order
    !         from nearest to farthest.  
    ! xi, yi, and zi are the coordinates of the point to which to interpolate.
    ! u are the values to be interpolated from points x,y,z to xi, yi, zi.
    ! Again, it is key that each array is placed in nearest-to-farthest order.
    ! Return value is u(n) interpolated to point xi, yi, zi.
    
    use ModRamMain, ONLY: Real8_
    

    implicit none
    
    ! Arguments and return value:
    integer, intent(in) :: n
    real(kind=Real8_), intent(in) :: x(n), y(n), z(n), u(n), xnew(3)
    real(kind=Real8_) :: WeightFit, SumWeights

    ! Internal variables
    real(kind=Real8_) ::  dv(n), xi, yi, zi
    integer :: i
    !-----------------------------------------------------------------------
    SumWeights = 0.0
    xi = xnew(1)
    yi = xnew(2)
    zi = xnew(3)
    WeightFit = 0.0
    do i = 1, n
       ! print*, 'WeightFit: value at grid point ', i, ' :', u(i)
       ! dv0(1) should be smallest
       dv(i) = (xi-x(i))**2 + (yi-y(i))**2 + (zi-z(i))**2 
    end do
    
    if (minval(dv) > 1.E-6) then ! Reference point far enough from neighbors
       do i=1, n
          SumWeights = SumWeights +  1.0 / dv(i)
          WeightFit  = WeightFit  + u(i) / dv(i)
       end do
       WeightFit = WeightFit / SumWeights
    else
       WeightFit = u(1) ! Nearest neighbor
    end if
    
    ! If the result is not within the bounds of surrounding points,
    ! just use nearest neighbor.  
    if( (WeightFit > maxval(u)) .or. (WeightFit < minval(u)) ) &
         WeightFit = u(1)

    !write(*,"(a, 3(1x, f11.7))") 'xNew = ', xi, yi, zi
    !write(*,"(a, 4(1x, f11.7))") 'xOld = ', x
    !write(*,"(a, 4(1x, f11.7))") 'yOld = ', y
    !write(*,"(a, 4(1x, f11.7))") 'zOld = ', z
    !write(*,"(a, 4(1x, f11.7))") 'dv   = ', dv
    !write(*,"(a, 4(1x, f11.7))") 'Bin  = ', u
    !write(*,*) 'Result = ', WeightFit

  END FUNCTION WeightFit
  
!=============================================================================

  subroutine test_Gcoul(verbose)
    ! Tests:
    ! 1. G function for coulomb collisions with zero input
    ! 2. As above for valid input
    use, intrinsic :: IEEE_arithmetic
    use ModRamMain,      ONLY: Real8_, test_neq_abs, nTestPassed, nTestRun

    implicit none
    logical, intent(in) :: verbose
    integer :: ii
    real(Real8_) :: test_err, test_val
    real(Real8_), dimension(3) :: expect = (/0.16221717321315, 0.213796646309081, &
                                             0.17504660173915945/)
    real(Real8_), dimension(3) :: answer
    logical :: failure, all_pass

    ! First test - zero input returns NaN
    test_val = Gcoul(0.0d0)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_gcoul (1): Gcoul(0) should equal NaN. Actual: ", test_val
    if (IEEE_IS_NAN(test_val)) then
      nTestPassed = nTestPassed + 1
    end if

    ! Second test - finite values return expected answers
    all_pass = .TRUE.
    do ii=1,3
      ! Gcoul of [0.5, 1, 1.5]
      answer(ii) = Gcoul(REAL(ii)/2.d0)
      call test_neq_abs(answer(ii), expect(ii), 1.0d-8, failure)
      if (failure) all_pass = .FALSE.
    end do
    if (verbose) then
      write(*,"(A23,3F12.8))") "test_gcoul (2): got = ", answer
      write(*,"(A23,3F12.8))") "           expected = ", expect
    end if
    nTestRun = nTestRun + 1
    if (all_pass) then
      nTestPassed = nTestPassed + 1
    end if

  end subroutine test_Gcoul

  subroutine test_NewtFit(verbose)
    ! Tests:
    ! 1. Newton Fit returns expected value on uniform surface
    ! 2. Returns nearest neighbor when point outside neighbors
    use ModRamMain,      ONLY: Real8_, test_neq_abs, nTestPassed, nTestRun

    implicit none
    logical, intent(in) :: verbose
    real(Real8_) :: test_err, test_val
    real(Real8_), dimension(4) :: x = (/0, 1, 2, 3/)
    real(Real8_), dimension(4) :: y = (/0, 0, 0, 0/)
    real(Real8_), dimension(4) :: z = (/0, 0, 0, 0/)
    real(Real8_), dimension(4) :: val = (/1, 1, 1, 1/)
    real(Real8_), dimension(4) :: val2= (/2, 3, 1, 0/)
    real(Real8_), dimension(3) :: rcon = (/0.5d0, 0.5d0, 0.0d0/)
    real(Real8_), dimension(3) :: rint = (/1.d0, 0.0d0, 0.0d0/)
    real(Real8_), dimension(3) :: rout = (/-1.0d0, -1.0d0, 0.0d0/)
    logical :: failure

    test_val = NewtFit(x, y, z, val, rcon)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_NewtFit (1): Evaluate constant. Expected: 1.0; Actual: ", test_val
    call test_neq_abs(test_val, 1.d0, 1.0d-8, failure)
    if (.not.failure) then
      nTestPassed = nTestPassed + 1
    end if

    test_val = NewtFit(x, y, z, val2, rint)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_NewtFit (2): Evaluate 1D interpolation at known point. Expected: 3.0; Actual: ", test_val
    call test_neq_abs(test_val, 3.d0, 1.0d-8, failure)
    if (.not.failure) then
      nTestPassed = nTestPassed + 1
    end if

    ! if outside range, return nearest - NOTE points are in distance order
    test_val = NewtFit(x, y, z, val2, rout)
    nTestRun = nTestRun + 1
    if (verbose) write(*,*) "test_NewtFit (3): Evaluate outside NNs. Expected: 2.0; Got: ", test_val
    call test_neq_abs(test_val, 2.d0, 1.0d-8, failure)
    if (.not.failure) then
      nTestPassed = nTestPassed + 1
    end if

  end subroutine test_NewtFit

end module ModRamFunctions
