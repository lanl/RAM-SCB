!------------------------------------------------------------------------------
!         Weimer 2005 ionospheric potential model provided by D. Weimer       !
!                Fortran version converted by Ben Foster at NCAR 2008         !
!------------------------------------------------------------------------------
module w05
  use nrtype
  implicit none

!
! Data read from W05scEpot.dat or W05scBpot.dat:
!
  integer,parameter :: csize=28, d1_pot=15, d2_pot=18
  integer :: ab(csize), ls(csize), ms(csize)
  real(dp) :: alschfits(d2_pot,csize), schfits(d1_pot,csize), ex_pot(2)
  integer :: maxl_pot,maxm_pot
!
! Data read from SCHAtable.dat:
!
  integer,parameter :: d1_scha=19, d2_scha=7, d3_scha=68
  real(dp) :: allnkm(d1_scha,d2_scha,d3_scha)
  integer :: maxk_scha,maxm_scha
  real(dp) :: th0s(d3_scha)
!
! Data read from W05scBndy.dat:
!
  integer,parameter :: na=6, nb=7
  real(dp) :: bndya(na),bndyb(nb),ex_bndy(2)

  real(dp) :: rad2deg,deg2rad           ! set by setmodel
  real(dp) :: bndyfitr                  ! calculated by setboundary
  real(dp) :: esphc(csize),bsphc(csize) ! calculated by setmodel
  real(dp) :: tmat(3,3),ttmat(3,3)      ! from setboundary
  integer,parameter :: mxtablesize=200
  real(dp) :: plmtable(mxtablesize,csize), &
          colattable(mxtablesize)
  real(dp) :: nlms(csize)

  contains
!-----------------------------------------------------------------------
  subroutine setmodel05(by,bz,tilt,swvel,swden)

    implicit none

    real(dp),intent(in) :: by,bz,tilt,swvel,swden
    integer :: i,j
    real(dp) :: bt,angle,stilt,stilt2,sw,swp,swe,c0,&
      rang,cosa,sina,cos2a,sin2a
    real(dp) :: cfits(d1_pot,csize),a(d1_pot)

    call read_potential_coeff05('W05_coeff.dat')
    
    rad2deg = 180./pi_d
    deg2rad = pi_d/180.

    bt = sqrt(by**2+bz**2)
    angle = atan2(by,bz)*rad2deg
    call setboundary05(angle,bt,tilt,swvel,swden)

    stilt = sin(tilt*deg2rad)
    stilt2 = stilt**2
    sw = bt*swvel/1000.
    swe = (1.-exp(-sw*ex_pot(2)))*sw**ex_pot(1)
    c0 = 1.
    swp = swvel**2 * swden*1.6726e-6
    rang = angle*deg2rad
    cosa = cos(rang)
    sina = sin(rang)
    cos2a = cos(2.*rang)
    sin2a = sin(2.*rang)
    if (bt < 1.) then ! remove angle dependency for IMF under 1 nT
      cosa = -1.+bt*(cosa+1.)
      cos2a = 1.+bt*(cos2a-1.)
      sina = bt*sina
      sin2a = bt*sin2a
    endif
    cfits = schfits ! schfits(d1_pot,csize) is in module read_data
    a = (/c0      , swe       , stilt      , stilt2     , swp, &
          swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa,         &
          swe*sina, stilt*sina, stilt2*sina, swp*sina,         &
          swe*cos2a,swe*sin2a/)

    esphc(:) = 0.
    do j=1,csize
       do i=1,int(d1_pot)
          esphc(j) = esphc(j)+cfits(i,j)*a(i)
       enddo
    enddo

  end subroutine setmodel05
!-----------------------------------------------------------------------
  subroutine setboundary05(angle,bt,tilt,swvel,swden)

    implicit none

    real(dp),intent(in) :: angle,bt,tilt,swvel,swden
    integer :: i
    real(dp) :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na)

! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
    xc = 4.2
    theta = xc*(deg2rad)
    ct = cos(theta)
    st = sin(theta)
!
    tmat(1,:) = (/ ct, 0._dp, st/) 
    tmat(2,:) = (/ 0., 1., 0./) 
    tmat(3,:) = (/-st, 0._dp, ct/)
!
    ttmat(1,:) = (/ct, 0._dp,-st/)
    ttmat(2,:) = (/ 0.,1., 0./)
    ttmat(3,:) = (/st, 0._dp, ct/)
!
    swp = swden*swvel**2*1.6726e-6 ! pressure
    tilt2 = tilt**2
    cosa = cos(angle*deg2rad)
    btx = 1.-exp(-bt*ex_bndy(1))
    if (bt > 1.) then
      btx = btx*bt**ex_bndy(2)
    else
      cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
    endif
    x = (/1._dp, cosa, btx, btx*cosa, swvel, swp/)
    c = bndya
    bndyfitr = 0.
    do i=1,na
      bndyfitr = bndyfitr+x(i)*c(i)
    enddo

  end subroutine setboundary05
!-----------------------------------------------------------------------
  subroutine epotval05(lat,mlt,fill,epot)

    implicit none

    real(dp),intent(in) :: lat,mlt,fill
    real(dp),intent(out) :: epot
    integer :: inside,j,m,mm,skip
    real(dp) :: z,phir,plm,colat,nlm
    real(dp) :: phim(2),cospm(2),sinpm(2)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
    call checkinputs05(lat,mlt,inside,phir,colat)
    if (inside == 0) then
      epot = fill
      return
    endif
!
! IDL code: 
! phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
! Here, maxm == maxm_pot == 2 (from read_data module), and phir is a scalar. The 
!   above IDL statement then becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where
!   phim will be dimensioned [1,2]
!
    phim(1) = phir
    phim(2) = phir*2.
!   write(6,"('epotval: phir=',1pe12.4,' phim=',2(1pe12.4))") phir,phim
    cospm(:) = cos(phim(:))
    sinpm(:) = sin(phim(:))
!
    z = 0.
    do j=1,csize
      if (skip == 1) then
        skip = 0
        cycle
      endif
      m = ms(j)
      if (ab(j)==1) then
        plm = scplm05(j,colat,nlm) ! scplm function is in this module

        skip = 0
        if (m == 0) then
          z = z+plm*esphc(j)
        else
          z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
          skip = 1
        endif

      endif ! ab(j)
    enddo
    epot = z 
  end subroutine epotval05
!-----------------------------------------------------------------------
  real(dp) function scplm05(index,colat,nlm)

    implicit none

    integer,intent(in) :: index
    real(dp),intent(in) :: colat
    real(dp),intent(out) :: nlm

    integer :: istat,i,j,l,m,skip
    real(dp) :: th0,out(1),colata(1),plm1
    real(dp) :: cth(mxtablesize)
    real(dp),save :: prevth0=1.e36
    integer,save :: tablesize

    scplm05 = 0.
    th0 = bndyfitr
    if (prevth0 /= th0) then
      tablesize = 3*nint(th0)
      if (tablesize > mxtablesize) then 
        write(6,"('>>> tablesize > mxtablesize: tablesize=',i5,&
          &' mxtablesize=',i5,' tn0=',e12.4)") tablesize,mxtablesize,th0
        stop 'tablesize'
      endif

      do i=1,tablesize
        colattable(i) = float(i-1)*(th0/float(tablesize-1))
        cth(i) = cos(colattable(i)*deg2rad)
      enddo

      prevth0 = th0
      nlms = 0. ! whole array init 
      do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        l = ls(j)
        m = ms(j)

        nlms(j) = nkmlookup05(l,m,th0) ! nkmlookup in this module

        call pm_n05(m,nlms(j),cth,plmtable(1:tablesize,j),tablesize)

        skip = 0
        if (m /= 0 .and. ab(j) > 0) then
          plmtable(1,j+1) = plmtable(1,j)
          nlms(j+1) = nlms(j)
          skip = 1
        endif

      enddo ! j=1,csize

    endif ! prevth0
    nlm = nlms(index)
    colata(1) = colat

    call interpol_quad05(plmtable(1:tablesize,index),colattable(1:tablesize),&
      colata,out)
    scplm05 = out(1)

  end function scplm05
!-----------------------------------------------------------------------
  subroutine pm_n05(m,r,cth,plmtable,tablesize)
    implicit none

    integer,intent(in) :: m,tablesize
    real(dp),intent(in) :: r
    real(dp),intent(in) :: cth(tablesize)
    real(dp),intent(out) :: plmtable(tablesize)

    integer :: i,k,ii
    real(dp) :: rm,rk,div,ans,xn
    real(dp),dimension(tablesize) :: a,x,tmp,table
!
    if (m == 0) then 
      a = 1. ! whole array op
    else
      do i=1,tablesize
        a(i) = sqrt(1.-cth(i)**2)**m
      enddo
    endif
    xn = r*(r+1.)
    x(:) = (1.-cth(:))/2.

    table = a ! whole array init

    k = 1
    pmn_loop: do         ! repeat-until loop in idl code
      do i=1,tablesize
        rm = float(m)
        rk = float(k)
        a(i) = a(i)*(x(i)*((rk+rm-1.)*(rk+rm)-xn)/(rk*(rk+rm)))
        table(i) = table(i)+a(i) ! "result" in idl code
      enddo

      k = k+1
      do i=1,tablesize
        div = abs(table(i))
        if (div <= 1.e-6) div = 1.e-6
        tmp(i) = abs(a(i)) / div
      enddo

      if (maxval(tmp) < 1.e-6) exit pmn_loop
    enddo pmn_loop
    ans = km_n05(m,r)

    plmtable(:) = table(:)*ans

  end subroutine pm_n05
!-----------------------------------------------------------------------
  real(dp) function km_n05(m,rn)
    implicit none

    integer,intent(in) :: m
    real(dp),intent(in) :: rn
    integer :: i,n
    real(dp) :: rm

    if (m == 0) then 
      km_n05 = 1.
      return
    endif
    
    rm = float(m)
    km_n05 = sqrt(2.*exp(lngamma05(rn+rm+1.)-lngamma05(rn-rm+1.))) / (2.**m*factorial05(m))

  end function km_n05
!-----------------------------------------------------------------------
  real(dp) function nkmlookup05(k,m,th0)

    implicit none
!
! Args:
    integer,intent(in) :: k,m
    real(dp),intent(in) :: th0
!
! Local:
    integer :: kk,mm
    real(dp) :: th0a(1),out(1)

    if (th0 == 90.) then
      nkmlookup05 = float(k)
      return
    endif
    th0a(1) = th0
    kk = k+1
    mm = m+1
    if (kk > maxk_scha) then
      write(6,"('>>> nkmlookup: kk > maxk: kk=',i4,' maxk=',i4)") kk,maxk_scha
      call interpol_quad05(allnkm(maxk_scha,mm,:),th0s,th0a,out)
    endif
    if (mm > maxm_scha) then
      write(6,"('>>> nkmlookup: mm > maxm: kk=',i4,' maxm=',i4)") kk,maxm_scha
      call interpol_quad05(allnkm(kk,maxm_scha,:),th0s,th0a,out)
    endif
    if (th0 < th0s(1)) then
      write(6,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,' th0s(1)=',e12.4)") &
        th0,th0s(1)
    endif

!   write(6,"('nkmlookup call interpol: kk=',i3,' mm=',i3,' th0=',e12.4,&
!     &' allnkm=',/,(6(1pe12.4)))") kk,mm,th0a,allnkm(kk,mm,:)

    call interpol_quad05(allnkm(kk,mm,:),th0s,th0a,out)

    nkmlookup05 = out(1)

  end function nkmlookup05
!-----------------------------------------------------------------------
  subroutine checkinputs05(lat,mlt,inside,phir,colat)
    implicit none
!
! Args:
    real(dp),intent(in) :: lat,mlt
    integer,intent(out) :: inside
    real(dp),intent(out) :: phir,colat
!
! Local:
    real(dp) :: lon,tlat,tlon,radii

    lon = mlt*15.
    call dorotation05(lat,lon,tlat,tlon)
    radii = 90.-tlat
    inside = 0
    if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
    phir = tlon*deg2rad
    colat = radii

  end subroutine checkinputs05
!-----------------------------------------------------------------------
  subroutine dorotation05(latin,lonin,latout,lonout)
    implicit none
!
! Args:
    real(dp),intent(in) :: latin,lonin
    real(dp),intent(out) :: latout,lonout
!
! Local:
    real(dp) :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
    integer :: i

    latr = latin*deg2rad
    lonr = lonin*deg2rad
    stc = sin(latr)
    ctc = cos(latr)
    sf = sin(lonr)
    cf = cos(lonr)
    a = ctc*cf
    b = ctc*sf
!
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
!   Pos is then (1,3)
!
    do i=1,3
      pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
    enddo
  
    latout = asin(pos(3))*rad2deg
    lonout = atan2(pos(2),pos(1))*rad2deg

!   write(6,"('dorotation: latin,lonin=',2f9.4,' latout,lonout=',2f9.4)") &
!     latin,lonin,latout,lonout
  
  end subroutine dorotation05
!-----------------------------------------------------------------------
  subroutine interpol_quad05(v,x,u,p)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
    implicit none
!
! Args:
    real(dp),intent(in) :: v(:),x(:),u(:)
    real(dp),intent(out) :: p(:)
!
! Local:
    integer :: nv,nx,nu,i,ix
    real(dp) :: x0,x1,x2
!
    nv = size(v)
    nx = size(x)
    nu = size(u)
    if (nx /= nv) then
      write(6,"('>>> interpol_quad: nx /= nv: nx=',i4,' nv=',i4)") nx,nv
      p(:) = 0.
      return
    endif
    do i=1,nu
      ix = value_locate05(x,u(i))
      if (ix <= 1.or.ix >= nx) then ! bug fix by btf 12/23/09
          p(i) = 0.
          cycle                       ! bug fix by btf 12/23/09
      endif
      x1 = x(ix)
      x0 = x(ix-1)
      x2 = x(ix+1)
      p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) + &
             v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) + &
             v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
      
    enddo
!   write(6,"('interpol_quad: nu=',i4,' p=',/,(1pe12.4)") nu,p

  end subroutine interpol_quad05
!-----------------------------------------------------------------------
  integer function value_locate05(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
!
    implicit none
!
! Args:
    real(dp),intent(in) :: vec(:),val
!
! Local:
    integer :: n,i
!
    value_locate05 = 0
    n = size(vec)
    if (val < vec(1)) return
    if (val > vec(n)) then
      value_locate05 = n
      return
    endif
    do i=1,n-1
      if (val >= vec(i) .and. val <= vec(i+1)) then
        value_locate05 = i
        return
      endif
    enddo

  end function value_locate05
!-----------------------------------------------------------------------
real(dp) function lngamma05(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
  implicit none
  real(dp),intent(in) :: xx
  real(dp) :: x,y,tmp,ser
  real(dp) :: cof(6) = (/76.18009172947146, -86.50532032941677, 24.01409824083091,&
                     -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5/)
  integer :: j
!
  y = xx
  x = xx
  tmp = x+5.5
  tmp = tmp-(x+0.5)*log(tmp)
  ser = 1.000000000190015
  do j=1,5
    y = y+1
    ser = ser+cof(j)/y
  enddo
  lngamma05 = -tmp+log(2.5066282746310005*ser/x)
end function lngamma05
!-----------------------------------------------------------------------
real(dp) function factorial05(n)
  implicit none
  integer,intent(in) :: n
  integer :: m
  if (n <= 0) then
    write(6,"('>>> factorial: n must be positive: n=',i4)") n
    factorial05 = 0.
    return
  endif
  if (n == 1) then
    factorial05 = 1.
    return
  endif
  factorial05 = float(n)
  do m = n-1,1,-1
    factorial05 = factorial05 * float(m)
  enddo
end function factorial05
!-----------------------------------------------------------------------
subroutine read_potential_coeff05(infile)

  ! read three data files combined from W04scEpot.dat, SCHAtable.data, 
  ! and W05scBndy.data
  implicit none

! Args:
  character(len=*),intent(in) :: infile 

! Local:
  character(len=16) :: fname
  integer :: i,j,lu=20
  integer :: csize_rd,d1_rd,d2_rd
  integer :: rd_na,rd_nb
! ----------------------------------
! Read ascii data file W05scEpot.dat
! ---------------------------------- 

  PRINT *,infile
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(28i3)") ab
  read(lu,"(3i3)") csize_rd,d1_rd,d2_rd
  if (csize_rd /= csize) then 
    write(6,"('>>> read_potential: file ',a,': incompatable csize: ',&
      &'csize_rd=',i4,' csize=',i4)") fname,csize_rd,csize
    stop 'csize'
  endif
  if (d1_rd /= d1_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d1: ',&
      &'d1_rd=',i4,' d1_pot=',i4)") fname,d1_rd,d1_pot
    stop 'd1'
  endif
  if (d2_rd /= d2_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d2: ',&
      &'d2_rd=',i4,' d2_pot=',i4)") fname,d2_rd,d2_pot
    stop 'd2'
  endif
  do i=1,csize
    read(lu,"(6e20.9)") alschfits(:,i)
  enddo
  read(lu,"(2f10.3)") ex_pot
  read(lu,"(28i3)") ls
  read(lu,"(2i3)") maxl_pot,maxm_pot
  read(lu,"(28i3)") ms
  do i=1,csize 
    read(lu,"(6e20.9)") schfits(:,i)
  enddo

!-----------------------------------
! Read ascii data file SCHAtable.dat
! ----------------------------------
  read(lu,"(a)") fname
  read(lu,"(2i3)") maxk_scha,maxm_scha
  do i=1,d3_scha
    do j=1,d2_scha
      read(lu,"(6e20.9)") allnkm(:,j,i)
    enddo
  enddo
  read(lu,"(8f10.4)") th0s

!-----------------------------------
! Read ascii data file W05scBndy.dat
!-----------------------------------
  read(lu,"(a)") fname
  read(lu,"(2i3)") rd_na,rd_nb
  if (rd_na /= na) then 
    write(6,"('>>> read_potential: file ',a,': incompatable na: ',&
      &'rd_na=',i4,' na=',i4)") fname,rd_na,na
    stop 'na'
  endif
  if (rd_nb /= nb) then 
    write(6,"('>>> read_potential: file ',a,': incompatable nb: ',&
      &'rd_nb=',i4,' nb=',i4)") fname,rd_nb,nb
    stop 'nb'
  endif
  read(lu,"(8e20.9)") bndya
  read(lu,"(8e20.9)") bndyb
  read(lu,"(8e20.9)") ex_bndy

  close(lu)

end subroutine read_potential_coeff05
!-----------------------------------------------------------------------
end module w05
