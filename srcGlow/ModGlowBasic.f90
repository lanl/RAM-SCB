module ModGlowBasic


  use cglow,only: jmax,nbins,npbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,pf,pc
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: pflux,pener,pdel,pia
  use cglow,only: cglow_init
  use cglow,only: data_dir

  implicit none
  save

  contains
!================================================================================    
  subroutine glowbasic_ram(idate_in, ut_in, glat_in, glong_in, ap_in, f107_in, &
       f107p_in, f107a_in, ef_in, ec_in, pedconductance, hallconductance, nE, &
       height, ionrate, eDensity, Pedcond, Hallcond, nHeight,&
       logec_diff, logef_diff, pf_in, pc_in, logpc_diff, logpf_diff)


! Modified from GLOW's original program: glowbasic.f90 in order to couple with RAM. 
! By Yiqun Yu 2022

! optional input variables: 
! logec_diff, logef_diff:                     if present, then use the full spectrum as input
! pf_in, pc_in, logpc_diff, logpf_diff: if present, then calculate proton impact

!====================================================================
! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file Glowlicense.txt.
! For more information see the file Glow.txt.

! Version 0.982, 2/2018
! Adapted from glowdriver by Stan Solomon, 2/2016

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
! Runs GLOW for designated inputs once, or multiple times.
! MPI and netCDF libraries not required.

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission

  implicit none

  integer, parameter :: Real8_ = selected_real_kind(12,100)
  integer, intent(in):: idate_in, nHeight
  real(real8_), intent(in) :: ut_in, glat_in, glong_in, ef_in, ec_in
  real(real8_), intent(in) :: ap_in, f107_in, f107p_in, f107a_in
  real(real8_), intent(out) :: pedconductance, hallconductance
  real(real8_), dimension(nHeight), intent(out) :: &
       height, ionrate, eDensity,Pedcond,Hallcond
  integer, intent(in) :: nE
  real(real8_), optional, intent(in) :: logec_diff(nE), logef_diff(nE), logpc_diff(nE), logpf_diff(nE)
  real(real8_), optional, intent(in) :: pf_in, pc_in ! proton flux/Ec 
  character(len=1024) :: iri90_dir
  
  real,allocatable :: z(:)                    ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:), zvn(:)          ! neutral wind components (not in use)
!  real,allocatable :: pedcond(:), hallcond(:) ! Pederson and Hall conductivities in S/m (mho)
  real,allocatable :: outf(:,:)               ! iri output (11,jmax)
  real :: rz12,stl,fmono,emono, pfmono,pemono
  real :: d(8), t(2), sw(25), oarr(30),logener, logphitop(nbins), logpflux(npbins)
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,itail, ierror
  integer :: instance,iostatus
  logical :: jf(12)
  data sw/25*1./

!  write(*,*)'present(logec_diff):',present(logec_diff)
  idate = idate_in
  ut = ut_in
  glat = glat_in
  glong= glong_in
  ap = int(ap_in)
  f107 = f107_in
  f107a = f107a_in
  f107p = f107p_in

  ef = ef_in
  ec = ec_in

  ! assume no proton precipitation for now
  if(present(pf_in))then
     pf = pf_in
     pc = pc_in
  else
     pf = 0.0
     pc = 0.0
  end if
!
! Initialize standard switches:
!
  iscale=1
  xuvfac=3.
  kchem=4
  jlocal=0
  itail=0
  fmono=0.
  emono=0.
  pfmono=0.
  pemono=0.
!
! Set data directories:
!
  data_dir    = 'IM/input_sce/glow_data/'
  iri90_dir   = 'IM/input_sce/glow_data/iri90/'
!
! Set number of altitude levels:
!
  jmax = nHeight
! jmax = 86
!
! Allocate local arrays:
!
  if (.not.allocated(z))then
     allocate(z(jmax))
     allocate(zun(jmax))
     allocate(zvn(jmax))
!     allocate(pedcond(jmax))
!     allocate(hallcond(jmax))
     allocate(outf(11,jmax))
  end if

! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)

  call cglow_init

  pedconductance = 0.
  hallconductance= 0.
  pedcond = 0.0
  hallcond= 0.0
!
! Call EGRID to set up electron energy grid (eV):
!
  call egrid (ener, del, nbins)
!
! Call PEGRID to set up proton energy grid:
!
  call pegrid (pener, pdel, npbins)
!
! Loop to call GLOW for designated inputs until end-of-file or any character on standard input:
!
!  do instance=1,10000
!
! Get input values:
!
!    write(6,"('Enter date, UT, lat, lon, F107a, F107, F107p, Ap, Ef, Ec')")
!    read(5,*,iostat=iostatus) idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
!    if (iostatus /= 0) stop
!
! Calculate local solar time:
!
    stl = ut/3600. + glong/15.
    if (stl < 0.) stl = stl + 24.
    if (stl >= 24.) stl = stl - 24.
!
! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
!

    call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                 z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Call MAXT to put auroral electron flux specified by namelist input into phitop array:
!
    
    phitop(:) = 0.
    if(ef>.0001 .and. ec>1)then ! only when the precipitaiton flux is large enough
       if (.not. present(logec_diff))then
          ! from ef and ec get a maxiwellian function (ef, ec) [ergs/cm^2/s, eV] --> phitop: /cm^2/s/eV
          call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)

          do i=1,nbins
             if(ener(i) .lt. 500 .or. ener(i) .gt. 50000.)then
                phitop(i) = 0.0
             end if
          end do

       else
          !interpolate in energy space (logE, logflux) (eV, /cm^2/s/eV)
          k = 0
          do i=1,nbins
             logener = log10(ener(i))

             if (logener .ge. logec_diff(2) .and. logener .le. logec_diff(nE))then
             !!! test: limit the energy range within the Robinson's energy range.
!             if(logener .ge. log10(500.) .and. logener .le. log10(50000.))then !  500 eV < E < 50 keV
                call lintp(logec_diff(2:), logef_diff(2:), nE-1, &         ! logef_diff: /cm^2/s/eV
                     logener, logphitop(i), ierror)
                phitop(i) = 10.**(logphitop(i))
             else
                phitop(i) = 0.0
                k = k + 1
             end if

          end do
          ! fill up the lower-E part, simply let it be same as the first available energy grid (test)
          !phitop(1:k) = phitop(k+1)
       end if

    end if

!
! Call MAXT to put auroral proton flux specified by input into pflux array:
!
    pflux(:) = 0.
    if (pf>.001 .and. pc>1.)then
       if (.not. present(logpc_diff))then
          call maxt (pf,pc,pener,pdel,npbins,itail,pfmono,pemono,pflux)
       else
          k = 0
          do i=1,npbins
             logener = log10(pener(i))
             if (logener .ge. logpc_diff(2) .and. logener .le. logpc_diff(nE))then
                call lintp(logpc_diff(2:), logpf_diff(2:), nE-1, &         ! logef_diff: /cm^2/s/eV
                     logener, logpflux(i), ierror)
                pflux(i) = 10.**(logpflux(i))
             else
                pflux(i) = 0.0
                k = k + 1
             end if

          end do
       end if
    end if

!
! Fill altitude array, converting to cm:
!
    zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!

!\
!!! if only the precipitaiton is considered, use the following if statement; 
!!! if the solar radiation is also included, no need to have the if statement
!/
   if (maxval(phitop) > 0.0)then
       call glow

! Call CONDUCT to calculate Pederson and Hall conductivities:
!
       do j=1,jmax

          call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
               zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
               pedcond(j), hallcond(j))

          ! Yiqun Yu 2017/02, only include the conductivity below 200km
          if (j .gt. 1 .and. z(j) .le. 200 .and. z(j) .gt. 75) then
             ! calculate height-integrated conductance
             pedconductance = pedconductance + pedcond(j)*(z(j)-z(j-1))*1.0e3    !S/m*(km*1.0e3)
             hallconductance = hallconductance + hallcond(j)*(z(j)-z(j-1))*1.0e3 !S/m*(km*1.0e3)
          end if
       enddo
    end if

    height   = z
    ionrate  = tir
    eDensity = ecalc

!
! Output section:
!
!    write(6,"(1x,i7,9f8.1)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
!    write(6,"('   Z     Tn       O        N2        NO      Ne(in)    Ne(out)  Ionrate      O+       O2+      NO+       N(2D)    Pederson   Hall')")
!    write(6,"(1x,0p,f5.1,f6.0,1p,12e10.2)") (z(j),ztn(j),zo(j),zn2(j),zno(j),ze(j), &
!      ecalc(j),tir(j),zxden(3,j),zxden(6,j),zxden(7,j),zxden(10,j),pedcond(j),hallcond(j),j=1,jmax)
!    write(6,"('   Z      3371    4278    5200    5577    6300    7320   10400    3644    7774    8446    3726    LBH     1356    1493    1304')")
!    write(6,"(1x,f5.1,15f8.2)")(z(j),(zeta(ii,j),ii=1,15),j=1,jmax)

!  enddo

!stop

  end subroutine glowbasic_ram
end module ModGlowBasic
