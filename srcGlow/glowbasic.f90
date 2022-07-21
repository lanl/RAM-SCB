program glowbasic

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.981, 6/2017

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

  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,zlbh
  use cglow,only: cglow_init
  use cglow,only: data_dir
  use cglow,only: ec_diff,ef_diff,ntimes_traj,glons,glats,mlats,mlons,mlts,nE,ec_all,ef_all,sigmap_r,sigmah_r
  use cglow, only: iyear_traj,imonth_traj,iday_traj,ihour_traj,iminute_traj,isecond_traj

  implicit none
  
  character(len=1024) :: iri90_dir,filename
  real, parameter  :: cpi = 3.141592654
  real,allocatable :: z(:)                    ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:), zvn(:)          ! neutral wind components (not in use)
  real,allocatable :: pedcond(:,:), hallcond(:,:) ! Pederson and Hall conductivities in S/m (mho)
  real,allocatable :: outf(:,:)               ! iri output (11,jmax)
  real,allocatable :: hallconductance(:),pedconductance(:)
  real,allocatable :: eDensity(:,:),eDensity_new(:,:),ionrate(:,:)
  real :: rz12,stl,fmono,emono, f107y, mlat,mlon
  real :: d(8), t(2), sw(25), oarr(30),logener(nbins),logphitop(nbins)
  real :: logec_diff(100),logef_diff(100), dlogE
  integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,itail,ierror,ie
  integer :: instance,iostatus,idoy,ndaymo,unittmp_
  logical :: jf(12), UseSpec, UseDmsp
  character(len=99) ::n_jmax
  data sw/25*1./

  ! false if only Eave and Fe is used for calculation (call maxt)
  UseSpec = .True.
  UseDmsp = .False.
!
! Initialize standard switches:
!
  iscale=1
  xuvfac=3.
  kchem =4
  jlocal=0
  itail = 0 ! itail: 1 -- add the low-energy tail (Meier et al 1989); itail: 0 -- no tail
  fmono=0.
  emono=0.
!
! Set data directories:
!
  data_dir    = 'data/'
  iri90_dir   = 'data/iri90/'
!
! Set number of altitude levels:
!
!  jmax = 102
  jmax = 112
!
!  jmax = 68

! Allocate local arrays:
!
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(outf(11,jmax))
!
! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
  call cglow_init
!
! Call EGRID to set up electron energy grid:
!
  call egrid (ener, del, nbins)

! read in the input file of precipitation, the location of the satellite trajectory (in glat, glon)
  if (UseSpec .and. .not. UseDmsp)then
     filename = 'inputs/GLOW_ram_flux_input_one_point.dat'
     call read_trajectory_flux(filename)
  elseif (UseSpec .and. UseDmsp)then
     filename = 'inputs/Glow_dmsp_input_flux.dat'
     call read_trajectory_flux(filename)
  else !!! test at one location with one Maxwellian input
     nTimes_traj = 50
     allocate(glats(nTimes_traj),glons(nTimes_traj), ec_all(ntimes_traj), ef_all(ntimes_traj))

     glats(:) = 70.
     glons(:) = 180.

     ef_all(:) = 1.  ! ergs/cm^2/s
     ec_all(1) = 0.5 !!! define a energy grid from 500 eV to 50 keV
     dlogE = (log10(50.) - log10(0.5))/nTimes_traj
     do ie=2, nTimes_traj
        ec_all(ie) = 10**(log10(0.5) + dlogE*ie)
     end do

  end if

!
! Loop to call GLOW for designated inputs until end-of-file or any character on standard input:
!
  allocate(pedcond(ntimes_traj,jmax))
  allocate(hallcond(ntimes_traj,jmax))
  allocate(pedconductance(ntimes_traj),hallconductance(ntimes_traj))
  allocate(ionrate(ntimes_traj,jmax), eDensity(ntimes_traj,jmax),eDensity_new(ntimes_traj,jmax))

  Pedconductance  = 0.0
  Hallconductance = 0.0
  pedcond  = 0.0
  hallcond = 0.0

  ionrate = 0.0
  eDensity= 0.0
  eDensity_new = 0.0

  do instance=1,nTimes_traj
!
! Get input values:
!

     glat  = glats(instance)
     glong = glons(instance)
!     call geomag(1,glong,glat,mlon,mlat)
!     if (glong .lt. 0)glong = glong + 360.

     if(UseSpec)then
        ut = ihour_traj(instance) * 3600.+iminute_traj(instance)*60+isecond_traj(instance)
        call moda(0, imonth_traj(instance), iday_traj(instance), idoy) ! moda in iri90.f
        idate = (iyear_traj(instance) - iyear_traj(instance)/100*100)*1000+idoy
        ! get the ap data from "apf107.dat"
        call apf_only(iyear_traj(instance), imonth_traj(instance), iday_traj(instance),&
             f107, f107p, f107a, f107y, ap)                  ! apf_only in irifun_2012.f

     else !!! test at one location with one Maxwellian input
        ut = 9*3600.
        idate = 2005256
        ap = 5
        f107  = 50
        f107a = 50
        f107p = 50
     end if

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
    if (ef_all(instance) > 0.0001 .and. ec_all(instance) > 1)then
       if(.not. UseSpec)then

          ec = ec_all(instance)*1.0e3 ! -> eV
          ef = ef_all(instance) ! ergs/cm^2/s

          ! from ef and ec get a maxiwellian function (ef, ec) [ergs/cm^2/s, eV] --> phitop: /cm^2/s/eV
!!! test: temporarily no limit on the ef and ec
          !if (ef>.001 .and. ec>1.) 
          write(*,*)'ec, ef:',ec, ef
          
          call maxt(ef,ec,ener,del,nbins,itail,fmono,emono,phitop)

          !----------------------------------------------------------------
          !!! for those out of 500 eV to 50 keV, force the flux to be zero (to be consistent with Robinson)
          !----------------------------------------------------------------
          do i=1, nbins
             if(ener(i) .lt. 500 .or. ener(i) .gt. 50000.)then
                phitop(i) = 0.0
             end if
          end do
          
       else
          ! interpolate in energy space (logE, logflux) (eV, /cm^2/s/eV)
          logener = log10(ener)
          k = 0
          
          do i=1,nbins
             ! NOTE: the input differential flux is in unit of /cm^2/s/sr/keV. 
             ! NOTE: "phitop" needs /cm^2/s/eV.
             logef_diff(1:nE) = log10(ef_diff(instance,:)*cpi/1000.) ! convert /cm^2/s/sr/keV to /cm^2/s/eV
             logec_diff(1:nE) = log10(ec_diff*1.0e3) ! convert to eV
             
             if (logener(i).gt.logec_diff(2) .and. logener(i) .le. logec_diff(nE)) then
                call lintp(logec_diff(2:nE), logef_diff(2:nE), nE-1,logener(i), logphitop(i), ierror)
                phitop(i) = 10.**(logphitop(i))
             else
                phitop(i) = 0.0
                k = k + 1
             end if
             
             !\
             ! make the spectrum within the Robinson's valid energy range (for comparison to Robinson's)
             ! comment out to take into account the full spectrum
             !/
             !if(ener(i) .lt. 500 .or. ener(i) .gt. 50000.)then
             !   phitop(i) = 0.0
             !end if
             
          end do
       end if

    end if
! Fill altitude array, converting to cm:
!
    zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!

    if(maxval(phitop) > 0.0)then
       call glow

!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
       do j=1,jmax
          call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
               zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
               pedcond(instance,j), hallcond(instance,j))
          !      ! Yiqun Yu 2017/02: only include the conductivity below 200 km...
          if (j .gt. 1 .and. z(j) .le. 200. .and. z(j) .ge. 75.) then
             ! calculate height-integrated conductance
             pedconductance(instance) = pedconductance(instance)   + pedcond(instance,j)*(z(j)-z(j-1))*1.0e3 !S/m*(km*1.0e3)
             hallconductance(instance) = hallconductance(instance) + hallcond(instance,j)*(z(j)-z(j-1))*1.0e3 !S/m*(km*1.0e3)
          end if
          !      write(*,*)"z(j), PedCond(instance,j):",z(j),z(j-1), PedConductance(instance)
          
       end do
    end if

   ionrate(instance, :) = tir(1:jmax)
   eDensity(instance,:) = ze(1:jmax)
   eDensity_new(instance,:) = ecalc(1:jmax)

!!
!! Output section:
!!
!    write(6,"(1x,i7,9f8.1)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
!    write(6,"('   Z     Tn       O        N2        NO      Ne(in)    Ne(out)  Ionrate      O+       O2+      NO+       N(2D)    Pederson   Hall')")
!    write(6,"(1x,0p,f5.1,f6.0,1p,12e10.2)") (z(j),ztn(j),zo(j),zn2(j),zno(j),ze(j), &
!      ecalc(j),tir(j),zxden(3,j),zxden(6,j),zxden(7,j),zxden(10,j),pedcond(j),hallcond(j),j=1,jmax)
!    write(6,"('   Z      3371    4278    5200    5577    6300    7320   10400    3644    7774    8446    3726    LBH     1356    1493    1304')")
!    write(6,"(1x,f5.1,15f8.2)")(z(j),(zeta(ii,j),ii=1,15),j=1,jmax)


 enddo
 
 ! output the conductivity for all the times
 unittmp_ = 1
 if (UseSpec .and. .not.UseDmsp)then
    open(unittmp_, file='glow_height_cond_ram_fullSpec.dat',status='unknown')
 elseif (.not. UseSpec .and. UseDmsp)then
    open(unittmp_, file='glow_height_cond_dmsp_max.dat',status='unknown')
 elseif (UseSpec .and. UseDmsp)then
    open(unittmp_, file='glow_height_cond_dmsp_fullSpec.dat',status='unknown')
 elseif (.not. UseSpec .and. .not. UseDMSP .and. itail == 1)then
    open(unittmp_, file='glow_height_cond_max_lowtail.dat',status='unknown')
 else
    open(unittmp_, file='glow_height_cond_max.dat',status='unknown')
 end if

! do j=1, jmax
!    if (z(j) .gt. 200.)exit
! end do

 write(unittmp_,*)'nHeight: ', jmax!j-1
 do i=1, nTimes_traj
    if(UseSpec)then
       write(unittmp_, *) &
            iyear_traj(i), imonth_traj(i),iday_traj(i),ihour_traj(i),&
            iminute_traj(i),isecond_traj(i), glats(i), glons(i), mlats(i), mlons(i),mlts(i), &
            sigmap_r(i), sigmah_r(i), pedconductance(i), hallconductance(i), ec_all(i), ef_all(i)
    else
       write(unittmp_, *) pedconductance(i), hallconductance(i), ec_all(i), ef_all(i)
    end if

    write(unittmp_,*)'Z   Ionrate   Ne(in)   Ne(out)   Pedersen   Hall'
    do j=1,jmax
!       if(z(j) .le. 200)&
            write(unittmp_, '(f6.2, 8e12.3)')z(j), ionrate(i,j), eDensity(i,j), eDensity_new(i,j), &
                 zxden(3,j),zxden(6,j),zxden(7,j), pedcond(i,j), hallcond(i,j)
    end do
 end do
 close(unittmp_)

if(UseSpec)then
   deallocate(ec_diff,ef_diff,glats,glons,mlats,mlons,mlts)
else
   deallocate(ec_all, ef_all)
end if
if(UseSpec)&
     deallocate(iyear_traj,imonth_traj,iday_traj,ihour_traj,iminute_traj,isecond_traj)
 
deallocate(z,zun,zvn,outf, pedcond, hallcond,Pedconductance,Hallconductance)

end program glowbasic
