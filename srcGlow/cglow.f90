module cglow

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.981, 6/2017

! Stan Solomon and Ben Foster, 1/2015
! Stan Solomon, 1/2016, 3/2016, consolidated with cxglow
! Stan Solomon, 6/2017, zeroed out arrays

! CGLOW Defines array dimensions and use-associated variables for the GLOW model.
! Replaces the header file glow.h and common blocks /CGLOW/, /CXSECT/, and /CXPARS/
! that were used in older versions of the model (v. 0.973 and earlier).

! For variable definitions, see subroutine GLOW and subroutine EXSECT.

! Old common blocks, for reference:

!     COMMON /CGLOW/
!    >    IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM,
!    >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
!    >    ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
!    >    ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
!    >    ZTN(JMAX), ZTI(JMAX), ZTE(JMAX),
!    >    PHITOP(NBINS), EFLUX(NF), EZERO(NF),
!    >    SZA, DIP, EFRAC, IERR,
!    >    ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX),
!    >    WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
!    >    ENER(NBINS), DEL(NBINS),
!    >    PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
!    >    PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX),
!    >    QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX),
!    >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
!    >    EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX),
!    >    ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW)

!     COMMON /CXSECT/ SIGS(NMAJ,NBINS), PE(NMAJ,NBINS), PIN(NMAJ,NBINS),
!    >                SIGA(NMAJ,NBINS,NBINS), SEC(NMAJ,NBINS,NBINS),
!    >                SIGEX(NEI,NMAJ,NBINS), SIGIX(NEI,NMAJ,NBINS),
!    >                IIMAXX(NBINS)

!     COMMON /CXPARS/ WW(NEI,NMAJ), AO(NEI,NMAJ), OMEG(NEI,NMAJ),
!    >                ANU(NEI,NMAJ), BB(NEI,NMAJ), AUTO(NEI,NMAJ),
!    >                THI(NEI,NMAJ),  AK(NEI,NMAJ),   AJ(NEI,NMAJ),
!    >                TS(NEI,NMAJ),   TA(NEI,NMAJ),   TB(NEI,NMAJ),
!    >                GAMS(NEI,NMAJ), GAMB(NEI,NMAJ)


  implicit none
  save

! Array dimensions, configurable:

  integer :: jmax                ! number of vertical levels

! Array dimensions, non-configurable:

!  integer,parameter :: nbins=190 ! number of energetic electron energy bins
  integer,parameter :: nbins=215 ! number of energetic electron energy bins
  integer,parameter :: npbins=27 ! number of energetic proton precipitation energy bins
  integer,parameter :: lmax=123  ! number of wavelength intervals for solar flux
  integer,parameter :: nmaj=3    ! number of major species
  integer,parameter :: nst=6     ! number of states produced by photoionization/dissociation
  integer,parameter :: nei=10    ! number of states produced by electron impact
  integer,parameter :: nex=12    ! number of excited/ionized species
  integer,parameter :: nw=15     ! number of airglow emission wavelengths
  integer,parameter :: nc=10     ! number of component production terms for each emission

! Directory containing data files needed by glow subroutines:

  character(len=1024) :: data_dir

  integer :: idate,iscale,jlocal,kchem,ierr
  real    :: ut,glat,glong,f107,f107a,f107p,ap,ef,ec,pf,pc
  real    :: xuvfac, sza, dip, efrac

  real,allocatable,dimension(:) ::             &                   ! (jmax)
    zz, zo, zn2, zo2, zno, zns, znd, zrho, ze, &
    ztn, zti, zte, eheat, tez, ecalc, tei, tpi, tir
  real,allocatable,dimension(:)     :: phitop, ener, del           ! (nbins) 
  real,allocatable,dimension(:)     :: pflux, pener, pdel          ! (npbins) 
  real,allocatable,dimension(:)     :: wave1, wave2, sflux         ! (lmax)
  real,allocatable,dimension(:,:)   :: pespec, sespec, uflx, dflx  ! (nbins,jmax)
  real,allocatable,dimension(:,:)   :: zmaj, zcol, pia, sion       ! (nmaj,jmax)
  real,allocatable,dimension(:,:,:) :: photoi, photod              ! (nst,nmaj,jmax)
  real,allocatable,dimension(:,:)   :: phono                       ! (nst,jmax)
  real,allocatable,dimension(:,:,:) :: aglw                        ! (nei,nmaj,jmax)
  real,allocatable,dimension(:,:)   :: zxden                       ! (nex,jmax)
  real,allocatable,dimension(:,:)   :: zeta                        ! (nw,jmax)
  real,allocatable,dimension(:,:,:) :: zceta                       ! (nc,nw,jmax)
  real,allocatable,dimension(:)     :: vcb                         ! (nw)
  real,allocatable,dimension(:,:)   :: zlbh                        ! (nc,jmax)
  real,allocatable,dimension(:,:)   :: sigs,pe,pin                 ! (nmaj,nbins)
  real,allocatable,dimension(:,:,:) :: sigex,sigix                 ! (nei,nmaj,nbins)
  real,allocatable,dimension(:,:,:) :: siga,sec                    ! (nei,nbins,nbins)
  integer,allocatable,dimension(:)  :: iimaxx                      ! (nbins)
  real,allocatable,dimension(:,:) :: &                             ! (nei,nmaj)
    ww,ao,omeg,anu,bb,auto,thi,ak,aj,ts,ta,tb,gams,gamb

  ! for input flux file reading
  real, allocatable,dimension(:,:) :: ef_diff, pf_diff
  real, allocatable,dimension(:) :: ec_diff,pc_diff, glats,glons,mlats,mlons,mlts
  integer, allocatable, dimension(:) :: iyear_traj,imonth_traj,iday_traj,&
       ihour_traj,iminute_traj,isecond_traj
  integer :: ntimes_traj,nE
  real, allocatable, dimension(:) :: ec_all, ef_all,pc_all,pf_all, sigmap_r, sigmah_r
  contains
!-----------------------------------------------------------------------

  subroutine cglow_init

! Allocate variable arrays:

    if (.not.allocated(zz))then
    allocate        &
      (zz   (jmax), &
       zo   (jmax), &
       zn2  (jmax), &
       zo2  (jmax), &
       zno  (jmax), &
       zns  (jmax), &
       znd  (jmax), &
       zrho (jmax), &
       ze   (jmax), &
       ztn  (jmax), &
       zti  (jmax), &
       zte  (jmax), &
       eheat(jmax), &
       tez  (jmax), &
       tei  (jmax), &
       tpi  (jmax), &
       tir  (jmax), &
       ecalc(jmax))

    allocate            &
      (zxden(nex,jmax), &
       zeta(nw,jmax),   &
       zceta(nc,nw,jmax), &
       vcb(nw),           &
       zlbh(nc,jmax))

    allocate          &
      (phitop(nbins), &
       ener  (nbins), &
       del   (nbins))

    allocate          &
      (pflux(npbins), &
       pener(npbins), &
       pdel(npbins))

    allocate          &
      (wave1(lmax),   &
       wave2(lmax),   &
       sflux(lmax))

    allocate               &
      (pespec(nbins,jmax), &
       sespec(nbins,jmax), &
       uflx  (nbins,jmax), &
       dflx  (nbins,jmax))

    allocate &
      (zmaj(nmaj,jmax), &
       zcol(nmaj,jmax), &
       pia (nmaj,jmax), &
       sion(nmaj,jmax)) 

    allocate &
      (aglw  (nei,nmaj,jmax), &
       photoi(nst,nmaj,jmax), &
       photod(nst,nmaj,jmax), &
       phono(nst,jmax))

    allocate             &
      (sigs(nmaj,nbins), &
       pe  (nmaj,nbins), &
       pin (nmaj,nbins))

    allocate                  &
      (sigex(nei,nmaj,nbins), &
       sigix(nei,nmaj,nbins))

    allocate                  &
      (siga(nei,nbins,nbins), &
       sec (nei,nbins,nbins))

    allocate(iimaxx(nbins))

    allocate           &
      (ww  (nei,nmaj), &
       ao  (nei,nmaj), &
       omeg(nei,nmaj), &
       anu (nei,nmaj), &
       bb  (nei,nmaj), &
       auto(nei,nmaj), &
       thi (nei,nmaj), &
       ak  (nei,nmaj), &
       aj  (nei,nmaj), &
       ts  (nei,nmaj), &
       ta  (nei,nmaj), &
       tb  (nei,nmaj), &
       gams(nei,nmaj), &
       gamb(nei,nmaj))

! Zero all allocated variable arrays:

       zz   (:)     =0.
       zo   (:)     =0.
       zn2  (:)     =0.
       zo2  (:)     =0.
       zno  (:)     =0.
       zns  (:)     =0.
       znd  (:)     =0.
       zrho (:)     =0.
       ze   (:)     =0.
       ztn  (:)     =0.
       zti  (:)     =0.
       zte  (:)     =0.
       eheat(:)     =0.
       tez  (:)     =0.
       tei  (:)     =0.
       tpi  (:)     =0.
       tir  (:)     =0.
       ecalc(:)     =0.
       zxden(:,:)   =0.
       zeta(:,:)    =0.
       zceta(:,:,:) =0.
       vcb(:)       =0.
       zlbh(:,:)    =0.
       phitop(:)    =0.
       ener  (:)    =0.
       del   (:)    =0.
       pflux (:)    =0.
       pener (:)    =0.
       pdel  (:)    =0.
       wave1(:)     =0.
       wave2(:)     =0.
       sflux(:)     =0.
       pespec(:,:)  =0.
       sespec(:,:)  =0.
       uflx  (:,:)  =0.
       dflx  (:,:)  =0.
       zmaj(:,:)    =0.
       zcol(:,:)    =0.
       pia (:,:)    =0.
       sion(:,:)    =0.
       aglw  (:,:,:)=0.
       photoi(:,:,:)=0.
       photod(:,:,:)=0.
       phono(:,:)   =0.
       sigs(:,:)    =0.
       pe  (:,:)    =0.
       pin (:,:)    =0.
       sigex(:,:,:) =0.
       sigix(:,:,:) =0.
       siga(:,:,:)  =0.
       sec (:,:,:)  =0.
       iimaxx(:)    =0.
       ww  (:,:)    =0.
       ao  (:,:)    =0.
       omeg(:,:)    =0.
       anu (:,:)    =0.
       bb  (:,:)    =0.
       auto(:,:)    =0.
       thi (:,:)    =0.
       ak  (:,:)    =0.
       aj  (:,:)    =0.
       ts  (:,:)    =0.
       ta  (:,:)    =0.
       tb  (:,:)    =0.
       gams(:,:)    =0.
       gamb(:,:)    =0.
    end if

  end subroutine cglow_init

!-----------------------------------------------------------------------

end module cglow
