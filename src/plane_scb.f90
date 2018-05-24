!============================================================================
!   plane()  --  calculate flux-tube content in the equatorial plane
!   Copyright (c) 2016, Los Alamos National Security, LLC
!   All rights reserved.
!******************************************************************************

    subroutine plane(yearIn, dayIn, secIn, Kp, ap, R, f10p7, dt, n, wpot)

        use ModPlane
        implicit none

        integer, intent(in) :: yearIn, dayIn
        real(kind=Real8_), intent(in) :: ap, R
        real(kind=Real8_), intent(in) :: secIn, KP, f10p7, dt
        real(kind=Real8_), intent(in) :: wpot(NR+1,NT)
        real(kind=Real8_), intent(inout) :: n(NL,0:*)

        real(kind=Real8_) :: dtH , time, gLong, gLat, mLong0, mLat, update
!        character*1 ans
        logical :: first = .true.
        integer :: IOcount = 0
        integer :: Tcount = 0
        save

        if (first) then
           year = yearIn
           day = dayIn
           sec = secIn
           time = sec/3600.
           update = time
           month = .5 + day/30.
           gLat = 45.                       ! get mLong at 0 mlt and
           gLong = 270.                     ! at initial time (American sector)
           call ggmxx(0, gLong, gLat, mLong0, mLat)
!
!  initialize grid
!
           call initGrid()

           call ngride1(L(1), NL, L(NL+1), L(0), CARTESIAN)
           call ogride1(NL)
           call ngride1(L(1), NL, L(NL+1), L(0), CARTESIAN)
           call ngridp1(mlt, NLT+1)
           call ogridp1(NLT+1)
           call ngridp1(mlt, NLT+1)
!
!  initial data
!
           call initContent(R)
           call volume()
           call TauFill(R, ap, f10p7, mLong0, mLat, Kp) ! get time scales
           update = time + 3.                  ! change fluxes every () hours
           first = .false.
        end if

        if (update .le. time) then
           update = time + 3.                  ! change fluxes every () hours
           call TauFill(R, ap, f10p7, mLong0, mLat, Kp) ! get time scales
        end if

        dtH = dt/2./3600.                             ! also convert to hr.
        call velocity(Kp, wpot)

        call advLtime(time, dtH, n)
        call nTOnsw(n)
        call advTtime(time, dtH)

        time = time + dtH

        call advTtime(time, dtH)
        call nswTOn(n)
        call advLtime(time, dtH, n)

        time = time + dtH

        return
    end subroutine plane

!============================================================================
!
!  advLtime()  --  advance one time step in L direction
!

    subroutine advLtime(time, dt, n)
        use ModPlane, ONLY: L, tau0, n0, uL, NL, NLT, Real8_
        implicit none

        real, intent(in) :: time, dt
        real(kind=Real8_), intent(inout) :: n(NL,0:*)
        real :: S(NL,0:NLT), null(1), sum(NL), bcB, bcE
        integer i, j
        save

        call IonSource(n, S)
        bcB = n0(0)/n0(1)
        do j = 0, NLT
            call veloce1(uL(1,j), NL, uL(NL+1,j), uL(0,j), dt)
            call sourcq1(NL, dt, 3, null, S(1,j), 0.0, 0.0)
            do i = 1, NL
                sum(i) = -4.*n(i,j)*uL(i,j)/L(i)
            end do
            call sourcq1(NL, dt, 3, null, sum, 0.0, 0.0)
            call sourcq1(NL, dt, 2, n(1,j), uL(1,j), uL(NL+1,j), uL(0,j))
            if (uL(NL+1,j) .ge. 0.0) then
                bcE = n0(NL+1)/n0(NL)
                call etbFCT1(n(1,j), n(1,j), NL, bcE, bcB)
            else
                bcE = 0.
                call etbFCT1(n(1,j), n(1,j), NL, bcE, bcB)
                n(NL,j) = 0.
            end if
        end do
         
        return
    end subroutine advLtime

!============================================================================
!
!  advTtime()  --  advance one time step in azimuthal direction
!

    subroutine advTtime(time, dt)
        use ModPlane, ONLY: mlt, nsw, ut, NL, NLT
        implicit none

        real, intent(in) :: time, dt
        integer i
        save

        do i = 1, NL
            call velocp1(ut(0,i), NLT+1, dt)
            call sourcp1(NLT+1, dt, 2, nsw(0,i), ut(0,i))
            call prbfct1(nsw(0,i), nsw(0,i), NLT+1)
        end do

        return
    end subroutine advTtime

!============================================================================
!
!  initContent()  --  initialize tube-content
!
    subroutine initContent(R)
        use ModPlane, ONLY: L, n0, NL, PI, day
        implicit none

        real, intent(in) :: R

        real x, fac, dy, a, Lc, nLc
        integer i   ! ,j
        save

!
!  calculate saturation levels [from, Carpenter and Anderson, 1992]
!
        dy = 2.*PI*(day + 9.)/365.
        fac = 0.15*(cos(dy) - 0.5*cos(2.*dy)) + 0.00127*R - 0.0635
        do i = 0, NL+1
            x = exp(-(L(i) - 2.)/1.5)
            n0(i) = 10**(3.9043 - 0.3145*L(i) + fac*x)
        end do

        a = 3.                                ! change to L^-a dependence
        Lc = 5.2                        ! at L = 5.2
        x = exp(-(Lc - 2.)/1.5)
        nLc = 10**(3.9043 - 0.3145*Lc + fac*x)
        fac = nLc * Lc**a
        do i = 0, NL+1
            if (L(i) .gt. Lc) n0(i) = fac * L(i)**(-a)
        end do
!
        return
    end subroutine initContent

!============================================================================
!
!  volume()  --  calculate normalized flux-tube volume (X B0/r0)
!
    subroutine volume()
        use ModPlane, ONLY: L, vol, NL, RE_km
        implicit none

        real a, x
        integer i
        save

        a = (1. + 250./RE_km)
        do i = 1, NL
            x = a/L(i)
            vol(i) = (32./35.) * L(i)**4 * sqrt(1. - x) &
                     * (1. + (1./2.)*x + (3./8.)*x**2 + (5./16.)*x**3)
        end do
        
        return
    end subroutine volume

!============================================================================
!
!  TauFill()  --  calculate refilling time constants (from ionospheric fluxes)
!
    subroutine TauFill(R, ap, f10p7, mLong, mLat, Kp)
        use ModRamPl_Ne, ONLY: UseSCB_nondipole, useFixedTau, Real8_
        use ModPlane, ONLY: mlt, L, vol, n0, tau0, NL, NLT, RE_cm, day, sec
        use ModIoUnit,      ONLY: UnitTmp_
        implicit none

        real, intent(in) :: R, ap, f10p7, mLong, mLat, Kp
        real :: ns(NL,0:NLT)

        real(kind=Real8_) :: n(NL,0:*)
        real :: S(NL,0:*)
        real :: flux, mLongt, vol_nd, gLong, gLat
        real, external :: flux10p7, volume_nondipole
        integer i, j
        save

        if(useFixedTau) then

!        vania addition, Jan 1997
           open(UnitTmp_,file='newtau.dat',status='old')
           do i=1,nl
              read(UnitTmp_,*) (tau0(i,j),j=0,nlt)
           enddo
           close(UnitTmp_)

        else

!        f10p7 = flux10p7(R)

! This initializes ns
!        call CARPENTER(ns,Kp,day,R)

! This return flux in units of cm^2/s

           do j = 0, NLT

! mlt is in hours from 0 to 24; multiply by 15 to convert to degrees.

              mLongt = mLong - sec/3600.*15. + mlt(j)*15.
              if(mLongt.gt.360) then
                 mLongt = mLongt - 360.
              endif

! convert back to geographic coordinates
              call ggmxx(1, gLong, gLat, mLongt, mLat)

              call fluxInf(gLat,gLong,mLat,mLongt,f10p7,ap,R,flux)

              do i = 1, NL

                 ns(i,j) = n0(i)

! Volume is already normalized by B1
                 if(UseSCB_nondipole) then
                    vol_nd = volume_nondipole(L(i),mlt(j))
                    tau0(i,j) = ns(i,j)*vol_nd*RE_cm/flux
                 else
                    tau0(i,j) = ns(i,j)*vol(i)*RE_cm/flux
                 endif
             
              end do
           end do

! Convert to hours
           tau0 = tau0/3600.

! Try adding minimum time scale; this is a bit of a hack

!       where(tau0<4) tau0 = 4.

        endif

 99        return

      entry IonSource(n, S)
        
        do i = 1, NL
            do j = 0, NLT
                S(i,j) = -(n(i,j) - n0(i)) / tau0(i,j)
!                if (S(i,j) .ne. 0.)        ! speed up refilling near saturation
!     &                    S(i,j) = S(i,j) / sqrt( abs(1. - n(i,j)/n0(i)) )
            end do
        end do

        return
    end subroutine TauFill

!
!  initGrid()  --  initialize grid (L in units of RE)
!  L=1.375, 1.5, 1.75 ... 10, 10.125 & MLT=0, 0.5, 1 ... 24
!
    subroutine initGrid()
        use ModPlane, ONLY: L, mlt, NL, NLT, DL, DLT
        implicit none

        integer i
        save

        L(1) = 1.5
        L(0) = L(1) - 0.5*DL
        do i = 2, NL
            L(i) = L(i-1) + DL
        end do
        L(NL+1) = L(NL) + 0.5*DL
        
        mlt(0) = 0.0
        do i = 1, NLT
            mlt(i) = mlt(i-1) + DLT
        end do

        return
    end subroutine initGrid

!============================================================================
!
!  fluxInf()  --  calculate maximum ionospheric flux at 1000 km
!
    subroutine fluxInf(gLat,gLong,mLat,mLong,f10p7,ap_,Rs,flux)
      use ModPlane, ONLY: RE_cm, RE_km, KB, year, day, month, sec
      use ModRamPl_Ne, ONLY: UseNRL_MSISE, UseIRI_2012
      implicit none

      real, intent(inout) :: flux
      real, intent(in) :: gLat, gLong, mLat, mLong, f10p7, ap_, Rs

! Note that I now define d(9) since NRLMSISE also calcules
!        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
      real :: altr, alt0, glt, glng, sLT, Tn, Ti, Te, nO, nH, nOp, &
           rf, r0, g, HOp, Ho, Hc, Hp, d(9), exoT(2), magbr
      real, dimension(7) :: ap
      real, external :: Rsun
      real :: MHP, MOP
      parameter(MHP=1.67e-24, MOP=2.67e-23)
      integer iyd, iMass, jMag
! ***Chris Jeffery, 6/13*** use magnetic lat long in iri
      data jMag/1/
        save
! New variables for IRI
      real, dimension(1:20,1:1000) :: OUTF
      real, dimension(1:100) :: OARR
      logical, dimension(50) :: JF
      real :: hours,heibeg, heiend, heistp
      integer :: mmdd

!        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!           - ARRAY CONTAINING:
!             (1) DAILY AP
!             (2) 3 HR AP INDEX FOR CURRENT TIME
!             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
!                 PRIOR TO CURRENT TIME
!             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
!                 PRIOR TO CURRENT TIME

! Since I don't have these let's just fill up ap with a constant value

      ap(1:7) = ap_

! Standard IRI 2012 run

      JF(:) = .true.
! Don't need Ne
      JF(1) = .false.
! These are standard
      JF(4:6) = .false.
      JF(21) = .false.
      JF(22) = .false.
      JF(23) = .false.
      JF(28:30) = .false.
      JF(33) = .false.
      JF(35) = .false.

      iyd = 1000*(year - 1900) + day
!                                        altr 50 km less than Richards and Torr
      altr = 450. + (750. - 450.)*(f10p7 - 70.)/(230. - 70.)
      glt = gLat                                        ! degrees
      glng = gLong                                        ! degrees
      sLT = amod(sec/3600. + glng/15., 24.)
! Use the sunspot number from the kp file
!      Rs = Rsun(f10p7)
      
      iMass = 0
      if(UseNRL_MSISE) then
         call gts7(iyd, sec, altr, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)      
      else
         call gts5(iyd, sec, altr, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      endif
      Tn = exoT(2)
      iMass = 16
      if(UseNRL_MSISE) then
         call gts7(iyd, sec, altr, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      else
         call gts5(iyd, sec, altr, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      endif
      nO = d(2)                                                     ! O
      if(UseIRI_2012) then
! Since day is day-of-year, I store this as a negative number in mmdd
         mmdd = -day
         hours = sec/3600.
         heibeg = altr
         heiend = altr
         heistp = 1

         call iri_sub(JF, jMag, mLat, mLong, year, mmdd, hours, heibeg, heiend, heistp, OUTF, OARR)
! Convert from m^(-3) to cm^(-3)
         nOp = OUTF(5,1)*1e-6
         Ti = OUTF(3,1)
         Te = OUTF(4,1)

      else
         call iri(altr, jMag, mLat, mLong, Rs, month, sLT, &
              nOp, Ti, Te, magbr)                                  ! Ti, Te
      endif      
      rf = RE_cm + altr*1.e5
      g = 980.*(RE_cm/rf)**2
      HOp = KB*(Ti + Te)/(MOP*g)
      Ho  = KB*Tn/(MOP*g)
      Hc = HOp*Ho / (HOp - Ho)
      Hp = ( HOp*Ho / (HOp + Ho) ) / 1.e5

      alt0 = altr - Hp*log( (2.e19/(nOp*nO)) * (Ti/Hc)**2 )
      iMass = 0
      if(UseNRL_MSISE) then
         call gts7(iyd, sec, alt0, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      else
         call gts5(iyd, sec, alt0, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      endif
      Tn = exoT(2)
      iMass = 1
      if(UseNRL_MSISE) then
         call gts7(iyd, sec, alt0, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      else
         call gts5(iyd, sec, alt0, glt, glng, sLT, &
              f10p7, f10p7, ap, iMass, d, exoT)
      endif
      nH = d(7)                                                    ! H
      if(UseIRI_2012) then
         heibeg = alt0
         heiend = alt0
         heistp = 1

         call iri_sub(JF, jMag, mLat, mLong, year, mmdd, hours, heibeg, heiend, heistp, OUTF, OARR)
! Convert from m^(-3) to cm^(-3)
         nOp = OUTF(5,1)*1e-6
         Ti = OUTF(3,1)
         Te = OUTF(4,1)

      else
         call iri(alt0, jMag, mLat, mLong, Rs, month, sLT, &
              nOp, Ti, Te, magbr)                                 ! O^+
      endif
      r0 = RE_cm + alt0*1.e5
      g = 980.*(RE_cm/r0)**2
      HOp = KB*(Ti+Te)/(MOP*g)

      flux = 2.5e-11*sqrt(Tn)*nH*nOp*HOp                           ! at alt0
      flux = flux * ( (RE_km+alt0)/(RE_km+1000.) )**3              ! at 1000 km

      return
    end subroutine fluxInf

!============================================================================
!
!  Rsun()  --  calculate sunspot number Rs from F_10.7
!
    real function Rsun(f10p7)
      implicit none

      real, intent(in) :: f10p7
      real :: a, b

      data a/.9261497/, b/57.36373/

      Rsun = (f10p7 - b)/a

      return
    end function Rsun

!============================================================================
!
!  flux10p7()  --  calculate F_10.7 from sunspot number Rs
!
    real function flux10p7(Rs)
      implicit none

      real, intent(in) :: Rs
      real :: a, b
      data a/.9261497/, b/57.36373/

      flux10p7 = a*Rs + b

      return
    end function flux10p7

!============================================================================
!
!  Bfield()  --  calculate Bfield from dipole or nondipole SCB file
!
    real function Bfield(L,mlto)
      use ModRamPl_Ne, ONLY: UseSCB_nondipole, hI_file_name, Nd, NLsh, Nmlt, Npa
      use ModIoUnit,      ONLY: UnitTmp_
      implicit none
      real, intent(in) :: L,mlto

      integer :: i,j,k,iL1,iL2,iT1,iT2
      real :: data_hold(6), Lsh(NLsh), mlt(Nmlt), Bz(NLsh,Nmlt), fluxVolume(NLsh,Nmlt)
      real :: L1,L2,T1,T2,Bz11,Bz12,Bz21,Bz22
      integer, dimension(1) :: temp

      if(UseSCB_nondipole) then

         open(UnitTmp_,file=hI_file_name,status='old')
         read(UnitTmp_,*)
         read(UnitTmp_,*)
         do i=1,NLsh
            do j=1,Nmlt
               read(11,*) Lsh(i), mlt(j), (data_hold(k),k=1,5), Bz(i,j), data_hold(6), fluxVolume(i,j)

               do k=1,Npa-1
                  read(UnitTmp_,*)
               enddo
            enddo
         enddo
         close(UnitTmp_)
! bilinear interpolation
  
         if(L<Lsh(1)) then
            Bfield = 0.3e-4 / L**3                         ! weber/m^2
         else

            temp = minloc( abs(L-Lsh) )
            iL1 = temp(1)
            if(Lsh(iL1)>L) then
               iL2 = iL1-1
            else
               iL2 = iL1+1

               if(iL2.gt.NLsh) then
                  iL2 = iL1-1
               endif
            endif

            temp = minloc( abs(mlt-mlto) )
            iT1 = temp(1)
            if(mlt(iT1)>mlto) then
               iT2 = iT1-1
            else
               iT2 = iT1+1
               if(iT2.gt.Nmlt) then
                  iT2 = 1
               endif
            endif

            L1 = Lsh(iL1)
            L2 = Lsh(iL2)

            T1 = mlt(iT1)
            T2 = mlt(iT2)

            Bz11 = Bz(iL1,iT1)
            Bz12 = Bz(iL1,iT2)
            Bz21 = Bz(iL2,iT1)
            Bz22 = Bz(iL2,iT2)

            Bfield = Bz11*(L2-L)*(T2-mlto) + Bz21*(L-L1)*(T2-mlto) + Bz12*(L2-L)*(mlto-T1) + Bz22*(L-L1)*(mlto-T1)
            Bfield = 1e-9*Bfield/(L2-L1)/(T2-T1)         ! weber/m^2
            
         endif
      else
         Bfield = 0.3e-4 / L**3                         ! weber/m^2
      endif

      return
    end function Bfield

!============================================================================
!
!  volume_nondipole()  --  calculate volume from nondipole SCB file
!
    real function volume_nondipole(Lo,mlto)
      use ModRamPl_Ne, ONLY: UseSCB_nondipole, hI_file_name, Nd, NLsh, Nmlt, Npa
      use ModPlane, ONLY: vol,L
      use ModIoUnit,      ONLY: UnitTmp_
      implicit none
      real, intent(in) :: Lo,mlto

      integer :: i,j,k,iL1,iL2,iT1,iT2
      real :: data_hold(6), Lsh(NLsh), mlt(Nmlt), Bz(NLsh,Nmlt), fluxVolume(NLsh,Nmlt)
      real :: L1,L2,T1,T2,vol11,vol12,vol21,vol22,volave
      integer, dimension(1) :: temp

      open(UnitTmp_,file=hI_file_name,status='old')
      read(UnitTmp_,*)
      read(UnitTmp_,*)
      do i=1,NLsh
         do j=1,Nmlt
            read(UnitTmp_,*) Lsh(i), mlt(j), (data_hold(k),k=1,5), Bz(i,j), data_hold(6), fluxVolume(i,j)

            do k=1,Npa-1
               read(UnitTmp_,*)
            enddo
         enddo
      enddo
      close(11)
! bilinear interpolation
  
      if(Lo<Lsh(1)) then
         
! The -1 is important; L is defined on (0,NL+1)
         temp = minloc( abs(Lo-L) ) - 1
         iL1 = temp(1)

         volume_nondipole = vol(iL1)

      else

         temp = minloc( abs(Lo-Lsh) )
         iL1 = temp(1)
         if(Lsh(iL1)>Lo) then
            iL2 = iL1-1
         else
            iL2 = iL1+1
            
            if(iL2.gt.NLsh) then
               iL2 = iL1-1
            endif
         endif
         
         temp = minloc( abs(mlt-mlto) )
         iT1 = temp(1)
         if(mlt(iT1)>mlto) then
            iT2 = iT1-1
         else
            iT2 = iT1+1
            if(iT2.gt.Nmlt) then
               iT2 = 1
            endif
         endif

         L1 = Lsh(iL1)
         L2 = Lsh(iL2)
         
         T1 = mlt(iT1)
         T2 = mlt(iT2)

         vol11 = fluxVolume(iL1,iT1)
         vol12 = fluxVolume(iL1,iT2)
         vol21 = fluxVolume(iL2,iT1)
         vol22 = fluxVolume(iL2,iT2)

         volume_nondipole = vol11*(L2-Lo)*(T2-mlto) + vol21*(Lo-L1)*(T2-mlto) + vol12*(L2-Lo)*(mlto-T1) + vol22*(Lo-L1)*(mlto-T1)
         volume_nondipole = volume_nondipole/(L2-L1)/(T2-T1)         ! weber/m^2

         
! The -1 is important; L is defined on (0,NL+1)
!         temp = minloc( abs(Lsh(1)-L) ) - 1
!         iL1 = temp(1)

!         volave = sum(fluxVolume(1,:))/Nmlt
!         volave = maxval(fluxVolume(1,:))
! This is a hack
!         volume_nondipole = volume_nondipole*vol(iL1)/volave

         vol11 = fluxVolume(iL1,1)
         vol21 = fluxVolume(iL2,1)
         
         volume_nondipole = (vol11*(L2-Lo) + vol21*(Lo-L1))/(L2-L1)

! This is better
         volume_nondipole = volume_nondipole*30000.      ! divide out equatorial B, 30000 nt at equator

      endif

      return
    end function volume_nondipole

!
!                     
!************************************************************
!*************** EARTH MAGNETIC FIELD ***********************
!************************************************************
!
!
    SUBROUTINE GGMXX(ART,LONG,LATI,MLONG,MLAT)
! CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT)
! FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
! AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
! LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.
      implicit none

      INTEGER, intent(in) :: ART
      REAL, intent(inout) :: MLONG,MLAT,LONG,LATI
        save
 
      real :: CBG, CBM, CI, CLG, CLM, faktor, SBG, SBM, SI, SLG, SLM, YLG, ZPI

      faktor = ATAN(1.0)*4./180.
      ZPI=FAKTOR*360.                    
      CBG=11.4*FAKTOR                              
      CI=COS(CBG)     
      SI=SIN(CBG)
      IF(ART.EQ.0) GOTO 10                         
      CBM=COS(MLAT*FAKTOR)                           
      SBM=SIN(MLAT*FAKTOR)                           
      CLM=COS(MLONG*FAKTOR)                          
      SLM=SIN(MLONG*FAKTOR)
      SBG=SBM*CI-CBM*CLM*SI                     
      LATI=ASIN(SBG)   
      CBG=COS(LATI)     
      SLG=(CBM*SLM)/CBG  
      CLG=(SBM*SI+CBM*CLM*CI)/CBG
        IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)                  
      LONG=ACOS(CLG)  
      IF(SLG.LT.0.0) LONG=ZPI-ACOS(CLG)            
      LATI=LATI/FAKTOR    
      LONG=LONG/FAKTOR  
      LONG=LONG-69.8    
      IF(LONG.LT.0.0) LONG=LONG+360.0                 
      RETURN          
10    YLG=LONG+69.8    
      CBG=COS(LATI*FAKTOR)                           
      SBG=SIN(LATI*FAKTOR)                           
      CLG=COS(YLG*FAKTOR)                          
      SLG=SIN(YLG*FAKTOR)                          
      SBM=SBG*CI+CBG*CLG*SI                        
      MLAT=ASIN(SBM)   
      CBM=COS(MLAT)     
      SLM=(CBG*SLG)/CBM                            
      CLM=(-SBG*SI+CBG*CLG*CI)/CBM
        IF(CLM.GT.1.) CLM=1.                 
      MLONG=ACOS(CLM)  
      IF(SLM.LT..0) MLONG=ZPI-ACOS(CLM)             
      MLAT=MLAT/FAKTOR    
      MLONG=MLONG/FAKTOR  
      RETURN          
    END SUBROUTINE GGMXX

!
!  nTOnsw()  --  copy n(i,j) to nsw(j,i)
!
    subroutine nTOnsw(n)
        use ModPlane, ONLY: nsw, NL, NLT, Real8_
        implicit none

        real(kind=Real8_), intent(in) :: n(NL,0:*)
        integer i, j
        save

        do i = 1, NL
            do j = 0, NLT
                nsw(j,i) = n(i,j)
            end do
        end do

        return
    end subroutine nTOnsw
      
!
!  nswTOn()  --  copy nsw(j,i) to n(i,j)
!
    subroutine nswTOn(n)
        use ModPlane, ONLY: nsw, NL, NLT, Real8_
        implicit none

        real(kind=Real8_), intent(inout) :: n(NL,0:*)
        integer i, j
        save
        
        do i = 1, NL
            do j = 0, NLT
                n(i,j) = nsw(j,i)
            end do
        end do

        return
    end subroutine nswTOn

!
!  velocity()  --  get velocity from electric potential - VaniaJ, 2005
!
    subroutine velocity(Kp, wpot)
        use ModPlane, ONLY: L, mlt, uL, ut, NR, NT, NL, NLT, DTR, PI, RE, Real8_
        implicit none

        real, intent(in) :: KP
        real(kind=Real8_), intent(in) :: wpot(NR+1,NT)

        integer i, j, j1, IER
        real :: A, B, phi, vr, vp, facL, facT, DL1, DPHI,Y
        real, ALLOCATABLE :: RPHI(:), LZ(:), RMLT(:), CWE(:,:)
        real, external :: Bfield
        logical corotation
        data corotation/.true./
        save

        ALLOCATE(LZ(NR+1), RMLT(NT), CWE(0:NL+1,0:NLT),RPHI(NT))

        DL1=(6.5-2.)/(NR-2)
        do i = 1, NR+1
         LZ(i) = 2.+(i-2)*DL1
        end do
        DPHI=2.*PI/(NT-1)      ! Grid size for local time [rad]
        DO J=1,NT
          RPHI(J)=(J-1)*DPHI        ! Magnetic local time in radian
          RMLT(J)=RPHI(J)*12./PI        ! Magnetic local time in hour
        END DO

        do i = 0, NL+1
         do j = 0, NLT
          call ELINTP2(LZ, RMLT, wpot, NR+1, NT, L(i), mlt(j), Y, IER)
          CWE(i,j) = Y*1e3        ! in Volts
!         write (10,15) L(i),mlt(j),CWE(i,j)/1e3
         end do
        end do
        close(10)
15      FORMAT(F5.2,F10.6,E13.4) 

!        A = 7.05e-6 / (1. - 0.159*Kp + 0.0093*Kp*Kp)**3
        facL = 3600. / RE                                ! per hour
        facT = 3600. * (24.0/(2.*PI)) / RE                ! per hour

        do i = 0, NL+1
           do j = 0, NLT-1
             B = Bfield(L(i),mlt(j))                             ! weber/m^2
             j1=j
             if (j.eq.0) j1=nlt
!                uL(i,j) = vr*cos(phi)                         ! dL/dt
!                ut(j,i) = vp*sin(phi)                         ! d(mlt)/dt
             uL(i,j) = (CWE(i,j1-1)-CWE(i,j+1))/2/RE/L(i)/B/DPHI*facL
            end do
            uL(i,nlt) = uL(i,0)
        end do

        do j = 0, NLT
           do i = 1, NL
             B = Bfield(L(i),mlt(j))                             ! weber/m^2
             phi = mlt(j)*15.*DTR
             ut(j,i) = (CWE(i+1,j)-CWE(i-1,j))/2/L(i)/B/DL1/RE*facT
             if (corotation) ut(j,i) = ut(j,i) + 1.
            end do
             B = Bfield(L(0),mlt(1))                             ! weber/m^2
             ut(j,0) = (CWE(1,j)-CWE(0,j))/L(0)/B/DL1/RE*facT
             if (corotation) ut(j,0) = ut(j,0) + 1.
             B = Bfield(L((NL+1)),mlt(1))                           ! weber/m^2
             ut(j,NL+1) = (CWE(NL+1,j)-CWE(NL,j))/L(NL)/B/DL1/RE*facT
             if (corotation) ut(j,NL+1) = ut(j,NL+1) + 1.
        end do

        DEALLOCATE(LZ,RMLT,CWE,RPHI)
        return
    end subroutine velocity

!============================================================================
!
! Routine CARPENTER is a plasmaspheric model of electron densities established
! by Carpenter and Anderson [JGR,1992]
!
    SUBROUTINE CARPENTER(NE,KPMAX,DAY,R)
      use ModPlane, ONLY: NL, NLT, Real8_
      implicit none
      
        real, intent(in) :: KPMAX, R
        real(kind=Real8_), intent(inout) :: NE(NL,0:NLT)
        integer, intent(in) :: DAY

        real :: L, LPPI, LPPO, XN
        real, external :: FUNC, RTBIS, PS, SPS
        
        COMMON /PARAS1/LPPI,X,C1,C2,T1,IDAY,XR
        real :: X, C1, C2, T1, XR
        integer :: IDAY, I, IER, J
        save

        IDAY=DAY
        XR=R
        LPPI=5.6-0.46*KPMAX
        DO 3 J=0,NLT
          T1=24./NLT*J
          IF(T1.LT.6.) THEN
            X=0.1
            C1=5800.
            C2=300.
          ENDIF
          IF(T1.GE.6.0.AND.T1.LE.19.) THEN
            X=0.1+0.011*(T1-6.)
            C1=-800.
            C2=1400.
          ENDIF
          IF(T1.GT.19..AND.T1.LT.20.) THEN    ! Day-night transition
            X=0.1-0.143*(T1-20.)
            C1=428600.
            C2=-21200.
          ENDIF
          IF(T1.GE.20.) THEN
            X=0.1
            C1=-1400.
            C2=300.
          ENDIF
!
! Find Lppo
        LPPO=RTBIS(FUNC,LPPI,10.,0.025,IER)
        DO 3 I=1,NL
          L=1.5+(I-1)*0.25
          IF(L.LE.LPPI) NE(I,J)=SPS(L,DAY,R)
          IF(L.GT.LPPI.AND.L.LE.LPPO) NE(I,J)=PS(L,LPPI,X,DAY,R)
          IF(L.GT.LPPO) THEN
            XN=PS(LPPO,LPPI,X,DAY,R)*(L/LPPO)**(-4.5)
            NE(I,J)=XN+1.-EXP((2.-L)/10.)
          ENDIF
3       CONTINUE
!
        RETURN
    end subroutine CARPENTER

!
! ************************** Functions ***********************************
!
! Function SPS calculate the electron density in the saturated plasmsphere
! segment 
    FUNCTION SPS(L,DAY,R)
      use ModPlane, ONLY: PI
      implicit none

        real, intent(in) :: L, R
        integer, intent(in) :: DAY

        real :: X1, X2, X3, XT, SPS

        X1=-0.3145*L+3.9043
        X2=0.15*(COS(2.*PI*(DAY+9)/365.)-0.5*COS(4.*PI*(DAY+9)/365.))
        X3=0.00127*R-0.0635
        XT=X1+(X2+X3)*EXP(-(L-2.)/1.5)
        SPS=10**XT
        RETURN
    END FUNCTION SPS
!
! Function PS calculate the electron density at the plasmapause segment
    FUNCTION PS(L,LPPI,X,DAY,R)
       implicit none

        real, intent(in) :: L, LPPI, X, R
        INTEGER DAY

        real, external :: SPS
        real :: PS

        PS=SPS(LPPI,DAY,R)*10.**((LPPI-L)/X)
        RETURN
    END FUNCTION PS
!
! Function FUNC is used to locate Lppo
    FUNCTION FUNC(L)
       implicit none

        real, intent(in) :: L

        real :: LPPI, X, C1, C2, T1, R, Y, FUNC
        real, external :: PS
        integer :: DAY
        COMMON /PARAS1/LPPI,X,C1,C2,T1,DAY,R

        Y=(C1+C2*T1)*L**(-4.5)+1.-EXP((2.-L)/10.)
        FUNC=Y-PS(L,LPPI,X,DAY,R)
        RETURN
    END FUNCTION FUNC
!
    FUNCTION RTBIS(FUNC,X1,X2,XACC,IER)
!
!    JMAX...max iterations
!
      implicit none

      real, external :: FUNC
      real, intent(in) :: X1, X2, XACC
      integer, intent(inout) :: IER

      integer :: J, JMAX
      PARAMETER (JMAX=40)

      real :: LPPI, X, C1, C2, T1, DAY, R, Z1, B1, SI1, yy, DX, F, FMID, XMID, RTBIS
      COMMON /PARAS1/LPPI,X,C1,C2,T1,DAY,R
      COMMON /CFUNY/Z1,B1,SI1
      COMMON /Cyy/yy

      IER = 0
      FMID=FUNC(X2)
      F=FUNC(X1)
      IF(FMID.EQ.0) THEN
        RTBIS=X2
        RETURN
      ENDIF        
      IF(F.EQ.0) THEN
        RTBIS=X1
        RETURN
      ENDIF        
!
!    not bracketed
!
      IF(F*FMID.GT.0.) THEN
        IER = 1
      ENDIF
!
!    orient F>0 for X+DX
!
      IF(F.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
!
!    bisection loop
!
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        FMID=FUNC(XMID)
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      IER = 2
      RETURN
    end function RTBIS

!============================================================================
!.........................Interpolation routines - V.J. 2005
!
    SUBROUTINE ELINTP(XX, YY, N, X, Y, IER)
!
      implicit none

      real, intent(in) :: XX(N), YY(N), X
      real, intent(inout) :: Y
      integer, intent(in) :: N
      integer, intent(inout) :: IER
      save

      integer :: JL, JU, JM, J
      real :: D

      IER = 0
!
!    initialize lower and upper values
!
      JL=1
      JU=N
!
!    if not done compute a midpoint
!
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
!
!    now replace lower or upper limit
!
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
!
!    try again
!
      GO TO 10
      ENDIF
!
!    this is J
!
      J=JL        ! If X.LE.XX(1) then J=1
!                  If X.GT.X(J).AND.X.LE.X(J+1) then J=J
!                  If X.GT.X(N) then J=N-1        
      D=XX(J+1)-XX(J)
      Y=(YY(J)*(XX(J+1)-X)+YY(J+1)*(X-XX(J)))/D
      RETURN
    END SUBROUTINE ELINTP

!============================================================================
!
    SUBROUTINE ELINTP2(X1A,X2A,YA,M,N,X1,X2,Y,IER)
!
      use ModPlane, ONLY: Real8_

      implicit none

      real, intent(in) :: X1, X2, Y, X1A(M),X2A(N)
      real(kind=Real8_), intent(in) :: YA(M,N)
      integer, intent(in) :: M, N
      integer, intent(inout) :: IER
      save
      
      real, ALLOCATABLE :: YTMP(:),YYTMP(:)
      integer :: J, K

      ALLOCATE(YTMP(N),YYTMP(M))

      IER = 0
!
!    do M evaluations of row constructed using 1-dimensional evaluator LINTP
!
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL ELINTP(X2A,YTMP,N,X2,YYTMP(J),IER)
        IER = 10 * IER
        IF (IER.EQ.10) RETURN
12    CONTINUE
!
!    evaluate it
!
      CALL ELINTP(X1A,YYTMP,M,X1,Y,IER)
      IER = IER * 10

      DEALLOCATE(YTMP,YYTMP)
      RETURN
   end subroutine ELINTP2
