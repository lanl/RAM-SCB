MODULE ModRamPlasmasphere
  use ModRamMain, ONLY: DP

  implicit none
  real(DP), parameter :: MHP=1.67e-24, MOP=2.67e-23

  contains
!==================================================================================================
  subroutine CARPENTER(KpMax,day,R)
     use ModRamConst,     ONLY: PI
     use ModRamGrids,     ONLY: nR, nT
     use ModRamVariables, ONLY: Lz, MLT, NECR

     implicit none

     real(DP), intent(in) :: KpMax, R, day

     integer,  parameter :: iterMax = 40
     real(DP), parameter :: iterAcc = 0.025

     integer  :: i, j, iter, iErr
     real(DP) :: T, L, C1, C2, X, XL, XD, XR, XT, XN, SPS, PPS, LPPO, LPPI, DX, Y
     real(DP), dimension(2) :: LPP, F
 


     ! Set input dependent variables
     LPPI = 5.6-0.46*KPMAX
     XD   = 0.15*(COS(2.*PI*(DAY+9)/365.)-0.5*COS(4.*PI*(DAY+9)/365.))
     XR   = 0.00127*R-0.0635

     ! Set initial values
     LPP(1) = LPPI
     LPP(2) = 10.
     F(1) = 0.
     F(2) = 0.

     ! Get MLT dependent variables
     DO j=1,nT
        T=MLT(j)
        IF (T < 6.) THEN
           X  = 0.1
           C1 = 5800.
           C2 = 300.
        ENDIF
        IF ((T >= 6.0).AND.(T <= 19.)) THEN
           X  = 0.1+0.011*(T-6.)
           C1 = -800.
           C2 = 1400.
        ENDIF
        IF ((T > 19.).AND.(T < 20.)) THEN    ! Day-night transition
           X  = 0.1-0.143*(T-20.)
           C1 = 428600.
           C2 = -21200.
        ENDIF
        IF (T >= 20.) THEN
           X  = 0.1
           C1 = -1400.
           C2 = 300.
        ENDIF

        ! Find LPPO for given MLT (between LPP(1) and LPP(2))
        do i = 1, 2
           XL = -0.3145*LPP(i)+3.9043
           XT = XL+(XD+XR)*EXP(-(LPP(i)-2.)/1.5)
           SPS = 10**XT
           PPS = SPS*10.**((LPPI-LPP(i))/X)
           Y = 1. + (C1+C2*T)*LPP(i)**(-4.5) - EXP((2.-LPP(i))/10.)
           F(i) = Y - PPS
        enddo
        if (abs(F(2)) <= 1e-10) then
           LPPO = LPP(2)
        elseif (abs(F(1)) <= 1e-10) then
           LPPO = LPP(1)
        else
           if (F(1) < 0) then
              LPPO = LPP(1)
              DX = LPP(2)-LPP(1)
           else
              LPPO = LPP(2)
              DX = LPP(1)-LPP(2)
           endif
           iter = 0
           do
              iter = iter + 1
              DX = DX*0.5
              LPP(2) = LPPO + DX
              XL = -0.3145*LPP(2)+3.9043
              XT = XL+(XD+XR)*EXP(-(LPP(2)-2.)/1.5)
              SPS = 10**XT
              PPS = SPS*10.**((LPPI-LPP(2))/X)
              Y = 1. + (C1+C2*T)*LPP(2)**(-4.5) - EXP((2.-LPP(2))/10.)
              F(2) = Y - PPS
              if (F(2) <= 0.) LPPO = LPP(2)
              if ((abs(DX) < iterAcc).or.(abs(F(2)) <= 1e-10).or.(iter > iterMax)) exit
           enddo
        endif

        ! Calculate electron density
        DO i=1,nR
           L = Lz(i)
           ! Saturated Plasmasphere Segment
           XL = -0.3145*L+3.9043
           XT = XL+(XD+XR)*EXP(-(L-2.)/1.5)
           SPS = 10**XT
           ! Plasmapause segment
           PPS = SPS*10.**((LPPI-L)/X)
           IF (L <= LPPI) THEN
              NECR(i,j) = SPS
           ELSEIF ((L > LPPI).AND.(L <= LPPO)) THEN
              NECR(I,J) = PPS
           ELSEIF (L > LPPO) THEN
              XN = PPS*(L/LPPO)**(-4.5)
              NECR(i,j) = XN+1.-EXP((2.-L)/10.)
           ENDIF
        ENDDO
     ENDDO

     RETURN

  end subroutine CARPENTER

!============================================================================
!
!  volume()  --  calculate normalized flux-tube volume (X B0/r0)
!
    subroutine volume_dipole()
        use ModRamConst,     ONLY: RE_km
        use ModRamGrids,     ONLY: nR
        use ModRamVariables, ONLY: Lz, flux_volume

        implicit none

        real a, x
        integer i
        save

        a = (1. + 250./RE_km)
        do i = 1, nR
            x = a/Lz(i)
            flux_volume(i,:) = (32./35.) * Lz(i)**4 * sqrt(1. - x) &
                             * (1. + (1./2.)*x + (3./8.)*x**2 + (5./16.)*x**3)
        end do
        
        return
    end subroutine volume_dipole

!============================================================================
!
!  TauFill()  --  calculate refilling time constants (from ionospheric fluxes)
!
    subroutine TauFill(day, R, ap, f10p7, Kp)
        use ModRamParams,    ONLY: useFixedTau
        use ModRamConst,     ONLY: RE_cm
        use ModRamGrids,     ONLY: nR, nT
        use ModRamVariables, ONLY: tau, flux_volume, gdLat, gdLon, nECR
        
        use ModIoUnit, ONLY: UnitTmp_

        implicit none

        real(DP), intent(in) :: day, R, ap, f10p7, Kp

        real(DP) :: flux
        integer  :: i, j

        if(useFixedTau) then  ! vania addition, Jan 1997
           open(UnitTmp_,file='newtau.dat',status='old')
           do i=1,nR
              read(UnitTmp_,*) (tau(i,j),j=0,nT)
           enddo
           close(UnitTmp_)
        else
        call CARPENTER(Kp,day,R) ! This initializes ne 

           do j = 1, nT
              do i = 1, nR
                 call fluxInf(gdLat(i,j),gdLon(i,j),f10p7,ap,R,flux)
                 tau(i,j) = nECR(i,j)*flux_volume(i,j)*RE_cm/flux ! This return flux in units of cm^2/s
              end do
           end do
           tau = tau/3600. ! Convert to hours

        endif

        return
    end subroutine TauFill

!============================================================================
!
!  fluxInf()  --  calculate maximum ionospheric flux at 1000 km
!
    subroutine fluxInf(lat,lon,f10p7,ap_,Rs,flux)
      use ModRamConst,  ONLY: KB, RE_cm, RE_km
      use ModRamTiming, ONLY: TimeRamNow

      use ModTimeConvert, ONLY: n_day_of_year

      implicit none

      real(DP), intent(inout) :: flux
      real(DP), intent(in)    :: lat, lon, f10p7, ap_, Rs

      integer  :: iyd, iMass, iyear, imonth, iday, ihour, imin, isec, doy
      real(DP) :: sec, hours, sLT
      real(DP) :: ap(7), d(9), exoT(2)
      real(DP) :: altr, alt0, Tn, Ti, Te, nO, nH, nOp, rf, g, Hc, Hp, HOp, Ho, r0
 
! New variables for IRI
      real(DP), dimension(1:20,1:1000) :: OUTF
      real(DP), dimension(1:100) :: OARR
      logical, dimension(50) :: JF
      integer :: mmdd

      iYear  = TimeRamNow%iYear
      iMonth = TimeRamNow%iMonth
      iDay   = TimeRamNow%iDay
      iHour  = TimeRamNow%iHour
      iMin   = TimeRamNow%iMinute
      iSec   = TimeRamNow%iSecond
      doy   = n_day_of_year(TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay)
      sec   = iSec + 60._dp*iMin + 3600._dp*iHour
      hours = sec/3600.
      mmdd  = -doy ! Since day is day-of-year, I store this as a negative number in mmdd

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
      JF(1) = .false. ! Don't need Ne
      JF(4:6) = .false.
      JF(21) = .false.
      JF(22) = .false.
      JF(23) = .false.
      JF(28:30) = .false.
      JF(33) = .false.
      JF(35) = .false.


      iyd = 1000*(iyear - 1900) + iday
      altr = 450. + (750. - 450.)*(f10p7 - 70.)/(230. - 70.) ! altr 50 km less than Richards and Torr
      sLT = amod(sec/3600. + lon/15., 24.)
      
      ! Get temperature at altr
      iMass = 0
      call gts7(iyd, sec, altr, lat, lon, sLT, f10p7, f10p7, ap, iMass, d, exoT) 
      Tn = exoT(2)

      ! Get Neutral Oxygen number density at altr
      iMass = 16
      call gts7(iyd, sec, altr, lat, lon, sLT, f10p7, f10p7, ap, iMass, d, exoT)
      nO = d(2) ! O

      ! Get O+ number density and electron and ion temperatures at altr
      call iri_sub(JF, 0, lat, lon, iyear, mmdd, hours, altr, altr, 1, OUTF, OARR)
      nOp = OUTF(5,1)*1e-6 ! Convert from m^(-3) to cm^(-3)
      Ti = OUTF(3,1)
      Te = OUTF(4,1)

      ! Compute new reference altitude
      rf = RE_cm + altr*1.e5
      g = 980.*(RE_cm/rf)**2
      HOp = KB*(Ti + Te)/(MOP*g)
      Ho  = KB*Tn/(MOP*g)
      Hc = HOp*Ho / (HOp - Ho)
      Hp = ( HOp*Ho / (HOp + Ho) ) / 1.e5
      alt0 = altr - Hp*log( (2.e19/(nOp*nO)) * (Ti/Hc)**2 )

      ! Get temperature at alt0
      iMass = 0
      call gts7(iyd, sec, alt0, lat, lon, sLT, f10p7, f10p7, ap, iMass, d, exoT)
      Tn = exoT(2)

      ! Get Hydrogen number density at alt0
      iMass = 1
      call gts7(iyd, sec, alt0, lat, lon, sLT, f10p7, f10p7, ap, iMass, d, exoT)
      nH = d(7) ! H

      ! Get O+ number density and electron and ion temperatures at alt0
      call iri_sub(JF, 0, lat, lon, iyear, mmdd, hours, alt0, alt0, 1, OUTF, OARR)
      nOp = OUTF(5,1)*1e-6  ! Convert from m^(-3) to cm^(-3)
      Ti = OUTF(3,1)
      Te = OUTF(4,1)

      ! Calculate ionospheric flux at alt0
      r0 = RE_cm + alt0*1.e5
      g = 980.*(RE_cm/r0)**2
      HOp = KB*(Ti+Te)/(MOP*g)
      flux = 2.5e-11*sqrt(Tn)*nH*nOp*HOp

      ! Extrapolate flux from alt0 t0 1000 km
      flux = flux * ( (RE_km+alt0)/(RE_km+1000.) )**3

      return
    end subroutine fluxInf

!
!  velocity()  --  get velocity from electric potential - VaniaJ, 2005 - ME, 2019
!
    subroutine velocity()

        use ModRamConst,     ONLY: PI, Re
        use ModRamGrids,     ONLY: nR, nT
        use ModRamTiming,    ONLY: DTs
        use ModRamVariables, ONLY: Lz, MLT, VT, BNES, MDR, DPHI, uL, uT

        implicit none

        integer  :: i, j, j1
        real(DP) :: facL, facT

!        A = 7.05e-6 / (1. - 0.159*Kp + 0.0093*Kp*Kp)**3
        facL = 3600. / RE                                ! per hour
        facT = 3600. * (24.0/(2.*PI)) / RE                ! per hour

        do i = 1, nR+1
           do j = 1, nT
             j1=j
             if (j.eq.0) j1=nt
             uL(i,j) = (VT(i,j1-1)-VT(i,j+1))/2/RE/Lz(i)/BNES(i,j)/DPHI*facL
            end do
        end do

        do j = 1, NT
           do i = 1, NR
             ut(j,i) = (VT(i+1,j)-VT(i-1,j))/2/Lz(i)/BNES(i,j)/MDR/RE*facT
           end do
        end do

        return
    end subroutine velocity

END MODULE ModRamPlasmasphere
