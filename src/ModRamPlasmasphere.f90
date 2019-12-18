MODULE ModRamPlasmasphere
  use ModRamMain, ONLY: DP

  implicit none
  real(DP), parameter :: MHP=1.67e-24, MOP=2.67e-23
  real(DP) :: R, DAY, HOUR

  contains
!==================================================================================================
  subroutine plasmasphere(dt)
     use ModRamParams,    ONLY: PlasmasphereModel
     use ModRamTiming,    ONLY: TimeRamNow
     use ModRamVariables, ONLY: Kp, uL, uT
     use ModRamGrids, ONLY: nR, nT

     real(DP), intent(in) :: dt

     integer  :: i, j, nmonth, iapda, isdate
     real(DP) :: f107_1au, f107_pd, f107_81, f107_365, RZ(3), IG(3), rsn, IAPO(7), &
                 dtL, dtT, dtmin, dts

     ! Need total number of days from TimeRamNow
     DAY = TimeRamNow%iMonth*30 + TimeRamNow%iDay
     HOUR = TimeRamNow%iHour + TimeRamNow%iMinute/60.0 + TimeRamNow%iSecond/3600.0

     ! These read in the ig_rz and apf107 files and compute values needed by iri and msis
     CALL TCON(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,DAY,RZ,IG,rsn,nmonth)
     R = Rz(3)

     select case(trim(PlasmasphereModel))
        case('Carpenter')
           call CARPENTER(Kp) ! Need to calculate last X hour Kp max

        case('Transport')
           ! Get AP and F10.7 information for IRI and MSIS models
           CALL APF_ONLY(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,f107_1au,f107_pd,f107_81,f107_365,iapda,isdate)
           CALL APFMSIS(isdate, HOUR, IAPO)

           ! Compute the ExB velocities
           call velocity

           ! Compute the refilling constant
           call TauFill(IAPO, f107_1au, Kp)

           ! Transport Plasmasphere Flux
           dtL = 10000000.0
           dtT = 10000000.0
           do i = 1, nR
              do j = 1, nT
                 dtL = min(dtL, 1._dp/abs(uL(i,j)))
                 dtT = min(dtT, 1._dp/abs(uT(j,i)))
              enddo
           enddo
           dtmin = 0.75*min(dtL,dtT)
           dtmin = max(dtmin, 1.0)
           dts = 0.0
           do
              dtmin = min(dtmin,dt - dts)
              call transport(dtmin)
              dts = dts + dtmin
              if (dts >= dt) exit
           enddo

        case default
           call CON_STOP('Unrecognized plasmasphere model - '//trim(PlasmasphereModel))

     end select

  end subroutine plasmasphere

!==================================================================================================
  subroutine transport(dt)
     use ModRamGrids, ONLY: nR, nT
     use ModRamVariables, ONLY: uL, uT, necr, tau, Lz
     real(DP), intent(in) :: dt

     integer  :: i, j, N, sgn, nX
     real(DP) :: X, FUP, R, LIMITER, CORR, spp, source, loss

     real(DP), allocatable :: CFL(:), xNE(:), xBND(:)

     nX = max(nR,nT)
     allocate(CFL(0:nX), xBND(0:nX), xNE(0:nX+2))

     !select case(TransportModel)
     !  case('superbee')
     ! Advect Radial
     do j = 1, nT
       CFL(1:nR) = uL(:,j)*dt
       CFL(0)    = CFL(1)
       xNE(1:nR) = necr(:,j)
       xNE(0)    = xNE(1)
       xNE(nR+1) = 1._dp
       xNE(nR+2) = 1._dp
       do i = 1, nR
          sgn = 1
          if (CFL(i) < 0._dp) sgn = -1
          X = xNE(i+1) - xNE(i)
          FUP = 0.5*(xNE(i+1) + xNE(i) - sgn*X)
          if (abs(X) < 1e-15) then
             XBND(i) = FUP
          else
             N = i + 1 - sgn
             R = (xNE(N) - xNE(N-1))/X
             if (R < 0._dp) then
                xBND(i) = FUP
             else
                LIMITER = max(min(2._dp*R,1._dp),min(R,2._dp))
                CORR = 0.5*(CFL(i) - sgn)*X
                xBND(i) = FUP - LIMITER*CORR
             endif
          endif
       enddo
       xBND(0) = xBND(1)
       do i = 1, nR
          necr(i,j) = necr(i,j) - (CFL(i)*xBND(i) - CFL(i-1)*xBND(i-1))
          if (necr(i,j) < 1e-15) necr(i,j) = 1E-15
       enddo
     enddo
   
     ! Advect Azimuthal
     do i = 1, nR
        CFL(1:nT) = uT(:,i)*dt
        CFL(0)    = CFL(nT-1)
        xNE(1:nT) = necr(i,:)
        xNE(0)    = xNE(nT-1)
        xNE(nT+1) = xNE(2)
        xNE(nT+2) = xNE(3)
        do j = 1, nT
           sgn = 1
           if (CFL(j) < 0._dp) sgn = -1
           X = xNE(j+1) - xNE(j)
           FUP = 0.5*(xNE(j+1) + xNE(j) - sgn*X)
           if (abs(X) < 1e-15) then
              XBND(i) = FUP
           else
              N = j + 1 - sgn
              R = (xNE(N) - xNE(N-1))/X
              if (R < 0._dp) then
                 xBND(j) = FUP
              else
                 LIMITER = max(min(2._dp*R,1._dp),min(R,2._dp))
                 CORR = 0.5*(CFL(j) - sgn)*X
                 xBND(j) = FUP - LIMITER*CORR
              endif
           endif
        enddo
        xBND(0) = xBND(nT-1)
        do j = 1, nT
           necr(i,j) = necr(i,j) - (CFL(j)*xBND(j) - CFL(j-1)*xBND(j-1))
           if (necr(i,j) < 1e-15) necr(i,j) = 1E-15
        enddo
     enddo
   
     ! Fill and/or empty
     do i = 1, nR
        spp = saturation(Lz(i))
        do j = 1, nT
           if (necr(i,j) < spp) then
              source = (spp - necr(i,j))/tau(i,j)
              necr(i,j) = necr(i,j) + source*dt
           elseif (necr(i,j) > spp) then
             ! For now we will have a fixed loss rate of 1.5 days
             loss = (spp - necr(i,j))/(1.5*24.0*3600.0)
             necr(i,j) = necr(i,j) - loss*dt
           endif
        enddo
     enddo

     !end select
     deallocate(CFL, xBND, xNE)
  
  end subroutine transport

!==================================================================================================
  subroutine CARPENTER(KpMax)
     use ModRamConst,     ONLY: PI
     use ModRamGrids,     ONLY: nR, nT
     use ModRamVariables, ONLY: Lz, MLT, NECR

     implicit none

     real(DP), intent(in) :: KpMax

     integer,  parameter :: iterMax = 40
     real(DP), parameter :: iterAcc = 0.025

     integer  :: i, j, iter, iErr
     real(DP) :: T, L, C1, C2, X, XN, SPS, PPS, LPPO, LPPI, DX, Y
     real(DP), dimension(2) :: LPP, F
 
     ! Set input dependent variables
     LPPI = 5.6-0.46*KPMAX

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
           SPS = saturation(LPP(i))
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
              SPS = saturation(LPP(2))
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
           SPS = saturation(L)
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

        !a = (1. + 250./RE_km)
        a = 1._dp
        do i = 1, nR
            x = a/Lz(i)
            flux_volume(i,:) = (32./35.) * Lz(i)**4 * sqrt(1. - x) &
                             * (1. + (1./2.)*x + (3./8.)*x**2 + (5./16.)*x**3)
        end do
        
        return
    end subroutine volume_dipole

!============================================================================
    function saturation(L)
      use ModRamConst, ONLY: PI

      real(DP), intent(in) :: L
      real(DP) :: saturation, XT, XD, XR, XL

      XD   = 0.15*(COS(2.*PI*(DAY+9)/365.)-0.5*COS(4.*PI*(DAY+9)/365.))
      XR   = 0.00127*R-0.0635
      XL = -0.3145*L+3.9043
      XT = XL+(XD+XR)*EXP(-(L-2.)/1.5)
      saturation = 10**XT

    end function saturation

!============================================================================
!
!  TauFill()  --  calculate refilling time constants
!
    subroutine TauFill(ap, f10p7, Kp)
        use ModRamParams,    ONLY: TauCalculation
        use ModRamConst,     ONLY: RE_cm
        use ModRamGrids,     ONLY: nR, nT
        use ModRamVariables, ONLY: tau, flux_volume, smLat, smLon, nECR, xRAM, &
                                   yRAM, zRAM, Lz, Phi
        use ModRamTiming,    ONLY: TimeRamNow

        use ModIoUnit,      ONLY: UnitTmp_
        use ModTimeConvert, ONLY: n_day_of_year

        use nrtype, ONLY: pi_d
        implicit none

        real(DP), intent(in) :: ap(7), f10p7, Kp

        real(DP) :: flux, lat, lon, xGSM, yGSM, zGSM, xGEO, yGEO, zGEO, gdLat, gdLon
        integer  :: i, j, iYear, iMonth, iDay, iHour, iMins, iSec, iDoY
        real(DP) :: MinVolume, MinL, FillDays, LFactor, TFactor
        integer :: MinLocation(2)

        call volume_dipole()     ! This initializes flux_volume, change to SCB flux tube volume

        select case(trim(TauCalculation))

           case('Test')
              FillDays = 4.5
              MinVolume = minval(flux_volume(:,:))
              MinLocation = minloc(flux_volume(:,:))
              MinL = Lz(MinLocation(1))
              do i = 1, nR
                 LFactor = 1._dp !(1._dp/Lz(i))**(0.3)
                 do j = 1, nT
                    TFactor = 1._dp !sin(Phi(j))
                    flux = 2._dp*LFactor*TFactor*MinVolume*MinL*saturation(MinL)/(FillDays*24._dp*3600._dp)
                    tau(i,j) = saturation(Lz(i))*flux_volume(i,j)/flux
                 enddo
              enddo

              iYear   = TimeRamNow%iYear
              iMonth  = TimeRamNow%iMonth
              iDay    = TimeRamNow%iDay
              iHour   = TimeRamNow%iHour
              iMins   = TimeRamNow%iMinute
              iSec    = TimeRamNow%iSecond
              iDoY    = n_day_of_year(iYear,iMonth,iDay)
              call RECALC_08(iYear,iDoY,iHour,iMins,iSec,-400._dp,0._dp,0._dp)

              do i = 1, nR

                 do j = 1, nT
                    CALL SMGSW_08(xRAM(1,i,j),yRAM(1,i,j),abs(zRAM(1,i,j)),xGSM,yGSM,zGSM,1)
                    CALL GEOGSW_08(xGEO,yGEO,zGEO,xGSM,yGSM,zGSM,-1)
                    gdLat = pi_d/2.0 - acos(zGEO/sqrt(xGEO**2+yGEO**2+zGEO**2))
                    gdLon = atan2(yGEO,xGEO)
                    call fluxInf(gdLat*180/pi_d,gdLon*180/pi_d,f10p7,ap,flux)
                    tau(i,j) = (saturation(Lz(i))*flux_volume(i,j)*RE_cm/flux)
                 end do

              end do

           case('Fixed')  ! vania addition, Jan 1997
              open(UnitTmp_,file='newtau.dat',status='old')
              do i=1,nR
                 read(UnitTmp_,*) (tau(i,j),j=0,nT)
              enddo
              close(UnitTmp_)
        
           case('Analytic')
              ! Analytical approximation for refilling time constants taken from DGCPM
              FillDays = 4.5
              MinVolume = minval(flux_volume(:,:))
              MinLocation = minloc(flux_volume(:,:))
              MinL = Lz(MinLocation(1))
              do i = 1, nR
                 LFactor = 1._dp !(1._dp/Lz(i))**(0.3)
                 do j = 1, nT
                    TFactor = 1._dp !sin(Phi(j))
                    flux = 2._dp*LFactor*TFactor*MinVolume*MinL*saturation(MinL)/(FillDays*24._dp*3600._dp)
                    tau(i,j) = saturation(Lz(i))*flux_volume(i,j)/flux
                 enddo
              enddo

           case('Rasmussen')
              ! Uses Ionospheric Fluxes to calculate a refilling time constant
              iYear   = TimeRamNow%iYear
              iMonth  = TimeRamNow%iMonth
              iDay    = TimeRamNow%iDay
              iHour   = TimeRamNow%iHour
              iMins   = TimeRamNow%iMinute
              iSec    = TimeRamNow%iSecond
              iDoY    = n_day_of_year(iYear,iMonth,iDay)
              call RECALC_08(iYear,iDoY,iHour,iMins,iSec,-400._dp,0._dp,0._dp)
   
              do i = 1, nR
                 do j = 1, nT
                    CALL SMGSW_08(xRAM(1,i,j),yRAM(1,i,j),abs(zRAM(1,i,j)),xGSM,yGSM,zGSM,1)
                    CALL GEOGSW_08(xGEO,yGEO,zGEO,xGSM,yGSM,zGSM,-1)
                    gdLat = pi_d/2.0 - acos(zGEO/sqrt(xGEO**2+yGEO**2+zGEO**2))
                    gdLon = atan2(yGEO,xGEO)
                    call fluxInf(gdLat*180/pi_d,gdLon*180/pi_d,f10p7,ap,flux) ! This return 'flux' in units of cm^2/s
                    tau(i,j) = saturation(Lz(i))*flux_volume(i,j)*RE_cm/flux
                 end do
              end do

           case default
              call CON_STOP('Unrecognized filling parameter method - '//trim(TauCalculation))

        end select

        return
    end subroutine TauFill

!============================================================================
!
!  fluxInf()  --  calculate maximum ionospheric flux at 1000 km
!
    subroutine fluxInf(lat,lon,f10p7,ap_,flux)
      use ModRamConst,  ONLY: KB, RE_cm, RE_km
      use ModRamTiming, ONLY: TimeRamNow

      use ModTimeConvert, ONLY: n_day_of_year

      implicit none

      real(DP), intent(inout) :: flux
      real(DP), intent(in)    :: lat, lon, f10p7, ap_(7)

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
      JF(34) = .false. ! Turn off messages
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
      alt0 = altr - Hp*log( (7.5e18/(nOp*nO)) * (Ti/Hc)**2 )
      ! Changed 2.0e19 to 7.5e18 as mentioned in Rasmussen et al. 1992 -ME

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
      !flux = flux * ( (RE_km+alt0)/(RE_km+1000.) )**3
      ! Why???? Both Richards and Torr (1985) and Rasmussen et al. (1992) say to
      ! just use flux at alt0 (z0 in their papers). Since the equations
      ! presented here are basically ver-batim from those papers we will not
      ! scale the value -ME

      return
    end subroutine fluxInf

!
!  velocity()  --  get velocity from electric potential - VaniaJ, 2005 - ME, 2019
!
    subroutine velocity()

        use ModRamConst,     ONLY: PI, Re
        use ModRamGrids,     ONLY: nR, nT
        use ModRamTiming,    ONLY: DTs
        use ModRamVariables, ONLY: Lz, RLz, MLT, VT, BNES, MDR, DPHI, uL, uT, EIP

        implicit none

        integer  :: i, j, jp, jm
        real(DP) :: facL, facT, facV, temp

!        A = 7.05e-6 / (1. - 0.159*Kp + 0.0093*Kp*Kp)**3
        facL = 1._dp / RE                         ! per second
        facT = (24.0/(2.*PI)) / RE                ! per second
        facV = 0.5/DPHI/MDR

        do i = 1, nR
           do j = 1, nT
             jp = j+1
             jm = j-1
             if (j.eq.1) jm = nt-1
             if (j.eq.nT) jp = 2
             uL(i,j) = (VT(i,jm)-VT(i,jp))/2/RE/Lz(i)/BNES(i,j)/DPHI*facL
             !uL(i,j) = facV*(VT(I,Jm) - VT(I,Jp))/BNES(I,J)/RLZ(i)
            end do
        end do
        uL(:,1) = uL(:,nT)

        do j = 1, NT
           ut(j,1) = (VT(2,j)-VT(1,j))/Lz(1)/BNES(1,j)/MDR/RE*facT
           !ut(j,1) = 2._dp*facV*(VT(2,j) - VT(1,j))/BNES(1,j)/RLz(i)
           do i = 2, NR
             ut(j,i) = (VT(i+1,j)-VT(i-1,j))/2/Lz(i)/BNES(i,j)/MDR/RE*facT
             !ut(j,i) = facV*(VT(i+1,j) - VT(i-1,j))/BNES(i,j)/RLz(i)
           end do
        end do
        uT(:,:) = uT(:,:) + 7.3E-5/DPHI ! Co-rotation
        uT(1,:) = uT(nT,:)

        return
    end subroutine velocity

END MODULE ModRamPlasmasphere
