!******************************************************************************
!    PROGRAM RAM: Ring current - Atmosphere interactions Model
!    Author: Vania K. Jordanova
!    Calculates the evolution of ring current/radiation belt distributions
!    solving the bounce-averaged kinetic equation considering convection, 
!    charge exchange, atmo loss at 200 km, and wave-particle interactions
!    uses a time-dependent convection potential and a non-dipole B field
!    [Jordanova 1995; Jordanova et al., 1996; 2016]
!
!    Numerical scheme for conserv terms: Lax-Wendroff + superbee flux limiter
!       
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************

	SUBROUTINE RAM_ALL (IS)

	use ModRamMain
	use ModTimeConvert
	use ModRamCouple, ONLY: SwmfPot_II
	use ModUtilities, ONLY: flush_unit
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamFunctions
	use ModRamIndices, ONLY: get_indices
	use Module1, ONLY: iDumpRAMFlux
	use ModRamIO, ONLY: gen_old_timetags,write_prefix,UseNewFmt,RamFileName

	implicit none

        real(kind=Real8_) :: XNN(NS,NR),XND(NS,NR),LNCN(NS,NR),LNCD(NS,NR), &
             F2r(NR,NT,NE,NPA),BNESPrev(NR+1,NT),BOUNISPrev(NR+1,NT,NPA), &
             FNISPrev(NR+1,NT,NPA)

	real(kind=Real8_) :: VTOL(NR+1,NT), VTN(NR+1,NT), AVS
 	real(kind=Real8_) :: VT_kV(NR+1,NT)

	CHARACTER*200 ST3, NameFileOut
	CHARACTER*4 ST4 
	CHARACTER*4 ST7no 
	CHARACTER*2 ST8
	CHARACTER*100 inifile
	CHARACTER*200 hIfile
	CHARACTER*100 HEADER
	CHARACTER*2 DayChar, charHr, charEl2 
	CHARACTER*4 charEl, charB

	character(len=200) :: NameEfile
	integer :: i, j, k, l, jw ! Iterators
	integer :: inBlank, iPA
	integer :: nHour, nFive, nFiveDay ! File numbering.
        integer :: nmonth, DAY ! day of year
!       For reading E-field files: UT, radius and phi of current record.
	real(kind=Real8_) :: tolv, RRL, PH, lambda
	type(TimeType) :: TimeNext
!	DT_hI, DtEfi, Dt_bc are timestep for SCB, E field & bound cond update
	integer, intent(in) :: IS ! Species to advect.
	
	SAVE

	T     = UTs
	S     = IS
	DT    = DTs

	if (DoInitOnly) then
	   call ReadPara
	   call Arrays
           IF(DoUseWPI) THEN
	     if (S.EQ.1) then
                CALL WAPARA_HISS	! Hiss diffusion coeff
                IF(DoUseBASdiff) then
                   print*, 'RAM-e: using BAS diff coeffic '
                   CALL WAPARA_BAS	! BAS chorus diff coeffic
                ELSE
                   print*, 'RAM-e: user-supplied diff coeffic '
                   CALL WAPARA_CHORUS	! user-supplied chorus diff coeffic
                ENDIF
             end if
	   ELSE
	     if (S.EQ.1) then
                  print*, 'RAM-e: using electron lifetimes '
	          CALL WAVEPARA1
	          CALL WAVEPARA2
             end if
           ENDIF
	   return
	end if

!	calls for given species S:
	CALL READPARA
	CALL ARRAYS
	CALL CEPARA
	CALL OTHERPARA

	ST3 = trim(PathRamIn)//'/'
	
!	Get timetags for file names.
	call gen_old_timetags(TimeRamStart, TimeRamNow, nHourTotal=nHour, &
          nFiveTotal=nFive, nFiveDay=nFiveDay)

!.......Preparation...CALL THIS ONLY ONCE TO INITIALIZE
	IF (ical.EQ.1) THEN
	   write(st1,'(i3.3)')nHour
	   if (boundary .EQ. 'SWMF') then ! bc files start from 0
	      write(ST7, '(a,a,i4.4)')event, '_', nFive
           else if ((boundary .eq. 'LANL') .or. (boundary .eq. 'PTM')) then
	      write(ST7, '(a,a,I4.4)') event, '_', nFiveDay
	   end if

	   !\
	   ! Electric field and indices.
	   !/
	   ! Create name of first electric field file
	   ! using current time of simulation.
	   call ram_gen_efilename(TimeRamNow, NameEfile)
	   ! Get first electric field values from file.
	   if(electric .ne. 'VOLS') call ram_get_electric(NameEfile, vtol)
	   ! Get indices. 
	   call get_indices(TimeRamNow%Time, Kp, f107)
	   ! Get time 5 mins in future then read next electric field file:
	   TimeNext = TimeRamNow
	   ! Don't iterate file for SCB-traced Efields.
	   if(electric .ne. 'IESC') TimeNext % Time = TimeNext % Time + DtEfi
	   call time_real_to_int(TimeNext)
	   
	   call ram_gen_efilename(TimeNext, NameEfile)
	   if(electric .ne. 'VOLS') call ram_get_electric(NameEfile, vtn)

	   IF (IsComponent.or.electric.EQ.'WESC'.or.electric.eq.'W5SC') THEN
	      VTN  = SwmfPot_II			! convection potential [V]
	      VTOL = SwmfPot_II
	   ENDIF

	   TOLV=T		! Set time of last file read

!.......Zero initial energy losses
	   LSDR(S)=0.
	   LSCHA(S)=0.
	   LSATM(S)=0.
	   LSCOE(S)=0.
	   LSCSC(S)=0.
	   LSWAE(S)=0.
	   ELORC(S)=0.
	   SETRC(S)=0.

!..Read initial files from 3D-Eq code: I(mu), h(mu), Beq(nT), Eind(V/m) & Hdns
	   IF (NameBoundMag .EQ. 'DIPL') THEN
              hIfile=trim(PathScbOut)//'hI_dipole.dat'
	   ELSE
	      if(UseNewFmt) then
		 hIFile=trim(PathScbOut)// &
        RamFileName('hI_output','dat',TimeRamNow)
	      else
                 write(hIFile,'(a,i4.4,a)')trim(PathScbOut)//'hI_output_', &
                      nFive,'.dat'
	      end if
	   END IF

	   OPEN(UNITTMP_,FILE=trim(hIfile),STATUS='OLD')
	   READ(UNITTMP_,'(A)') HEADER
	   READ(UNITTMP_,'(A)') HEADER
	   DO I=2,NR+1
	      DO J=1,NT
		 DO L=1,NPA
		 if (UseEfind) then
		    READ(UNITTMP_,40) RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L), &
           FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J),HDNS(I,J,L),EIR(I,J),EIP(I,J)
		 else
		    READ(UNITTMP_,40) RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L), &
           FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J),HDNS(I,J,L)
		    EIR(I,J)=0.
		    EIP(I,J)=0.
		 endif
		    dIdt(I,J,L)=0.
		    dIbndt(I,J,L)=0.
		 ENDDO
		     IF (J.LT.8.or.J.GT.18) THEN
		      DO L=15,2,-1
			 if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
     			  FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
			 endif
		      ENDDO
		     ENDIF
		 BNES(I,J)=BNES(I,J)/1e9 	! to convert in [T]
		 dBdt(I,J)=0.
	      ENDDO
	   ENDDO
	   CLOSE(UNITTMP_)
	   DO J=1,NT	 	 	 	! use dipole B at I=1
	      BNES(1,J)=0.32/LZ(1)**3/1.e4
	      dBdt(1,J) = 0.
	      EIR(1,J) = 0.
	      EIP(1,J) = 0.
	      DO L=1,NPA
		 FNHS(1,J,L) = FUNT(MU(L))
		 FNIS(1,J,L) = FUNI(MU(L))
		 BOUNHS(1,J,L)=FUNT(COSD(PAbn(L)))
		 BOUNIS(1,J,L)=FUNI(COSD(PAbn(L)))
		 HDNS(1,J,L)=HDNS(2,J,L)
		 dIdt(1,J,L)=0.
		 dIbndt(1,J,L)=0.
	      ENDDO
	   ENDDO
40	   FORMAT(F8.2, F10.1, 3X, I2, 5X, 6(3X, E12.4))

!.......Set up initial and boundary conditions
	   CALL INITIAL(LNCN,LNCD,XNN,XND)
	   CALL GEOSB(nFive)

!.......Read in initial plasmasphere density NECR
           IF(DoUsePlane_SCB) THEN
              OPEN(UNITTMP_,FILE='ne_full.dat',STATUS='OLD') ! Kp=1-2 (quiet)
	      READ(UNITTMP_,'(A)') HEADER
              READ(UNITTMP_,*) ((NECR(I,J),I=1,NL),J=0,NLT)  ! L= 1.5 to 10
	      CLOSE(UNITTMP_)
           END IF

        END IF			!!!! END OF INITIALIZATION !!!!

!.......START THE CALCULATION, write output every DtOutput = 1 hr
	IF(MOD(INT(T),INT(DtOutput)).eq.0.and.S.eq.4) write(ST1,'(I3.3)') nHour
	IF(MOD(INT(T),INT(DtOutput)).EQ.0) CALL WRESULT(LNCN,LNCD,XNN,XND)

        AVS=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE  ! Voll-Stern parameter

	DO I=1,NR+1			! E-field 5 min update
	   DO J=1,NT
	    if (electric .ne. 'VOLS') then
	       VT(I,J)=VTOL(I,J)+(VTN(I,J)-VTOL(I,J))*(T-TOLV)/DtEfi
	    else
	       VT(I,J)=AVS*(LZ(I)*RE)**2*SIN(PHI(J)-PHIOFS)	! [V]
	    endif
	   ENDDO
	ENDDO

! Added capability to couple to plasmaspheric density model
        VT_kV = VT/1e3 ! potential in kV for call to plane

        IF((DoUsePlane_SCB).and.(S.eq.4)) THEN
           ! Need total number of days from TimeRamNow
           DAY = TimeRamNow%iMonth*30 + TimeRamNow%iDay
           CALL APFMSIS(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,TimeRamNow%iHour,IAPO)
           CALL TCON(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,DAY,RZ,IG,rsn,nmonth)
           CALL PLANE(TimeRamNow%iYear,DAY,T,KP,IAPO(2),RZ(3),F107,2.*DT,NECR,VT_kV)
        END IF

! Save flux every 5 min
	DO I = 2, NR
	   DO K = 2, NE
	      DO L = 2, NPA
		 DO J = 1, NT-1
		    FLUX(S,I,J,K,L) = F2(S,I,J,K,L)/FFACTOR(I,K,L)/FNHS(I,J,L)
		 ENDDO
	      ENDDO
	   ENDDO
	ENDDO

! Dump RAM flux to disk in NetCDF format (time unlimited variable)
        IF( (MOD(T,DT_hI) < 1e-3) .and. (DoWriteFlux) ) THEN
           CALL netcdf_flux_dump(FLUX(S,:,:,:,:), NR, LZ, NT, PHI, NE, EKEV, &
                NPA, MU, S, UTs)
        ENDIF

        IF( (MOD(T,DT_hI) < 1e-3) .and. (DoWriteDensity) .and. (S.eq.4) ) THEN
           CALL netcdf_density_dump(NECR, NL, NLT+1, UTs)
        ENDIF

!.......Update boundary conditions and wave diffusion matrix
	   IF (MOD(INT(T),INT(Dt_bc)).EQ.0 .and. ical.GT.1) THEN
	      IF (S.EQ.4) THEN
		 write(ST7no,'(I4.4)') nFive
		 if (boundary .EQ. 'SWMF') ST7 = event//'_'//ST7no
		 if (boundary .EQ. 'LANL') &
		 print*, 'RAM: calling GEOSB at time (hr) = ', T/3600.
	      ENDIF
	      CALL GEOSB(nFive)
              IF (DoUseWPI) THEN
                 IF (S.EQ.1) THEN
                    IF(.not.DoUseBASdiff) then
                       CALL WAPARA_CHORUS	! other chorus diff coeffic
                    ENDIF
                 ENDIF
              ENDIF
	   ENDIF

!.......Update h, I, Beq, Eind, Hdns at every 5 min after initial start
	 IF (NameBoundMag .NE. 'DIPL') THEN
	 IF(MOD(INT(T),INT(DT_hI)).EQ.0 .and. ical.GT.1) THEN ! new 3D-Eq file
	    IF (S.EQ.4) THEN  
	       write(ST4,'(I4.4)') nFive
	       if(UseNewFmt)then
		  hIFile=trim(PathScbOut)//RamFileName('hI_output','dat', &
         TimeRamNow)
		  else
		     hIfile = trim(PathScbOut)//'hI_output_'//ST4//'.dat'
		  end if
	       inBlank = index(hIfile, ' ') - 1
	       OPEN(UNITTMP_,FILE=hIfile(1:inBlank),STATUS='OLD')
	       call flush(6)
	       READ(UNITTMP_,'(A)') HEADER
	       READ(UNITTMP_,'(A)') HEADER
	       DO I=2,NR+1
		  DO J=1,NT
		     BNESPrev(I,J)=BNES(I,J)
		     DO L=1,NPA
			FNISPrev(I,J,L)=FNIS(I,J,L)
			BOUNISPrev(I,J,L)=BOUNIS(I,J,L)
			if (UseEfind) then
			READ(UNITTMP_,40)RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L),&
     	FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J),HDNS(I,J,L),EIR(I,J),EIP(I,J)
			else
			READ(UNITTMP_,40)RRL,PH,IPA,FNHS(I,J,L),BOUNHS(I,J,L),&
     	FNIS(I,J,L),BOUNIS(I,J,L),BNES(I,J),HDNS(I,J,L)
			   EIR(I,J)=0.
			   EIP(I,J)=0.
			endif
			dIdt(I,J,L)=(FNIS(I,J,L)-FNISPrev(I,J,L))/DT_hI
			dIbndt(I,J,L)=(BOUNIS(I,J,L)-BOUNISPrev(I,J,L))/DT_hI
		     ENDDO
		     IF (J.LT.8.or.J.GT.18) THEN
			DO L=15,2,-1
			   if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
			      FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
			   endif
			ENDDO
		     ENDIF
		     BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
		     dBdt(I,J)=(BNES(I,J)-BNESPrev(I,J))/DT_hI
		  ENDDO
	       ENDDO
	       CLOSE(UNITTMP_)
	       DO J=1,NT	! use dipole B at I=1
		  BNES(1,J)=0.32/LZ(1)**3/1.e4
		  dBdt(1,J) = 0.
		  EIR(1,J) = 0.
		  EIP(1,J) = 0.
		  DO L=1,NPA
		     FNHS(1,J,L) = FUNT(MU(L))
		     FNIS(1,J,L) = FUNI(MU(L))
		     BOUNHS(1,J,L)=FUNT(COSD(PAbn(L)))
		     BOUNIS(1,J,L)=FUNI(COSD(PAbn(L)))
		     HDNS(1,J,L)=HDNS(2,J,L)
		     dIdt(1,J,L)=0.
		     dIbndt(1,J,L)=0.
		  ENDDO
	       ENDDO
	    ENDIF
	 ENDIF			! endif read new 3D-Eq file
	ENDIF			! endif dipole case

!.......Update electric field and indices at every DtEfi = 5 min
	 if(MOD(INT(T),INT(DtEfi)).eq.0.and.ical.gt.1) then
	    IF (S.EQ.4) THEN
	    ! Set "new" values of indices and E-field to "old" values. 
	    VTOL = VTN
	    print*,'RAM: updating E field at time ', T/3600.
!       Generate name of next electric field file:
	    TimeNext = TimeRamNow
	    if(electric .ne. 'IESC') then
	       TimeNext % Time = TimeNext % Time + DtEfi
	    else
	       TimeNext % Time = TimeNext % Time
	    end if
	    call time_real_to_int(TimeNext)

	    call ram_gen_efilename(TimeNext, NameEfile)
	    ! Open this file and save contents.
	    if(electric .ne. 'VOLS') call ram_get_electric(NameEfile, vtn)
	    ! No interpolation for IESC case.
	    if(electric .eq. 'IESC') vtol=vtn

!       Set time of "old" file for interpolation:
 1002	    TOLV=T
	    if(IsComponent .OR. electric=='WESC' .or. electric=='W5SC') then
	       VTOL = SwmfPot_II
	       VTN  = SwmfPot_II
	    endif

	  ENDIF
	 endif

!.......Call routines to calculate the changes of distribution function
!.......considering drifts, charge exchange and atmospheric losses
	   CALL DRIFTR
	   CALL DRIFTP
  	   CALL DRIFTE
	   CALL DRIFTMU
	   CALL SUMRC
	   LSDR(S)=LSDR(S)+ELORC(S)
	   if (S.GT.1) then
	     CALL CHAREXCHANGE
 	     CALL SUMRC
	     LSCHA(S)=LSCHA(S)+ELORC(S)
	   else
             IF (DoUseWPI) THEN
                CALL WPADIF	! pitch angle diffusion
             ELSE
	       CALL WAVELO	! electron lifetimes
             ENDIF
	     CALL SUMRC
	     LSWAE(S)=LSWAE(S)+ELORC(S)
	   endif
	   CALL ATMOL
	   CALL SUMRC
	   LSATM(S)=LSATM(S)+ELORC(S)
!.......time splitting
	   CALL ATMOL
	   CALL SUMRC
	   LSATM(S)=LSATM(S)+ELORC(S)
	   if (S.GT.1) then
	     CALL CHAREXCHANGE
	     CALL SUMRC
	     LSCHA(S)=LSCHA(S)+ELORC(S)
	   else
             IF (DoUseWPI) THEN
                CALL WPADIF
             ELSE
	       CALL WAVELO
             ENDIF
	     CALL SUMRC
	     LSWAE(S)=LSWAE(S)+ELORC(S)
	   endif
	   CALL DRIFTMU
	   CALL DRIFTE
	   CALL DRIFTP
	   CALL DRIFTR
	   CALL SUMRC	
	   LSDR(S)=LSDR(S)+ELORC(S)

! Interpolate diffusivity as funtion of Kp, if necessary
           if(DoUseWPI.and.DoUseBASdiff.and.DoUseKpDiff) then
              call WAPARA_Kp()
           end if

	   IF(MOD(INT(T),INT(DT_hI)).EQ.0) CALL ANISCH	! pressure anis calcs

	 RETURN
	 END

 
!**************************  END OF MAIN  ********************************


!*************************************************************************
!                             READPARA
!		Read parameters: Species, Storm event
!*************************************************************************
	SUBROUTINE READPARA

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
 	save

!	Output file names
	ST0='ram'
        IF(S.EQ.1) ST2='_e'
        IF(S.EQ.2) ST2='_h'
        IF(S.EQ.3) ST2='he'
        IF(S.EQ.4) ST2='_o'

	RETURN
	END

!**************************************************************************
!				ARRAYS
!			Set up all the arrays
!**************************************************************************
	SUBROUTINE ARRAYS

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamFunctions
	implicit none
	real(kind=Real8_) :: degrad, camlra, elb, rw, rwu
	real(kind=Real8_) :: clc, spa
	integer :: i, j, k, l, iml, ic, ip
	real(kind=Real8_) :: CONE(NR+4),RLAMBDA(NPA),PA(NPA),MUBOUN
	CHARACTER*80 TITLE
 	save

!.......Grid size of L shell
	DL1 = (RadiusMax - RadiusMin)/(nR - 1)
	IF(MOD(DL1,0.25_8).NE.0) THEN
	  WRITE(6,*) 'RAM: Error : DL is not a multiple of 0.25 '
	  STOP
	END IF

	degrad=pi/180.
	amla(1)=0.		! Magnetic latitude grid in degrees
	DO I=2,6
	 amla(i)=amla(i-1)+0.2
	ENDDO
	DO I=7,24
	 amla(i)=amla(i-1)+0.5
	ENDDO
	DO I=25,Slen
	 amla(i)=amla(i-1)+2.
	ENDDO

	IR1=DL1/0.25
	DR=DL1*RE               ! Grid size for Z=RO  
	DO I=1,NR+1
	 LZ(I)=2.+(I-2)*DL1
	 Z(I)=RE*LZ(I)
	 DO IML=1,Slen
          camlra=amla(iml)*degrad
	  BE(I,IML)=0.32/LZ(I)**3*SQRT(1.+3.*SIN(camlra)**2) &	!in [gauss]
        /COS(camlra)**6
	  ENDDO
	END DO

	DPHI=2.*PI/(NT-1)      ! Grid size for local time [rad]
	IF(MOD(NLT,NT-1).NE.0) THEN
	  WRITE(6,*) ' Error : NT-1 is not a factor of NLT '
	  STOP
	END IF

	DO J=1,NT
	  PHI(J)=(J-1)*DPHI	! Magnetic local time in radian
	  MLT(J)=PHI(J)*12./PI	! Magnetic local time in hour
	END DO
	IP1=(MLT(2)-MLT(1))/0.5
	
	DO I=1,NS
	  RMAS(I)=MP*M1(I)	! rest mass of each species (kg)
	END DO
	
!.......Calculate Kinetic Energy EKEV [keV] at cent, RW depends on NE
	ELB=0.1	 		! Lower limit of energy in keV
	IF(ELB.EQ.0.01) THEN
	  WE(1)=2.8E-3		!  |_._|___.___|____.____|______.______|
	  RW=1.36		!    .     <   DE   >    <      WE     >
	END IF                  !   EKEV                EBND
	IF(ELB.EQ.0.1) THEN	! relativistic
	  WE(1)=3E-2
	  RW=1.27
	END IF
	IF(ELB.EQ.1) THEN
	  WE(1)=0.31
	  RW=1.16
	END IF

	EKEV(1)=ELB+0.5*WE(1)
	GREL(1)=1.+EKEV(1)*1000.*Q/RMAS(S)/CS/CS
	V(1)=CS*SQRT(GREL(1)**2-1.)/GREL(1)
	EBND(1)=ELB+WE(1)
	GRBND(1)=1.+EBND(1)*1000.*Q/RMAS(S)/CS/CS
	VBND(1)=CS*SQRT(GRBND(1)**2-1.)/GRBND(1)
 	DO K=1,NE-1
	  WE(K+1)=WE(K)*RW                   ! WE(K) [keV] is a power series
	  EBND(K+1)=EBND(K)+WE(K+1)          ! E[keV] at bound of grid
	  DE(K)=0.5*(WE(K)+WE(K+1))
	  EKEV(K+1)=EKEV(K)+DE(K)	     ! E[keV] at cent of grid
	  GREL(K+1)=1.+EKEV(K+1)*1000.*Q/RMAS(S)/CS/CS
	  V(K+1)=CS*SQRT(GREL(K+1)**2-1.)/GREL(K+1)   ! Veloc [m/s] at cent
	  GRBND(K+1)=1.+EBND(K+1)*1000.*Q/RMAS(S)/CS/CS
	  VBND(K+1)=CS*SQRT(GRBND(K+1)**2-1.)/GRBND(K+1) ! Veloc [m/s] at bound
	END DO
	DE(NE)=0.5*WE(NE)*(1.+RW)

!.......CONE - pitch angle loss cone in degree
	DO I=1,NR
	  CLC=(RE+HMIN)/Z(I)
	  CONE(I)=ASIND(SQRT(CLC**3/SQRT(4.-3.*CLC)))
	END DO 
	CONE(NR+1)=2.5		! to calcul PA grid near 0 deg
	CONE(NR+2)=1.5
	CONE(NR+3)=1.
	CONE(NR+4)=0.

!.......PA is equatorial pitch angle in deg - PA(1)=90, PA(NPA)=0.
!       MU is cosine of equatorial PA
        PA(1)=90.
	MU(1)=0.
	PA(NPA)=0.
	MU(NPA)=1.
	RWU=0.98
	WMU(1)=(MU(NPA)-MU(1))/32
!		                    ! |_._|___.___|____.____|______.______| 
	DO L=1,46	            !   MU    <  DMU   >    <     WMU     >
	  WMU(L+1)=WMU(L)*RWU
	  DMU(L)=0.5*(WMU(L)+WMU(L+1))
	  MU(L+1)=MU(L)+DMU(L)
	  PA(L+1)=ACOSD(MU(L+1))
	END DO
	PA(48)=18.7
	MU(48)=COSD(PA(48))
	DMU(47)=(MU(48)-MU(47))     
	IC=2
	DO L=48,NPA-1
	  PA(L+1)=CONE(IC)
          IF(L.EQ.49) THEN
	   PA(50)=16.
	  ELSE
	   IC=IC+1
	  ENDIF
	  MU(L+1)=COSD(PA(L+1))
	  DMU(L)=(MU(L+1)-MU(L))       ! Grid size in cos pitch angle
          WMU(L)=2.*(DMU(L-1)-0.5*WMU(L-1))
          IF(L.GT.55)WMU(L)=0.5*(DMU(L)+DMU(L-1))
	END DO
	DMU(NPA)=DMU(NPA-1)
	WMU(NPA)=DMU(NPA-1)
	DO L=1,NPA-1
	 MUBOUN=MU(L)+0.5*WMU(L)           
	 PAbn(L)=ACOSD(MUBOUN)		! PA at boundary of grid
	ENDDO
	PAbn(NPA)=0.
	
	 
!.......Determine the range of NPA such that PA is outside the loss cone:
!	UPA is upper boundary for pitch angle for given Z
	DO I=1,NR
	   UPA(I) = NPA ! SZ, otherwise UPA = 0 for small enough loss cones
	 DO L=NPA,1,-1
	  IF(PA(L).LE.CONE(I)) UPA(I) = L     ! F(UPA)=0. - in loss cone
	 END DO
	END DO

!**  calculate pitch angles for mlat
	DO 15 I=1,NR
	 DO 15 IML=1,Slen
      	  DO 15 IP=1,NPA
           spa=SQRT(SIND(PAbn(ip))**2*BE(i,iml)/BE(i,1))
           IF (spa.GT.1.0) spa=1.0
           ZRpabn(i,ip,iml)=ASIN(spa)
           IF (spa.EQ.1.0) THEN
            ZRpabn(i,ip,iml)=-1.0
	   END IF
15	CONTINUE

!.......FFACTOR is ratio of F2 in conservative space to flux
!	E* are factors to calculate temperature anisotropy
	DO I=1,NR
	 DO K=1,NE
	  DO L=2,NPA				
	   FFACTOR(I,K,L)=LZ(I)*LZ(I)*GREL(K)/SQRT(GREL(K)**2-1.)*MU(L)
	if (ffactor(i,k,l).le.0) print*,'i,k,l,ffactor=',i,k,l,ffactor(i,k,l)
	  ENDDO
	  FFACTOR(I,K,1)=FFACTOR(I,K,2)
	 END DO
	END DO

	DO K=1,NE
	 ERNH(K)=WE(K)*GREL(K)/SQRT((GREL(K)-1.)*(GREL(K)+1.))	! [1/cm3]
	 EPP(K)=ERNH(K)*EKEV(K)
	 FACGR(K)=GREL(K)*SQRT((GREL(K)-1.)*(GREL(K)+1.))
	END DO

!.......to keep F constant at boundary 
	CONF1=((LZ(NR)+DL1)/LZ(NR))**2
	CONF2=((LZ(NR)+2.*DL1)/LZ(NR))**2

3	FORMAT(A)
4	FORMAT(2F7.3,2X,1PE12.5)
	
	RETURN
	END

!**************************************************************************
!			       CEPARA
!       Routine calculates the charge exchange and atmosphere loss rates
!**************************************************************************
	SUBROUTINE CEPARA

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
 	save
	integer :: l, k, i, j
	real(kind=Real8_) :: x, y, alpha, taub
	
	IF (S.EQ.1) GO TO 21
!.......Calculate charge exchange crosss-section of ring species S with H
!         and then the charge exchange decay rate ACHAR
        GO TO (6,7,8) S-1

6	DO 9 L=2,NPA
	DO 9 K=2,NE
	DO 9 I=2,NR
	DO 9 J=1,NT
	     X=LOG10(EKEV(K))
	     IF(X.LT.-2.) X=-2.
	     Y=-18.767-0.11017*X-3.8173e-2*X**2-0.1232*X**3-5.0488e-2*X**4
             ALPHA=10.**Y*V(K)*HDNS(I,J,L)*DT  
!............10**Y is cross section of H+ in m2
9	ACHAR(I,J,K,L)=EXP(-ALPHA)
	GO TO 21

7	DO 11 L=2,NPA
	DO 11 K=2,NE
	DO 11 I=2,NR
	DO 11 J=1,NT
	   X=LOG10(EKEV(K))
	   IF(X.LT.-2.) X=-2.
	   Y=-20.789+0.92316*X-0.68017*X**2+0.66153*X**3-0.20998*X**4
           ALPHA=10.**Y*V(K)*HDNS(I,J,L)*DT	 
!..........10**Y is cross sect of He+ in m2
11	ACHAR(I,J,K,L)=EXP(-ALPHA)
	GO TO 21

8	DO 19 L=2,NPA
	DO 19 K=2,NE
	DO 19 I=2,NR
	DO 19 J=1,NT
	   X=LOG10(EKEV(K))
	   IF(X.LT.-2.) X=-2.
	   Y=-18.987-0.10613*X-5.4841E-3*X**2-1.6262E-2*X**3 &
         -7.0554E-3*X**4
           ALPHA=10.**Y*V(K)*HDNS(I,J,L)*DT	 
!...........10**Y is cross sect of O+ in m2
19	ACHAR(I,J,K,L)=EXP(-ALPHA)

!.......Calculate the losses due to the collis with atmosphere
21	DO 22 K=2,NE
	DO 22 I=2,NR
	TAUB=2*Z(I)/V(K)
22	ATLOS(I,K)=EXP(-DT/TAUB)

	RETURN
	END

!**************************************************************************
!				OTHERPARA
!       		Calculate drift parameters
!**************************************************************************
	SUBROUTINE OTHERPARA

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamFunctions
	implicit none
	save
	real(kind=Real8_) :: RA(NS), MUBOUN, QS
	integer :: i, j, k, l

	DATA RA/1.,.77,.2,.03/    ! Proportions of species in the plasmasph
!.......Electric field offset in radians and particle charge
	 PHIOFS=0*PI/12.
	 QS=1.
	 IF (S.EQ.1) QS=-1.
!...Parameters used in calculating radial and azimuthal drifts at boundaries
	DO I=1,NR
          DO J=1,NT
	    VR(I)=DT/DR/(Z(I)+0.5*DR)/2/DPHI
!...........Kp dependent part of azimuthal drift 
	    P1(I)=DT/DPHI/2/DR/Z(I)
	  END DO				       	
	  DO K=1,NE
	   P2(I,K)=DT*EKEV(K)*1000*(GREL(K)+1)/GREL(K)/Z(I)**2/DPHI/QS
	  END DO				    
	END DO

!.......Pitch angle and energy time derivatives at boundaries of grids
	DO I=1,NR
	  DO L=1,NPA-1
	   MUBOUN=MU(L)+0.5*WMU(L)
	   MUDOT(I,L)=(1.-MUBOUN**2)*DT/2/MUBOUN/Z(I)
	  END DO
	  MUDOT(I,NPA)=0.
	  DO K=1,NE
            EDOT(I,K)=EBND(K)*DT/Z(I)*(GRBND(K)+1)/GRBND(K)/2.
          END DO
	END DO
	
	RETURN
	END

!**************************************************************************
!			       WAVEPARA1
!       	Life time due to WPI inside plasmasphere
!**************************************************************************
	SUBROUTINE WAVEPARA1

	use ModRamMain
	implicit none
	save
	integer :: i, ii, j, k
	real(kind=Real8_):: TAU_WAVE,xE,xL,xlife
	real(kind=Real8_):: rEa(5),rL(8),rlife(5,8),clife(5)
	DATA rEa/0.2,0.5,1.0,1.5,2.0/
	DATA rL/5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.65/

!...... life time (as a function of energy and L)
	rlife(1,1)=6.80
	rlife(1,2)=16.44
	rlife(1,3)=13.75
	rlife(1,4)=17.38
	rlife(1,5)=53.08
	rlife(1,6)=187.06
	rlife(1,7)=93.72
	rlife(1,8)=101571.57
	rlife(2,1)=23.38
	rlife(2,2)=55.98
	rlife(2,3)=43.43
	rlife(2,4)=31.75
	rlife(2,5)=38.20
	rlife(2,6)=104.90
	rlife(2,7)=164.86
	rlife(2,8)=185.67
	rlife(3,1)=343.16
	rlife(3,2)=475.15
	rlife(3,3)=99.87
	rlife(3,4)=62.46
	rlife(3,5)=98.82
	rlife(3,6)=134.95
	rlife(3,7)=171.96
	rlife(3,8)=73.63
	rlife(4,1)=619.62
	rlife(4,2)=356.89
	rlife(4,3)=139.64
	rlife(4,4)=130.32
	rlife(4,5)=210.25
	rlife(4,6)=283.46
	rlife(4,7)=359.03
	rlife(4,8)=159.19
	rlife(5,1)=1062.13
	rlife(5,2)=381.88
	rlife(5,3)=210.37
	rlife(5,4)=231.97
	rlife(5,5)=370.61
	rlife(5,6)=498.14
	rlife(5,7)=638.07
	rlife(5,8)=473.75

!.......Calculates the losses due to the w-p interaction
	DO 22 K=2,NE
	DO 22 II=2,NR

	 xE=EKEV(K)/1000.
	 xL=LZ(II)

	 if(xL.ge.1.65.and.xL.le.5.0) then
	  do i=8,2,-1
	   if(xL.ge.rL(i).and.xL.lt.rL(i-1)) then
	    do j=1,5
	     clife(j)= &
           (log10(rlife(j,i-1))-log10(rlife(j,i)))/(rL(i-1)-rL(i)) &
           *(xL-rL(i))+log10(rlife(j,i))
	     clife(j)=10.**clife(j)
	    enddo
	    goto 822
	   endif
	  enddo
	 elseif(xL.gt.5.0) then
	  do j=1,5
	   clife(j)= &
         (log10(rlife(j,1))-log10(rlife(j,2)))/(rL(1)-rL(2)) &
         *(xL-rL(1))+log10(rlife(j,1))
	   clife(j)=10.**clife(j)
	  enddo
	 elseif(xL.lt.1.65) then
	  do j=1,5
	   clife(j)= &
         (log10(rlife(j,7))-log10(rlife(j,8)))/(rL(7)-rL(8)) &
         *(xL-rL(8))+log10(rlife(j,8))
	   clife(j)=10.**clife(j)
	  enddo
	 endif
822	continue
	if(xE.ge.0.2.and.xE.lt.2.0) then
	 do i=1,4
	  if(xE.ge.rEa(i).and.xE.lt.rEa(i+1)) then
	   xlife=(log10(clife(i+1))-log10(clife(i))) &
         /(log10(rEa(i+1))-log10(rEa(i))) &
         *(log10(xE)-log10(rEa(i)))+log10(clife(i))
	   xlife=10.**xlife
	   goto 823
	  endif
	 enddo 
	elseif(xE.lt.0.2) then
	 xlife=(log10(clife(2))-log10(clife(1))) &
       /(log10(rEa(2))-log10(rEa(1))) &
       *(log10(xE)-log10(rEa(1)))+log10(clife(1))
	 xlife=10.**xlife
	 goto 823
	elseif(xE.ge.2.0) then
	 xlife=(log10(clife(5))-log10(clife(4))) &
       /(log10(rEa(5))-log10(rEa(4))) &
       *(log10(xE)-log10(rEa(5)))+log10(clife(5))
	 xlife=10.**xlife
	 goto 823
	endif
823	continue

	tau_wave=xlife*60.*60.*24.  ! day -> sec

22	WALOS1(II,K)=tau_wave

	RETURN
	END


!**************************************************************************
!			       WAVEPARA2
!       Another life time due to diffusion not everywhere strong
!**************************************************************************
	SUBROUTINE WAVEPARA2

	use ModRamMain
	use ModRamFunctions
	implicit none
	save

	integer :: i, k
	real(kind=Real8_):: TAU_WAVE,EMEV,R1,R2,CONE(NR+4),CLC

	DO 22 K=2,NE
	DO 22 I=2,NR
	 EMEV=EKEV(K)*0.001
	 R1=0.08*EMEV**(-1.32)
	 R2=0.4*10.**(2.*LZ(I)-6.+0.4*log10(29.*EMEV))
	 tau_wave=min(R1,R2)
	 tau_wave=1.0/tau_wave
	 tau_wave=tau_wave*60.*60.*24.  ! day -> sec
	 WALOS2(I,K)=tau_wave
22	CONTINUE

!.......CONE - in degree
	DO I=1,NR
	  CLC=(RE+HMIN)/Z(I)
	  CONE(I)=ASIND(SQRT(CLC**3/SQRT(4.-3.*CLC)))
	END DO 
	CONE(NR+1)=2.5
	CONE(NR+2)=1.5
	CONE(NR+3)=1.
	CONE(NR+4)=0.

	DO I=2,NR
	DO K=2,NE
	 WALOS3(I,K)=64.0*LZ(I)*RE/35./(1-0.25)/ &
       SIN(CONE(I)*PI/180.)/SIN(CONE(I)*PI/180.)/V(K)
	ENDDO
	ENDDO

	RETURN
	END


!*************************************************************************
!                              WAPARA_HISS
!       Routine reading normalized Energy & PA hiss diffusion coeff
!**************************************************************************
        SUBROUTINE WAPARA_HISS

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
	save
	integer :: i, ix, kn, l
        CHARACTER*80 HEADER
	CHARACTER*3 ST4
	CHARACTER*2 ST3

	DO 10 I=1,NR
	 fpofc(1)=2.
	 write(ST4,'(I3.3)') INT(LZ(I)*100)
         DO 15 IX=1,NCF
	 write(ST3,'(I2.2)') INT(fpofc(ix))

         OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/whis_L'//ST4//'_'// &
              ST3//ST2//'.aan',STATUS='old')
	 READ(UNITTMP_,20) HEADER

          DO 100 KN=1,ENG
           read(UNITTMP_,17) ENOR(KN)
           read(UNITTMP_,27)(ndaaj(i,kn,l,ix),l=1,npa)
100       CONTINUE
          ndvvj = 0
	  IF (IX.LT.NCF) fpofc(ix+1)=fpofc(ix)+4.

          DO KN=1,ENG
          DO L=1,NPA
	   if (ndaaj(i,kn,l,ix).lt.1e-20) ndaaj(i,kn,l,ix)=1e-20
	  ENDDO
	  ENDDO
15	 CONTINUE
10      CONTINUE

17	FORMAT(E13.4)
20	FORMAT(A80)
27	FORMAT(80(1PE12.3))
	CLOSE(UNITTMP_)

        RETURN
        END


!*************************************************************************
!                              WAPARA_CHORUS
!       Routine reading bounce-aver PA wave diffusion coeff
!**************************************************************************
        SUBROUTINE WAPARA_CHORUS

	use ModRamMain
	implicit none
	save
	integer :: i, j, kn, l, ikp
	real(kind=Real8_):: RLDAA(ENG,NPA),RUDAA(ENG,NPA)
	CHARACTER*1 ST3
	CHARACTER*80 HEADER

	ikp=INT(KP)
	IF (ikp.gt.4) ikp=4
	write(ST3,'(I1.1)') ikp

	OPEN(UNIT=22,FILE=trim(PathRamIn)//'/wlowcho_Kp'//ST3// &
          ST2//'.aan',STATUS='old')
	OPEN(UNIT=23,FILE=trim(PathRamIn)//'/wuppcho_Kp'//ST3// &
          ST2//'.aan',STATUS='old')
	 READ(22,20) HEADER
	 print*,'in WAPARA_CHORUS: ',trim(HEADER)
	 READ(23,20) HEADER

	DO 10 I=1,NR
	 DO 10 J=1,NT
	 READ(22,20) HEADER
	 READ(23,20) HEADER

          DO 15 KN=1,ENG
           read(22,17) ECHOR(KN)
           read(23,17) ECHOR(KN)
           read(22,27)(RLDAA(kn,l),l=1,npa)
           read(23,27)(RUDAA(kn,l),l=1,npa)
 15	  CONTINUE
	  
          DO KN=1,ENG
          DO L=1,NPA
	   if (RLDAA(kn,l).lt.1e-20) RLDAA(kn,l)=1e-20
	   if (RUDAA(kn,l).lt.1e-20) RUDAA(kn,l)=1e-20
	   BDAAR(i,j,kn,l)=(RLDAA(kn,l)+RUDAA(kn,l))

	  ENDDO
	  ENDDO
10      CONTINUE

17	FORMAT(E13.4)
20	FORMAT(A80)
27	FORMAT(80(1PE12.3))
	CLOSE(22)
	CLOSE(23)

        RETURN
        END


! *************************************************************************
!                              WAPARA_Kp
!       Interpolate Kp-dependent diffusion cofficients
!**************************************************************************
        SUBROUTINE WAPARA_Kp()

	use ModRamMain
	implicit none

        integer :: i1,i2

        if(Kp.gt.maxval(Kp_chorus)) then
           CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,NKpDiff)

        else
           i1 = minloc(abs(Kp-Kp_chorus),dim=1)
           if(Kp.lt.Kp_chorus(i1)) then
              i1 = i1-1
           end if
           i2 = i1+1

           ! Linear interpolation of Kp
           CDAAR(:,:,:,:) = (Kp-Kp_chorus(i1))*CDAAR_chorus(:,:,1:35,:,i2) + &
               (Kp_chorus(i2)-Kp)*CDAAR_chorus(:,:,1:35,:,i1)
           CDAAR = CDAAR/(Kp_chorus(i2) - Kp_chorus(i1))

        end if

        CDAAR(:,25,:,:) = CDAAR(:,1,:,:)

        return
        end SUBROUTINE WAPARA_Kp


! *************************************************************************
!                              WAPARA_BAS
!       Routine reading normalized Energy & PA wave diffusion coeff
!**************************************************************************
        SUBROUTINE WAPARA_BAS

        use ModRamMain, only : NR_Dxx,NT_Dxx,NE_Dxx,NPA_Dxx,LZ,EKEV, &
             CDAAR,CDAAR_chorus,DoUseKpDiff,NKpDiff,NR,NT,NE,ENG,NPA, &
             RCHOR_Dxx,TCHOR_Dxx,ECHOR_Dxx,PACHOR_Dxx,BDAAR,MU,Real8_
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamFunctions
	implicit none

        integer :: i,j,k,l,IER,nkp,nloop
        CHARACTER*32 :: H1,H2,H3,nchar
        CHARACTER*64 :: fname
        real(kind=Real8_) :: Dxx_hold(NR_Dxx,NE_Dxx,NPA_Dxx)
        real(kind=Real8_) :: PA(NPA)

        write(*,*) "Starting WAPARA_BAS"

	DO L=1,NPA
           PA(L)=ACOSD(MU(L))
        END DO

        if(DoUseKpDiff) then
           nloop = NKpDiff
        else
           nloop = 1
        endif

        do nkp=1,nloop

           write(nchar,'(i1)') nkp-1
           fname = 'bav_diffcoef_chorus_rpa_Kp'//trim(nchar)//'.PAonly.dat'
           OPEN(UNIT=UNITTMP_,FILE=fname,STATUS='old')

           ! First skip over header
           do i=1,12,1
              read(UNITTMP_,*)
           end do

           do i=1,NR_Dxx
              do j=1,NT_Dxx
                 do k=1,NE_Dxx
                    read(UNITTMP_,'(A19,F6.4,A9,F6.4,A12,F8.2)') H1,  &
                    RCHOR_Dxx(i), H2, TCHOR_Dxx(j), H3, ECHOR_Dxx(k)
                    
                    do l=1,NPA_Dxx
                       read(UNITTMP_,'(F15.6,3E18.6)') PACHOR_Dxx(l),  &
                       CDAAR_chorus(i,j,k,l,nkp)
                    end do
                    
                    ! skip over blank lines
                    do l=1,4
                       read(UNITTMP_,*)
                    end do
                    
                 end do
              end do
           end do
           
           CLOSE(UNIT=UNITTMP_)

! Interpolate onto L-shell, energy and pitch angle, assuming time is normal

           if((NT.ne.25).and.(NE.ne.35)) then
              write(*,*) "We have a problem... assuming NT=25 & NE=35"
              stop
           end if

           if((NR_Dxx.eq.NR).and.(NPA.eq.NPA_Dxx)) then
            write(*,*) "No interpolation of diffusion coeff for", fname
           end if

        end do   ! end NKpDiff loop

! Initialization
        CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,1)

        write(*,*) "Finished WAPARA_BAS"

        RETURN
        END


!**************************************************************************
!				 INITIAL
!   		Initial set up of distribution function
!**************************************************************************
	SUBROUTINE INITIAL(LNCN,LNCD,XNN,XND)

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamIO,     ONLY: write_prefix
	use ModRamFunctions
	
	implicit none
	save
	integer :: i, j, k, l, inBlank, ik, ier, ne1
	real(kind=Real8_):: esum, weight, e1, y, yz, y10
	real(kind=Real8_):: XNN(NS,NR),XND(NS,NR),LNCN(NS,NR),LNCD(NS,NR),N, &
         F2r(NR,NT,36,NPA),FI(20,93),LI(20),EI(93),LIN(13),EIN(4),NI(13,4)
	CHARACTER*80 HEADER
	CHARACTER*200 ST4
	CHARACTER*8 DU

	ST4 = trim(PathRamIn)
	inBlank = index(ST4, ' ') - 1

	N=0				! total # dens
	ESUM=0				! total energy
	DO I=1,NR
	  ENERN(S,I)=0.
	  ENERD(S,I)=0.
	END DO

	if (.not.DoUseVAPini) go to 25	! quiet time initial conditions

!......Initial flux Frbsp*[1/cm2/s/sr/keV] from RBSP data
	if(.not.IsRestart)then
	 call write_prefix
	 write(*,*)'RBSP initial conditions for species ', st2
	 OPEN(UNIT=UNITTMP_, &
	  FILE=ST4(1:inBlank)//'Frbspmar13'//ST2//'.in',STATUS='OLD')
	  READ(UNITTMP_,10) HEADER
10	  FORMAT(A80)
          ne1=91
          if (s.gt.2) ne1=72
	  READ(UNITTMP_,*) DU,(EI(K),K=1,ne1)		! energy grid in (keV)
	  DO I=1,20
	   READ(UNITTMP_,*) LI(I),(FI(I,K),K=1,ne1)	! input is in log(FI)
	  ENDDO
	 CLOSE(UNITTMP_)          
	 DO K=1,ne1
	   EI(K)=LOG10(EI(K))
	 ENDDO
!.......Read the quiet time PA anisotropy of [Garcia & Spjeldvik,1985]
         OPEN(UNIT=UNITTMP_,FILE=ST4(1:inBlank)//'sspd85.in',STATUS='OLD')
          READ(UNITTMP_,10) HEADER
          READ(UNITTMP_,*) (LIN(I),I=1,13)
          READ(UNITTMP_,10) HEADER
          READ(UNITTMP_,*) (EIN(K),K=1,4)
          READ(UNITTMP_,10) HEADER
          DO I=1,13
           READ(UNITTMP_,*) (NI(I,K),K=1,4)
          ENDDO
         CLOSE(UNITTMP_)
!.......Two-dimensional interpolation to get F at midnight
	 DO I=1,NR                               ! F2(R=1)=0 - bound cond
	  DO K=NE,1,-1
	    E1=LOG10(EKEV(K))
	    CALL LINTP2(LI,EI,FI,20,ne1,LZ(I),E1,Y,IER)
	    CALL LINTP2(LIN,EIN,NI,13,4,LZ(I),EKEV(K),YZ,IER)
!...........where LINTP2 is a 2-D interpolation routine in RECIPES2.F
	     Y10=10.**Y                         ! flux input
	     DO L=1,NPA				! F2 in conserv space
	       F2(S,I,1,K,L)=Y10*(1.-MU(L)**2)**(YZ/2.)*FFACTOR(I,K,L)
	     END DO
	  END DO
	 END DO
!.......symmetric initial RC
	 DO 16 I=1,NR
	 DO 16 K=1,NE
	 DO 16 L=1,NPA
	  DO J=2,NT
	   F2(S,I,J,K,L)=F2(S,I,1,K,L)
	  ENDDO
16	 CONTINUE   
         go to 40
	end if

!.......read the function F2 from initial files
25	if(.not.IsRestart)then
	   call write_prefix
	   write(*,*)'Loading initial conditions for species ', st2
	   if (nT==25) then
	      open(unit=UNITTMP_, &
            file=ST4(1:inBlank)//'f2inirc'//ST2//'.dat',status='old') 
	   else 
	      call CON_stop('RAM: f2ini file does not match this grid.')
	   end if
		DO I=1,20
		   DO J=1,NT
		      DO IK=1,36
			 READ (UNITTMP_,30) (F2r(I,J,IK,L),L=1,NPA)
		      ENDDO
		      DO K=1,NE
			 DO L=1,NPA
			    F2(S,I,J,K,L)=F2r(I,J,K,L)*FNHS(I,J,L)/FUNT(MU(L))
			 END DO
		      ENDDO
		   ENDDO
		ENDDO
 30		FORMAT(75(1PE11.3))
		close(UNITTMP_)
	     end if

!.......Calculate the total number of particles and energy of this specie
40	DO I=2,NR
	   DO J=2,NT
	      DO K=2,NE
		DO L=2,UPA(I)-1
	         WEIGHT=F2(S,I,J,K,L)*WE(K)*WMU(L)
	         IF(MLT(J).LE.6.OR.MLT(J).GE.18.) THEN
		  XNN(S,I)=XNN(S,I)+WEIGHT		     	
	          ENERN(S,I)=ENERN(S,I)+EKEV(K)*WEIGHT
		 ELSE
		  XND(S,I)=XND(S,I)+WEIGHT		     	
	          ENERD(S,I)=ENERD(S,I)+EKEV(K)*WEIGHT
		 ENDIF		  
		END DO
	      ENDDO
	   ENDDO
	   ESUM=ESUM+ENERN(S,I)+ENERD(S,I)
	   N=N+XNN(S,I)+XND(S,I)
	ENDDO

	FACTOR=3.4027E10*DR*DPHI

!.......Initial loss is zero
	DO I=1,NR
	 LNCN(S,I)=0.
	 LNCD(S,I)=0.
	 LECN(S,I)=0.
	 LECD(S,I)=0.
	END DO

	RETURN
	END

!**************************************************************************
!				WRESULT
!       	Routine printing all results at time T
!**************************************************************************
	SUBROUTINE WRESULT(LNCN,LNCD,XNN,XND)

	use ModRamMain
	use ModIoUnit, ONLY: UNITTMP_, io_unit_new
	use ModRamIO,  ONLY: UseNewFmt, RamFileName
	use ModRamFunctions

	implicit none
	save
	integer :: i, j, k, l, jw, iw
	real(kind=Real8_):: weight, esum, csum, psum, precfl
        real(kind=Real8_)::LNCN(NS,NR),LNCD(NS,NR),XNN(NS,NR),XND(NS,NR), &
             XNNO(NR),XNDO(NR)
     	real(kind=Real8_)::F(NR,NT,NE,NPA),NSUM,FZERO(NR,NT,NE), &
           ENO(NR),EDO(NR),AVEFL(NR,NT,NE),BARFL(NE)
	character(len=23) :: StringDate
	CHARACTER*4 ST3
	CHARACTER*100 ST4, NameFileOut

	ST4=trim(PathRamOut)
        write(StringDate,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2)'.',i3.3)") &
        TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
        TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond, &
        floor(TimeRamNow%FracSecond*1000.0)

!.......Print RAM output at every hour

	DO I=2,NR
	 XNNO(I)=XNN(S,I)
	 XNDO(I)=XND(S,I)
	 XNN(S,I)=0.
	 XND(S,I)=0.
	 ENO(I)=ENERN(S,I)
	 EDO(I)=ENERD(S,I)
	 ENERN(S,I)=0.
	 ENERD(S,I)=0.
	 DO K=2,NE
	  DO L=1,NPA
	   DO J=1,NT-1
	    IF (L.EQ.1) THEN
	     F(I,J,K,1)=F2(S,I,J,K,2)/FFACTOR(I,K,2)/FNHS(I,J,2)
	    ELSE
	     F(I,J,K,L)=F2(S,I,J,K,L)/FFACTOR(I,K,L)/FNHS(I,J,L)
	    ENDIF
	    if(f(i,j,k,l).lt.0) then
	       iUnit7 = io_unit_new()
	     open(UNIT=iUnit7,FILE=trim(ST4)//'fch'//ST0//ST1//ST2//'.dat', &
           STATUS='UNKNOWN',ACCESS='APPEND')
	     write(iUnit7,*)' f=',f(i,j,k,l),' i=',i,' j=',j,' k=',k,' l=',l, &
           '  T=',T/3600.
!       '
	     f2(S,i,j,k,l)=1E-15
	     f(i,j,k,l)=1E-15
	     close(iUnit7)
            endif
	    if(f2(S,i,j,k,l).lt.1E-5) f2(S,i,j,k,l)=1E-5
	    if(f(i,j,k,l).lt.1E-5) f(i,j,k,l)=1E-5
	    IF(L.LT.UPA(I)) THEN
	     WEIGHT=F2(S,I,J,K,L)*WE(K)*WMU(L)
	     IF(MLT(J).LE.6.OR.MLT(J).GE.18.) THEN
	      XNN(S,I)=XNN(S,I)+WEIGHT
	      ENERN(S,I)=ENERN(S,I)+EKEV(K)*WEIGHT
	     ELSE
	      XND(S,I)=XND(S,I)+WEIGHT
	      ENERD(S,I)=ENERD(S,I)+EKEV(K)*WEIGHT
	     ENDIF	      
	    ENDIF
	   END DO
	   F(I,NT,K,L)=F(I,1,K,L)
	  END DO
	 END DO
	 LNCN(S,I)=XNNO(I)-XNN(S,I)
	 LECN(S,I)=ENO(I)-ENERN(S,I)
	 LNCD(S,I)=XNDO(I)-XND(S,I)
	 LECD(S,I)=EDO(I)-ENERD(S,I)
	END DO

!.......Write the trapped equatorial flux [1/s/cm2/sr/keV]
	if (NT.EQ.49) JW=6
	if (NT.EQ.25) JW=3
        IW=4
	if(UseNewFmt) then
	   NameFileOut=trim(PathRamOut)//RamFileName('ram'//st2,'t',TimeRamNow)
	else
	   NameFileOut=trim(ST4)//ST0//ST1//ST2//'.t'
	end if
	open(unit=UnitTMP_, file=trim(NameFileOut), status='UNKNOWN')
	
        DO I=4,NR,IW
        DO 25 J=1,NT-1,JW
	   WRITE(UNITTMP_,32) StringDate,LZ(I),KP,MLT(J)
	   DO 27 K=4,NE-1
27	     WRITE(UNITTMP_,30) EKEV(K),(F(I,J,K,L),L=2,NPA-2)
25       CONTINUE
	END DO
	close(UNITTMP_)

30      FORMAT(F7.2,72(1PE11.3))
32	FORMAT(' EKEV/PA, Date=',a,' L=',F6.2,' Kp=',F6.2,' MLT=',F4.1)
40	FORMAT(20(3X,F8.2))

!.......Write the total precipitating flux [1/cm2/s]
        IF(DoUseWPI) THEN
        OPEN(UNIT=UNITTMP_,FILE=trim(ST4)// &
            ST0//ST1//ST2//'.tip',STATUS='UNKNOWN')
	   WRITE(UNITTMP_,71) T/3600,KP
	   DO I=2,NR
	      DO J=1,NT
		 PRECFL=0.
		 DO K=2,NE	! 0.15 - 430 keV
		    AVEFL(I,J,K)=0.
		    DO L=UPA(I),NPA
		       AVEFL(I,J,K)=AVEFL(I,J,K)+F(I,J,K,L)*WMU(L)
		    ENDDO
		    AVEFL(I,J,K)=AVEFL(I,J,K)/(MU(NPA)-MU(UPA(I)))
		    PRECFL=PRECFL+AVEFL(I,J,K)*PI*WE(K)
		 ENDDO
		 WRITE(UNITTMP_,70) LZ(I),PHI(J),PRECFL
	      END DO
	   END DO
	close(UNITTMP_)
        ENDIF

70	FORMAT(F5.2,F10.6,E13.4) 	
71      FORMAT(2X,3HT =,F8.0,2X,4HKp =,F6.2,2X,'  Total Precip Flux [1/cm2/s]')

!.......Write the electric potential [kV]
	   if(UseNewFmt)then
	    NameFileOut=trim(PathRamOut)//RamFileName('efield','in',TimeRamNow)
	   else
	    NameFileOut=trim(PathRamOut)//'efield_'//ST1//'.in'
	   end if
	   OPEN(UNIT=UNITTMP_,FILE=NameFileOut,STATUS='UNKNOWN')
     	   WRITE(UNITTMP_,*)'UT(hr)    Kp   F107   VTmax   VTmin   Date=', &
             StringDate
           WRITE(UNITTMP_,31) T/3600.,KP,F107,maxval(VT)/1e3,minval(VT)/1e3
	   WRITE (UNITTMP_,*) ' L     MLT        Epot(kV)'
	    DO I=1,NR+1
	    DO J=1,NT
              WRITE (UNITTMP_,22) LZ(I),MLT(J),VT(I,J)/1e3
	    END DO
	    END DO
22	   FORMAT(F5.2,F10.6,E13.4)
31	   FORMAT(F6.2,1X,F6.2,2X,F6.2,1X,F7.2,1X,F7.2,1X,1PE11.3)
	   CLOSE(UNITTMP_)

!.......Write the plasmaspheric electron density [cm-3]
          IF(DoUsePlane_SCB) THEN
	  if(UseNewFmt)then
	   NameFileOut=trim(PathRamOut)//RamFileName('plasmne','in',TimeRamNow)
	  else
	   NameFileOut=trim(PathRamOut)//'plasmne_'//ST1//'.in'
	  end if
	   OPEN(UNIT=UNITTMP_,FILE=NameFileOut,STATUS='UNKNOWN')
           WRITE(UNITTMP_,96) T/3600,KP,IAPO(2),RZ(3),StringDate
	   DO  I=1,NR+1
	     K=(I-2)*IR1+3
	   DO  J=1,NT
	     L=(J-1)*IP1	  
             WRITE(UNITTMP_,70) LZ(I),PHI(J),NECR(K,L)
	   ENDDO
	   ENDDO
96	   FORMAT(2HT=,F6.2,4H Kp=,F5.2,4H AP=,F7.2,4H Rs=,F7.2,7H Date= ,A23,&
           ' Plasmasphere e- density [cm-3]')
	   CLOSE(UNITTMP_)
           ENDIF

!.......Write the RAM flux [1/s/cm2/sr/keV] at given pitch angle
           IF(DoUseWPI) THEN
           OPEN(UNIT=UNITTMP_,FILE=trim(PathRamOut)//'outm'//ST0//ST2//ST1 &
           //'.dat',STATUS='UNKNOWN')
 	   DO I=1,NR
	   DO 822 J=1,NT-1
	   DO 822 K=2,NE
	    WRITE(UNITTMP_,31) T/3600.,LZ(I),KP,MLT(J),EKEV(K),F(I,J,K,27)
822	   CONTINUE
	   ENDDO
	   CLOSE(UNITTMP_)
           ENDIF

        RETURN
        END

!************************************************************************
!			SUMRC     
!       	Calculate the total energy & energy loss 
!************************************************************************
	SUBROUTINE SUMRC

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_

	implicit none
	save
	integer :: i, j, k, l
	real(kind=Real8_):: enold
	real(kind=Real8_):: WEIGHT

	ELORC(S)=0.
	ENOLD=SETRC(S)
	SETRC(S)=0.
	DO I=2,NR
	 DO K=2,NE
	  DO L=2,NPA
	   DO J=1,NT-1
	    WEIGHT=F2(S,I,J,K,L)*WE(K)*WMU(L)
	    SETRC(S)=SETRC(S)+EKEV(K)*WEIGHT
	   END DO
	  END DO
	 END DO
	END DO
	ELORC(S)=ENOLD-SETRC(S)

        RETURN
        END

!************************************************************************
!			     DRIFTR
!       	Calculate changes due to radial drift
!************************************************************************
	SUBROUTINE DRIFTR

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
	save
      	integer :: UR, i, j, j0, j1, k, l, n
	real(kind=Real8_):: p4, qs, x, fup, r, corr, cgr1, cgr2, cgr3
	real(kind=Real8_):: CGR(NR,NT,NE,NPA),RGR(NR,NT,NE,NPA), DtsNextLocal
	real(kind=Real8_):: F(NR+2),FBND(NR),C(NR,NT),LIMITER,CR(NR,NT)

	DtsNextLocal = 100000.0
	QS=1.
	IF (S.EQ.1) QS=-1.
	DO I=1,NR
	DO J=1,NT
	  J0=J-1
	  IF (J.EQ.1) J0=NT-1
	  J1=J+1
	  IF (J.EQ.NT) J1=2
	  CR(I,J)=VR(I)*(VT(I,J0)+VT(I+1,J0)-VT(I,J1)-VT(I+1,J1))/ &
     	(BNES(I,J)+BNES(I+1,J))+ &
        (EIP(I,J)+EIP(I+1,J))/(BNES(I,J)+BNES(I+1,J))*DT/DR
	ENDDO
	ENDDO

	DO 1 K=2,NE
	  P4=DT*EKEV(K)*1000.0*(GREL(K)+1)/GREL(K)/DPHI/DR/QS
	DO 1 L=2,NPA
	DO 1 J=1,NT			!	< FBND(I) >
	  J0=J-1			!  |____.____|____.____|____.____|
	  IF (J.EQ.1) J0=NT-1		!  <   F(I)  >
	  J1=J+1			!      F - average in cell(i,j,k,l)
	  IF (J.EQ.NT) J1=2
	 DO I=1,NR
	  F(I)=F2(S,I,J,K,L)
	  CGR1=FNIS(I+1,J1,L)+FNIS(I,J1,L)-FNIS(I+1,J0,L)-FNIS(I,J0,L)
	  CGR2=BNES(I+1,J1)+BNES(I,J1)-BNES(I+1,J0)-BNES(I,J0)
	  CGR3=CGR1+ &
        (FNIS(I+1,J,L)+FNIS(I,J,L)-2*FNHS(I+1,J,L)-2*FNHS(I,J,L))* &
        CGR2/2./(BNES(I+1,J)+BNES(I,J))
	  CGR(I,J,K,L)=CGR3/(FNHS(I,J,L)+FNHS(I+1,J,L))*P4/2./ &
        (BNES(I,J)+BNES(I+1,J))/(Z(I)+0.5*DR)

	  C(I,J)=CR(I,J)+CGR(I,J,K,L)
	  RGR(I,J,K,L)=C(I,J)*DR/DT
	  DTsNextLocal = min( DTsNextLocal, FracCFL*DT/abs(C(I,J)))

	  ISGN(I,J)=1
	  IF (C(I,J).NE.ABS(C(I,J))) ISGN(I,J)=-1
	 END DO

	  IF (ISGN(NR,J).EQ.1) THEN
	   FBND(1)=0.
	   FBND(NR)=F(NR)
	   UR=NR-1
	  ELSE	   
	   FBND(1)=F(2)
	   UR=NR
	   F(NR+1)=FGEOS(S,J,K,L)*CONF1*FNHS(NR,J,L)
	   F(NR+2)=FGEOS(S,J,K,L)*CONF2*FNHS(NR,J,L)
	  END IF

	  DO I=2,UR
             X=F(I+1)-F(I)
	     FUP=0.5*(F(I)+F(I+1)-ISGN(I,J)*X)
	     IF (ABS(X).LE.1.E-27) FBND(I)=FUP
	     IF (ABS(X).GT.1.E-27) THEN
              N=I+1-ISGN(I,J)
	      R=(F(N)-F(N-1))/X
	      IF (R.LE.0) FBND(I)=FUP
	      IF (R.GT.0) THEN
	       LIMITER=MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
	       CORR=-0.5*(C(I,J)-ISGN(I,J))*X      
	       FBND(I)=FUP+LIMITER*CORR
	      END IF
	     END IF
	  END DO

!........update the solution for next time step
	  DO I=2,NR
     	    F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(I,J)*FBND(I)+C(I-1,J)*FBND(I-1)
	    if(f2(s,i,j,k,l).lt.0) then
	     f2(S,i,j,k,l)=1E-15
            endif
	  END DO
1	CONTINUE
	DtsNext = min(DtsNext, DtsNextLocal)

	RETURN
	END

!************************************************************************
!			     DRIFTP
!       	Calculate changes due to azimuthal drift
!************************************************************************
	SUBROUTINE DRIFTP

	use ModRamMain
	use ModIoUnit, ONLY: UNITTMP_
	implicit none
	save
	integer :: i, iSign, j, j1, k, l, n
	real(kind=Real8_) :: x, fup, r, corr, ome, DtsNextLocal
	real(kind=Real8_):: C(NR,NT,NE,NPA),VPA(NR,NT,NE,NPA), &
      AGR(NR,NT,NE,NPA),GPA1,GPA2,GPA
	real(kind=Real8_):: FBND(NT),F(NT),LIMITER

	DtsNextLocal = 100000.0
	OME=7.3E-5			! Earth's angular velocity [rad/s]
	DO 1 L=2,NPA
	DO 1 K=2,NE
  	DO 1 I=2,NR
	   DO J=1,NT
	     F(J)=F2(S,I,J,K,L)
	   END DO

	   DO J=2,NT
	      J1=J+1
	      IF (J.EQ.NT) J1=2
	    GPA1=FNIS(I,J,L)+FNIS(I,J1,L)+(FNIS(I+1,J1,L)+FNIS(I+1,J,L)- &
          FNIS(I-1,J,L)-FNIS(I-1,J1,L))*Z(I)/2./DR
          GPA2=Z(I)/4./DR* &
               (FNIS(I,J,L)+FNIS(I,J1,L)-2*FNHS(I,J,L)-2*FNHS(I,J1,L))* &
               (BNES(I+1,J1)+BNES(I+1,J)-BNES(I-1,J)-BNES(I-1,J1))/ &
               (BNES(I,J)+BNES(I,J1))
	    GPA=GPA1+GPA2
	  VPA(I,J,K,L)=P2(I,K)*GPA/(FNHS(I,J,L)+FNHS(I,J1,L))/ &
        (BNES(I,J)+BNES(I,J1))*DPHI/DT*Z(I)

            C(I,J,K,L)=((VT(I+1,J)+VT(I+1,J1)-VT(I-1,J)-VT(I-1,J1))*P1(I)- &
                 P2(I,K)*GPA/(FNHS(I,J,L)+FNHS(I,J1,L)) - &
                 (EIR(I,J1)+EIR(I,J))/Z(I)*DT/DPHI)/ &
                 (BNES(I,J)+BNES(I,J1))+OME*DT/DPHI

	      AGR(I,J,K,L)=C(I,J,K,L)*DPHI/DT
	      DTsNextLocal = min(DTsNextLocal, FracCFL*DT/abs(C(I,J,K,L)))
	      ISIGN=1
	      IF(C(I,J,K,L).NE.ABS(C(I,J,K,L))) ISIGN=-1
	      X=F(J1)-F(J)
	      FUP=0.5*(F(J)+F(J1)-ISIGN*X)
	      IF (ABS(X).LE.1.E-27) FBND(J)=FUP
	      IF (ABS(X).GT.1.E-27) THEN
	        N=J+1-ISIGN
		IF (N.GT.NT) N=N-NT+1 
 	        R=(F(N)-F(N-1))/X
	        IF (R.LE.0) FBND(J)=FUP 
	        IF (R.GT.0) THEN
	          LIMITER=MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
		  CORR=-0.5*(C(I,J,K,L)-ISIGN)*X      
	          FBND(J)=FUP+LIMITER*CORR
	        END IF
	      END IF
	  END DO
	  C(I,1,K,L)=C(I,NT,K,L)
	  FBND(1)=FBND(NT)

	  DO J=2,NT
	    F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(I,J,K,L)*FBND(J)+ &
          C(I,J-1,K,L)*FBND(J-1)
	    if(f2(s,i,j,k,l).lt.0) then
	     f2(S,i,j,k,l)=1E-15
            endif
	  END DO
	  F2(S,I,1,K,L)=F2(S,I,NT,K,L)
1	CONTINUE
	DtsNext = min(DtsNext, DtsNextLocal)

	RETURN
	END

!**************************************************************************
!			DRIFTE
!		Calculate energization along the drift path 
!**************************************************************************
	SUBROUTINE DRIFTE

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
	save
	integer :: i, isign, j, j0, j2, k, l, n
	real(kind=Real8_)::ezero,gpa,gpr1,gpr2,gpr3,gpp1,gpp2,edt1,qs, &
      drdt, dpdt, dbdt1, didt1, x, fup, r, corr, ome, DtsNextLocal
	real(kind=Real8_)::GRZERO, DRD1(NR,NT),DPD1(NR,NT),GPR(NR,NT,NPA), &
      GPP(NR,NT,NPA),DRD2(NR,NT,NPA),DPD2(NR,NT,NPA), &
      EGR(NR,NT,NE,NPA)
	real(kind=Real8_)::FBND(NE),F(0:NE+2),C(NE),LIMITER

	DtsNextLocal=10000.0
	QS=1.
	IF (S.EQ.1) QS=-1.
	OME=7.3E-5
	EZERO=EKEV(1)-WE(1)
	GRZERO=1.+EZERO*1000.*Q/RMAS(S)/CS/CS
	F(NE+1)=0.
	F(NE+2)=0.
	DO 1 J=1,NT
	  J0=J-1
	  IF (J.EQ.1) J0=NT-1
	  J2=J+1
	  IF (J.EQ.NT) J2=2
	DO 1 I=2,NR
	   DRD1(I,J)=(EIP(I,J)*Z(I)-(VT(I,J2)-VT(I,J0))/2./DPHI)/BNES(I,J)
	   DPD1(I,J)=OME*Z(I)+((VT(I+1,J)-VT(I-1,J))/2/DR-EIR(I,J))/BNES(I,J)

	DO 1 L=2,NPA
	   GPA=(1.-FNIS(I,J,L)/2./FNHS(I,J,L))/BNES(I,J)
	   GPR1=GPA*(BNES(I+1,J)-BNES(I-1,J))/2./DR
           GPR2=-FNIS(I,J,L)/FNHS(I,J,L)/Z(I)
           GPR3=-(FNIS(I+1,J,L)-FNIS(I-1,J,L))/2./DR/FNHS(I,J,L)
	   GPR(I,J,L)=GPR1+GPR2+GPR3
	   GPP1=GPA*(BNES(I,J2)-BNES(I,J0))/2./DPHI
	   GPP2=-(FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI/FNHS(I,J,L)
	   GPP(I,J,L)=GPP1+GPP2
	   DRD2(I,J,L)=(FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI+ &
         (FNIS(I,J,L)-2*FNHS(I,J,L))*(BNES(I,J2)-BNES(I,J0))/ &
         4/BNES(I,J)/DPHI
	   DPD2(I,J,L)=FNIS(I,J,L)+(FNIS(I+1,J,L)-FNIS(I-1,J,L))* &
         Z(I)/2/DR+Z(I)*(FNIS(I,J,L)-2*FNHS(I,J,L))/4/DR* &
         (BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)

	 DO K=2,NE
	  F(K)=F2(S,I,J,K,L)
	 END DO

	 F(1)=F(2)*GREL(1)/GREL(2)*SQRT((GREL(2)**2-1)/(GREL(1)**2-1))
	 F(0)=F(1)*GRZERO/GREL(1)*SQRT((GREL(1)**2-1)/(GRZERO**2-1))

	 DO K=1,NE

	  EDT1=EBND(K)*1e3*(GRBND(K)+1)/2/GRBND(K)/FNHS(I,J,L)/ &
        Z(I)/BNES(I,J)/QS
	  DRDT=DRD1(I,J)+EDT1*DRD2(I,J,L)*Z(I)
	  DPDT=DPD1(I,J)-EDT1*DPD2(I,J,L)
	  dBdt1=dBdt(I,J)*(1.-FNIS(I,J,L)/2./FNHS(I,J,L))*Z(I)/BNES(I,J)
	  dIdt1=-dIdt(I,J,L)*Z(I)/FNHS(I,J,L)
	  C(K)=EDOT(I,K)*(GPR(I,J,L)*DRDT+GPP(I,J,L)*DPDT + &
        dBdt1 + dIdt1)

	  EGR(I,J,K,L)=C(K)/DT
	  DTsNextLocal = min(DTsNextLocal, FracCFL*DT*DE(K)/abs(C(K)))
	  ISIGN=1
	  IF(C(K).NE.ABS(C(K))) ISIGN=-1
	  X=F(K+1)-F(K)
	  FUP=0.5*(F(K)+F(K+1)-ISIGN*X)
	  IF (ABS(X).LE.1.E-27) FBND(K)=FUP
	  IF (ABS(X).GT.1.E-27) THEN
	    N=K+1-ISIGN
 	    R=(F(N)-F(N-1))/X
	    IF (R.LE.0) FBND(K)=FUP
	    IF (R.GT.0) THEN
	      LIMITER=MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
	      CORR=-0.5*(C(K)/DE(K)-ISIGN)*X   
	      FBND(K)=FUP+LIMITER*CORR
	    END IF
	  END IF  
	 END DO

	 DO K=2,NE
          F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(K)/WE(K)*FBND(K)+C(K-1)/WE(K)* &
               FBND(K-1)
	    if(f2(s,i,j,k,l).lt.0) then
	     f2(S,i,j,k,l)=1E-15
            endif
	 END DO
1	CONTINUE
	DtsNext = min(DtsNext, DtsNextLocal)

	RETURN
	END

!**************************************************************************
!			DRIFTMU
!	Calculate pitch angle changes along the drift path
!**************************************************************************
	SUBROUTINE DRIFTMU

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
	save
	integer :: i, j, j0, j1, k, l, n
	real(kind=Real8_)::gmr1, gmr2, gmr3, gmp1, gmp2, qs, &
      drdm, dpdm, dbdt2, dibndt2, x, fup, r, corr, ome, DtsNextLocal
	real(kind=Real8_)::CMUDOT(NR,NT,NPA),EDT(NPA),GMR(NR,NT,NPA), &
      GMP(NR,NT,NPA),DRM2(NR,NT,NPA),DPM2(NR,NT,NPA),UGR(NR,NT,NE,NPA)
	real(kind=Real8_)::FBND(NPA),F(NPA+1),C(NPA),LIMITER,DRM1(NR,NT), &
      DPM1(NR,NT)
	integer :: URP(NPA),ISGM(NPA)

	DtsNextLocal = 10000.0
	QS=1.
	IF (S.EQ.1) QS=-1.
	OME=7.3E-5
	DO 1 K=2,NE	
	DO 1 J=1,NT
	  J0=J-1
	  IF (J.EQ.1) J0=NT-1
	  J1=J+1
	  IF (J.EQ.NT) J1=2
	DO 1 I=2,NR
	  DRM1(I,J)=(EIP(I,J)*Z(I)-(VT(I,J1)-VT(I,J0))/2/DPHI)/BNES(I,J)
	  DPM1(I,J)=OME*Z(I)+((VT(I+1,J)-VT(I-1,J))/2/DR-EIR(I,J))/BNES(I,J)

	  C(1)=0.
	  FBND(1)=0.
	  UGR(I,J,K,1)=0.
	  C(NPA)=0.
	  FBND(NPA)=0.

	 DO L=2,NPA
	  F(L)=F2(S,I,J,K,L)
	  CMUDOT(I,J,L)=MUDOT(I,L)*BOUNIS(I,J,L)/BOUNHS(I,J,L)
	  GMR1=(BNES(I+1,J)-BNES(I-1,J))/4/DR/BNES(I,J)
	  GMR2=1/Z(I)
	  GMR3=(BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))/2/DR/BOUNIS(I,J,L)
	  GMR(I,J,L)=GMR1+GMR2+GMR3
	  GMP1=(BNES(I,J1)-BNES(I,J0))/4/DPHI/BNES(I,J)
	  GMP2=(BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI/BOUNIS(I,J,L)
	  GMP(I,J,L)=GMP1+GMP2
	  EDT(L)=EKEV(K)*1e3*(GREL(K)+1)/2/GREL(K)/BOUNHS(I,J,L)/ &
        Z(I)/BNES(I,J)/QS
	  DRM2(I,J,L)=(BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI+ &
        (BOUNIS(I,J,L)-2*BOUNHS(I,J,L))*(BNES(I,J1)-BNES(I,J0))/ &
        4/BNES(I,J)/DPHI
	  DPM2(I,J,L)=BOUNIS(I,J,L)+(BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))* &
        Z(I)/2/DR+(BOUNIS(I,J,L)-2*BOUNHS(I,J,L))*Z(I)/4/DR* &
        (BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)

	  DRDM=DRM1(I,J)+EDT(L)*DRM2(I,J,L)*Z(I)
	  DPDM=DPM1(I,J)-EDT(L)*DPM2(I,J,L)
	  dBdt2=dBdt(I,J)/2./BNES(I,J)*Z(I)
	  dIbndt2=dIbndt(I,J,L)*Z(I)/BOUNIS(I,J,L)
	  C(L)=-CMUDOT(I,J,L)*(GMR(I,J,L)*DRDM+GMP(I,J,L)*DPDM + &
        dBdt2 +dIbndt2)

	  UGR(I,J,K,L)=C(L)/DT
          DTsNextLocal=min(DTsNextLocal,FracCFL*DT*DMU(L)/max(1e-32,abs(C(L))))
	  ISGM(L)=1
	  IF(C(L).NE.ABS(C(L))) ISGM(L)=-1
	  IF(ISGM(L).EQ.1) THEN
	   URP(L)=NPA-1
	  ELSE
	   URP(L)=NPA-2
	  ENDIF
	 END DO
	 F(1)=F(2)
	 FBND(NPA-1)=F(NPA)

	 DO L=2,NPA-2
	  X=F(L+1)-F(L)
	  FUP=0.5*(F(L)+F(L+1)-ISGM(L)*X)
	  IF (ABS(X).LE.1.E-27) FBND(L)=FUP
	  IF (ABS(X).GT.1.E-27) THEN
	    N=L+1-ISGM(L)
 	    R=(F(N)-F(N-1))/X
	    IF (R.LE.0) FBND(L)=FUP
	    IF (R.GT.0) THEN
	      LIMITER=MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
 	      CORR=-0.5*(C(L)/DMU(L)-ISGM(L))*X
	      FBND(L)=FUP+LIMITER*CORR
	    END IF
	  END IF  
	 END DO

	 DO L=2,NPA-1
	  F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(L)/WMU(L)*FBND(L)+ &
        C(L-1)/WMU(L)*FBND(L-1)
	    if(f2(s,i,j,k,l).lt.0) then
	     f2(S,i,j,k,l)=1E-15
            endif
	 END DO
	 F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)*FNHS(I,J,NPA)*MU(NPA)/ &
       FNHS(I,J,NPA-1)/MU(NPA-1)
1	CONTINUE
	DtsNext = min(DtsNext, DtsNextLocal)

	RETURN
	END

!************************************************************************
!			CHAREXCHANGE 
!       	Calculate the decay due to charge exchange
!************************************************************************
	SUBROUTINE CHAREXCHANGE

	use ModRamMain
	implicit none
 	save
	integer :: i, j, k, l
	
	DO 10 K=2,NE
	DO 10 J=1,NT
	DO 10 I=2,NR
	 DO 1 L=2,NPA
1	  F2(S,I,J,K,L)=F2(S,I,J,K,L)*ACHAR(I,J,K,L)
10	CONTINUE

	RETURN
	END


!************************************************************************
!			ATMOSPHERIC LOSSES
!       	Calculate loss to the atmosphere due to collisions
!************************************************************************
	SUBROUTINE ATMOL

	use ModRamMain
	use ModIoUnit,    ONLY: UNITTMP_
	implicit none
 	save
	integer :: i, j, k, l

	DO 10 K=2,NE
	DO 10 J=1,NT
	DO 10 I=2,NR
	 DO 2 L=UPA(I),NPA
2	  F2(S,I,J,K,L)=F2(S,I,J,K,L)*ATLOS(I,K)**(1/FNHS(I,J,L))
10	CONTINUE

	RETURN
	END

!************************************************************************
!			WAVELO
!       Calculate loss due to waves everywhere using electron lifetimes
!************************************************************************
	SUBROUTINE WAVELO

	use ModRamMain
	implicit none
	save
	integer :: i, j, k, l, j1, i1
	real(kind=Real8_):: TAU_LIF,RLpp(NT),Bw,Kpmax

	Bw=30.
	IF(KP.GE.4.0) Bw=100.                                 ! pT
	Kpmax=KP	 ! need Kpmax from previous 12 hrs, TO FIX!!!
	DO J=1,NT
	  J1=(J-1)*IP1
	  RLpp(J)=5.39-0.382*KPmax	 ! PP from Moldwin et al. [2002]
	DO I=2,NR
	  I1=(I-2)*IR1+3
	  IF (DoUsePlane_SCB.and.NECR(I1,J1).gt.50.) RLpp(J)=LZ(I) 
	ENDDO
	ENDDO

	DO 10 K=2,NE
	DO 10 I=2,NR
	DO 10 J=1,NT
	DO 10 L=2,NPA
	 IF(LZ(I).LE.RLpp(J)) THEN
	   TAU_LIF=WALOS1(I,K)*((10./Bw)**2)
	 ELSEIF (LZ(I).GT.RLpp(J)) THEN
           IF(EKEV(K).LE.1000.) THEN
            TAU_LIF=WALOS2(I,K)*(1+WALOS3(I,K)/WALOS2(I,K))
	    if (ekev(k).le.1.1) then
	     tau_lif=tau_lif*37.5813*exp(-1.81255*ekev(k))
	    else if (ekev(k).gt.1.1.and.ekev(k).le.5.) then
	     tau_lif=tau_lif*(7.5-1.15*ekev(k))
	    else
	     tau_lif=tau_lif
	    endif
	   ELSEIF(EKEV(K).GT.1000.) THEN
	    TAU_LIF=5.*3600*24/KP
	   ENDIF
	  ENDIF
	  F2(S,I,J,K,L)=F2(S,I,J,K,L)*EXP(-DT/TAU_LIF)
10	CONTINUE

	RETURN
	END


!*************************************************************************
!			 	WPADIF
!     Routine calculates the decay of the distribution function
!        due to WPI pitch angle diffusion using implicit scheme
!*************************************************************************
	SUBROUTINE WPADIF

	use ModRamMain
	implicit none
	save
	integer :: i, j, k, l
	real(kind=Real8_):: F(NPA),AN,BN,GN,RP,DENOM,RK(NPA),RL(NPA),FACMU(NPA)

!        write(*,*) "starting wpadiff"

	DO 1 J=1,NT
	DO 1 I=2,NR
	DO 1 K=2,NE
	  DO L=2,NPA
	     FACMU(L)=FNHS(I,J,L)*MU(L)
	     F(L)=F2(S,I,J,K,L)/FACMU(L)
	  END DO
	  FACMU(1)=FNHS(I,J,1)*MU(1)

	  F(1)=F(2)	   			! lower b.c.
	  RK(1)=0.	   			
	  RL(1)=-1.	   			
	  DO L=2,NPA-1
	   AN=(ATAW(I,J,K,L)+ATAC(I,J,K,L))/DMU(L)       ! Hiss & chorus
	   GN=(ATAW(I,J,K,L-1)+ATAC(I,J,K,L-1))/DMU(L-1) !  "
           AN=AN*DT/FACMU(L)/WMU(L)
           GN=GN*DT/FACMU(L)/WMU(L)
	   BN=AN+GN   		                        
	if (abs(-1-bn).lt.(abs(an)+abs(gn))) then
	 open(20,file=trim(PathRamOut)//'diffcf_e.dat',status='unknown', &
          position='append')
	 write(20,*) ' T=',T/3600,' hr'
	 write(20,*) 'i=',i,' j=',j,' k=',k,' l=',l
	 write(20,*) 'an=',AN,' -1-bn=',(-1-BN),' gn=',GN
	 close(20)
	endif
	   RP=F(L)
	   DENOM=BN+GN*RL(L-1)+1
	   RK(L)=(RP+GN*RK(L-1))/DENOM
	   RL(L)=-AN/DENOM
	  END DO	   

	  F2(S,I,J,K,NPA-1)=RK(NPA-1)/(1+RL(NPA-1))	! upper b.c.
	  DO L=NPA-2,1,-1
	   F2(S,I,J,K,L)=RK(L)-RL(L)*F2(S,I,J,K,L+1)
	  END DO
	  F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)
	 DO L=1,NPA
	  F2(S,I,J,K,L)=F2(S,I,J,K,L)*FACMU(L)
	 ENDDO
1	CONTINUE

	RETURN
	END

! *************************************************************************
!				 GEOSB
!   			Boundary conditions set up
!**************************************************************************
	SUBROUTINE GEOSB(nFive)

	use ModRamMain
	use ModRamCouple, ONLY: FluxBats_IIS, generate_flux, TypeMHD, &
      MhdDensPres_VII,FluxBats_anis
	use ModIoUnit,    ONLY: UNITTMP_
	use ModRamGeomlt, ONLY: get_geomlt_flux
	USE ModConst, ONLY: cProtonMass
	use ModRamIO, ONLY: UseNewFmt, RamFileName, gen_old_timetags

	implicit none
	save
	integer, intent(IN) :: nFive
	integer :: ij, ik, j, jw, k, l, nLines, nFiveDay
	real(kind=Real8_) :: bexp, ahe0, ahe1, gexp, doy, azir, &
      fracComposition
	character(len=200) :: NameFluxFile
        CHARACTER*90 HEADER
	real(kind=Real8_)::RELAN(NTL,NEL),FLAN(NT,NEL),FluxLanl(nT,nE)

	! Create Young et al. composition ratios.
	BEXP=(4.5E-2)*EXP(0.17*KP+0.01*F107)
	AHE0=6.8E-3
	AHE1=0.084
	GEXP=0.011*EXP(0.24*KP+0.011*F107)
	GEXP=BEXP*(AHE0/GEXP+AHE1)
	select case(s)
	case(2)! Fraction H+
	   fracComposition=4./(4.+BEXP+2.*GEXP)
	case(3)! Fraction He+
	   fracComposition=2.*GEXP/(4.+BEXP+2.*GEXP)
	case(4)! Fraction O+
	   fracComposition=BEXP/(4.+BEXP+2.*GEXP)
	case default ! Electrons.
	   fracComposition=1.0
	end select

	! Zero out FGEOS at this species.
	FGEOS(S,:,:,:)=0. 

!.......Read LANL flux (1/cm2/s/sr/keV) assumed isotropic
	IF ((boundary .eq. 'LANL') .OR. (boundary .eq. 'PTM')) THEN
	      ! LANL interpolated flux files.
	    if(s==1)then
	       call get_geomlt_flux('elec', FluxLanl)
	    else
	       call get_geomlt_flux('prot', FluxLanl)
	    end if

! ..... fix corrupted bc, NAN check doesn't work? limit ion flux:
            if (S.gt.1) then
	    do ij=2,nT; do ik=1,nE;
	       if (FluxLanl(iJ,iK).gt.1e7) FluxLanl(iJ,iK)=FluxLanl(iJ-1,iK)
               if (FluxLanl(iJ,iK).gt.1e7) then
                 FluxLanl(iJ,iK)=1e7
                 write(*,*) ' in GEOSB: limit ion flux to 1e7'
               endif
	    end do; end do
            do ik=1,nE
               FluxLanl(1,iK)=FluxLanl(nT,iK)
            end do
            end if
! ..... corrupted bc, end of fix

	    ! Adjust flux for composition
	    FluxLanl=FluxLanl*fracComposition
	    do ik=1,NE; do ij=1,nT; do L=2,upa(NR)-1;
	       FGEOS(s,iJ,iK,L)=FluxLanl(iJ,iK) * FFACTOR(NR,IK,L)
	    end do; end do; end do
	 
	   ELSE	IF (boundary .EQ. 'SWMF') THEN
!.......Read SWMF flux (1/cm2/s/sr/keV) assumed isotropic
	   if (IsComponent) then ! If SWMF component, get flux from Bats.
	      write(*,*) 'RAM_SCB: Getting flux from BATS-R-US'
	      do iK=1, nE
		 do iJ=1,nT
		    do L=2,UPA(NR)-1
		       if(.not.DoAnisoPressureGMCoupling)then
			 FGEOS(s,iJ,iK,L)=FluxBats_IIS(iK,iJ,s)*FFACTOR(NR,IK,L)
		       else
!                        anisotropy coupling is pitch angle dependent.
			 FGEOS(s,iJ,iK,L)=FluxBats_anis(iK,L,iJ,s)* &
        			FFACTOR(NR,iK,L)
		       endif
		    end do
		 end do
	      end do

	   else			! Open file if running stand alone.
	      if(DoMultiBcsFile) then
		 write(NameFluxFile, "(a,'/',a,'_',i1,'.swf')") &
        trim(PathSwmfOut), ST7, S
		 nlines = 25
	      else
		 write(NameFluxFile, "(a,'/',a,'.swf')") trim(PathSwmfOut), ST7
		 nlines=49
	      endif	      
	      write(*,'(2a)')'RAM_SCB: Loading flux from ',trim(NameFluxFile)
	      OPEN(UNIT=UNITTMP_,FILE=trim(NameFluxFile),STATUS='OLD')
	      READ(UNITTMP_,10) HEADER
	      DO IJ=1,nLines
	         READ(UNITTMP_,*) DOY,AZIR,(RELAN(IJ,IK),IK=1,NE)
	      ENDDO
	      CLOSE(UNITTMP_)
10	FORMAT(A90)
	      
	      DO IK=1,NE
		 DO IJ=1,NT
		    IF (NT.EQ.49) JW=IJ
		    if ((NT.EQ.25) .and. DoMultiBcsFile) then 
		       JW=IJ 
		    else 
		       JW=2*IJ-1
		    end if
		    FLAN(IJ,IK) = RELAN(JW,IK)
!       same mass density at geo as without composition
		    IF ((S.EQ.2) .and. .not.DoMultiBcsFile) &
           FLAN(IJ,IK) = FLAN(IJ,IK)/(1+16*BEXP+4*GEXP) 
		    IF ((S.EQ.3) .and. .not.DoMultiBcsFile) &
           FLAN(IJ,IK) = FLAN(IJ,IK)*GEXP/(2+32*BEXP+8*GEXP) 
		    IF ((S.EQ.4) .and. .not.DoMultiBcsFile) &
           FLAN(IJ,IK) = FLAN(IJ,IK)*BEXP/(4+64*BEXP+16*GEXP)
		    DO L=2,UPA(NR)-1
		      FGEOS(S,IJ,IK,L)=FLAN(IJ,IK)*FFACTOR(NR,IK,L)
		    ENDDO
		 ENDDO
	      ENDDO
	   end if
	END IF

!.......Write interpolated fluxes to file.
	if(UseNewFmt)then
	   NameFluxFile=trim(PathRamOut)//RamFileName('Dsbnd/ds'//St2, &
         'dat',TimeRamNow)
	else
	   call gen_old_timetags(TimeRamStart, TimeRamNow, nFiveDay=nFiveDay)
	   write(NameFluxFile,'(a,i4.4,a)')trim(PathRamOut)//'Dsbnd/ds_', &
         nFiveDay,ST2//'.dat'
	end if
	OPEN(UNIT=UNITTMP_,FILE=NameFluxFile, STATUS='UNKNOWN')
	WRITE(UNITTMP_,20)'EKEV FGEOSB [1/cm2/s/ster/keV] T=',T/3600,Kp,F107
20	FORMAT(A33,F8.3,5H, Kp=,F6.2,7H, F107=,F8.2)
70	DO K=2,NE
	   WRITE(UNITTMP_,50) EKEV(K),(FGEOS(S,J,K,2)/FFACTOR(NR,K,2),J=1,NT)
	END DO
	CLOSE(UNITTMP_) 
50	FORMAT(F7.2,50(1PE11.3))	

        RETURN
        END

!*************************************************************************
!			 	ANISCH
!     		Calculate pressure in equatorial plane
!*************************************************************************
	SUBROUTINE ANISCH

	use ModRamMain
        use ModRamFunctions
	use ModIoUnit,    ONLY: UNITTMP_

	implicit none
	save
	integer :: i, iwa, j, k, klo, l, iz, ier, kn, i1, j1
	real(kind=Real8_):: cv, rfac, rnhtt, edent, pper, ppar, rnht, Y, &
           eden,sume,suma,sumn,ernm,epma,epme,anis,epar,taudaa,taudea,taudee, &
           gausgam,anist,epart,fnorm,xfrl,Bw,esu,omega,er1,dx
        real(kind=Real8_)::MUBOUN,DWAVE(NPA),KEVERG,CMRA(SLEN),BWAVE(NR,NT), &
           AVDAA(NPA),TAVDAA(NPA),DAA(NE,NPA,Slen),DUMP(ENG,NCF), &
           DN(2,25),EPO(2,25),AII(2,25),XFR(NR,NT),XFRe(NCF),ALENOR(ENG), &
           DUME(ENG,NCF),DVV(NE,NPA,Slen),AVDVV(NPA),TAVDVV(NPA),WCDT(NR,NT), &
           XFRT(NR,NT),PA(NPA),DAMR(NPA,NCO),DAMR1(NPA),GREL_new(Ny),BOUNHS_(Nx)
	INTEGER::MINP(NT),MAXP(NT),KHI(5)
	CHARACTER*80 HEADER
	DATA khi/6, 10, 25, 30, 35/ ! ELB=0.1 keV -> 0.4,1,39,129,325 keV 

	IF(DoUseWPI.and.MOD(INT(T),3600).EQ.0.) THEN
	 OPEN(UNIT=UNITTMP_,FILE=trim(PathRamOut)//ST0//ST1//ST2//'.wans', &
	       STATUS='UNKNOWN')
	 WRITE(UNITTMP_,56) T/3600,KP 				! total
56	 FORMAT(2HT=,F8.3,2X,3HKp=,F3.1,/,2X, &
     	 'L   PHI     ANIS   EDEN[keV/cm3]   RNHT[1/cm3]  PPER[keV/cm3]', &
     	 '  PPAR[keV/cm3]')
	ENDIF

	cv=CS*100		! speed of light in [cm/s]
	esu=Q*3E9		! elementary charge in [esu]
	RFAC=4*pi/cv
	gausgam=1.E-5
	khi(5)=NE
!.......calculate ring current parameters
	DO 15 I=2,NR
          I1=(I-2)*IR1+3
	DO 15 J=1,NT
          J1=(J-1)*IP1
          IF (DoUsePlane_SCB) THEN
             XNE(I,J)=NECR(I1,J1)
          ELSE
             if (S.EQ.1.and.I.eq.2.and.J.eq.1) &
               write (*,*) " Need to specify Ne if using WPI"
          ENDIF
	  klo=2
	  PPERT(S,I,J)=0.
	  PPART(S,I,J)=0.
	  RNHTT=0.
	  EDENT=0.
	  do iwa=1,5
	    PPER=0.
	    PPAR=0.
	    RNHT=0.
	    EDEN=0.
	    DO K=klo,khi(iwa)
	      F2(S,I,J,K,1)=F2(S,I,J,K,2)
	      SUME=0.
	      SUMA=0.
	      SUMN=0.
	      DO L=1,UPA(I)-1
	        ERNM=WMU(L)/FFACTOR(I,K,L)/FNHS(I,J,L)
	        EPMA=ERNM*MU(L)*MU(L)
	        EPME=ERNM-EPMA
	        SUME=SUME+F2(S,I,J,K,L)*EPME
	        SUMA=SUMA+F2(S,I,J,K,L)*EPMA
	        SUMN=SUMN+F2(S,I,J,K,L)*ERNM
	      ENDDO
	      PPER=PPER+EPP(K)*SUME
	      PPAR=PPAR+EPP(K)*SUMA
	      RNHT=RNHT+ERNH(K)*SUMN
	      EDEN=EDEN+ERNH(K)*EKEV(K)*SUMN
	    ENDDO
	    ANIS=PPER/2./PPAR-1.
	    EPAR=2*PPAR/RNHT
	    RNHT=RNHT*RFAC
	    EDEN=EDEN*RFAC
	    PPAR=2*RFAC*PPAR
	    PPER=RFAC*PPER
!.......write anisotropy, kT parallel, RC density
	    klo = khi(iwa)+1
	    EDENT=EDENT+EDEN
	    PPERT(S,I,J)=PPERT(S,I,J)+PPER
	    PPART(S,I,J)=PPART(S,I,J)+PPAR
	    RNHTT=RNHTT+RNHT
	  enddo
	  ANIST=PPERT(S,I,J)/PPART(S,I,J)-1.
	  EPART=PPART(S,I,J)/RNHTT
	  IF(DoUseWPI.and.MOD(INT(T),3600).EQ.0.) WRITE(UNITTMP_,551)  &
           LZ(I),PHI(J),ANIST,EDENT,RNHTT,PPERT(S,I,J),PPART(S,I,J)
15	CONTINUE
	CLOSE(UNITTMP_)
551	FORMAT(F5.2,F10.6,10(1PE12.3))

	IF (S.GT.1) GO TO 100			! if ions
	IF(MOD(INT(T),INT(Dt_bc)).EQ.0.and.DoUseWPI) THEN
!.......zero PA diffusion coefficients
30	DO 17 I=1,NR
	DO 17 J=1,NT
	DO 17 K=1,NE
	DO 17 L=1,NPA
	 ATAW(I,J,K,L)=0.0			! hiss
	 ATAW_emic(I,J,K,L)=0.0			! EMIC
	 ATAC(I,J,K,L)=0.0			! chorus
17	CONTINUE

!.......linear interpolation for CHORUS b-aver diff coefficients
	OPEN(UNIT=UNITTMP_,FILE=trim(PathRamOut)//ST0//ST1//ST2//'.wdc', &
          STATUS='UNKNOWN')
	WRITE(UNITTMP_,7) T/3600,KP
7	FORMAT(2X,3HT =,F8.3,2X,4HKp =,F6.2,' CHORUS diff coeff')

	DO L=1,NPA
	  PA(L)=ACOSD(MU(L))
	  DO IZ=1,NCO
	     DAMR(L,IZ)=0.
	  ENDDO
	  DAMR1(L)=0.
	ENDDO
	DO 41 I=2,NR
	 DO 41 J=1,NT
	  IF (XNE(I,J).LE.50.) THEN			! outside pp
	    xfrl=CS*SQRT(XNE(I,J)*RMAS(S)*40*PI)/10./BNES(I,J) ! Fpe/Fcyc
	    fnorm=1			! b-av, no Bw-dep
	    write(UNITTMP_,550)
	    write(UNITTMP_,553) LZ(I),MLT(J),fnorm,XNE(I,J),' Fp/Fc=',xfrl
! interpolate PA diffusion coefficients for implicit scheme
            DO K=2,NE
              DO L=1,NPA
                IF (DoUseBASdiff) THEN
                  DAMR1(L)=LOG10(CDAAR(I,J,K,L))
                ELSE
                  DAMR1(L)=LOG10(BDAAR(I,J,K,L))
                ENDIF
              ENDDO
              DO L=1,NPA
                MUBOUN=MU(L)+0.5*WMU(L)
                CALL LINTP(PA,DAMR1,NPA,PAbn(L),Y,IER)
                taudaa=10.**Y*fnorm			! <Daa/p2> [1/s]
                if (taudaa.gt.1e0) then
                   print*,'taudaa=',taudaa,' L=',LZ(I),' MLT=',MLT(J)
                   taudaa=1e-1
                endif
                if (taudaa.lt.1e-30) taudaa=1e-30
                DWAVE(L)=taudaa* &
                     (1.-MUBOUN*MUBOUN)*MUBOUN*BOUNHS(I,J,L)
                ATAC(I,J,K,L)=DWAVE(L)      ! call WPADIF twice, implicit
                if (l.eq.2.or.l.eq.6.or.l.eq.10.or.l.eq.15.or.l.eq.20.or. &
                     l.eq.25.or.l.eq.30.or.l.eq.40.or.l.eq.50) then
                     write(UNITTMP_,552) ekev(k),PAbn(L),taudaa,1/taudaa
                endif
              END DO
            END DO
         ENDIF
41	CONTINUE
	CLOSE(UNITTMP_)

!	go to 100                           ! no hiss
!.......Wave amplitude for plasmaspheric hiss
	 Bw=30.
	 IF(KP.GE.4.0) Bw=100.                                 ! pT
	 OPEN(UNIT=UNITTMP_,FILE=trim(PathRamOut)//ST0//ST1//ST2//'.wdh', &
          STATUS='UNKNOWN')
	 WRITE(UNITTMP_,8) T/3600,KP
8	 FORMAT(2X,3HT =,F8.3,2X,4HKp =,F6.2,' HISS diff coeff')

!.......linear interpolation for HISS normalized diff coefficients
	DO KN=1,ENG 
	  ALENOR(KN)=LOG10(ENOR(KN))
	  DO IZ=1,NCF
	     DUMP(KN,IZ)=0.
	  ENDDO
	ENDDO
	DO 401 I=2,NR
	 DO 401 J=1,NT
	  IF (XNE(I,J).GT.50.) THEN			! inside pp
	    omega=esu*10*BNES(I,J)/(RMAS(S)*cv)		! Fcyc, cgs
	    xfrl=CS*SQRT(XNE(I,J)*RMAS(S)*40*PI)/10./BNES(I,J) ! Fpe/Fcyc
	    if (xfrl.gt.18) xfrl=18.		! upper interp bound
	    if (xfrl.lt.2) xfrl=2.		! lower interp bound
	    fnorm=omega*(Bw*1e-3)**2*gausgam**2/1e8/BNES(I,J)/BNES(I,J)
	    write(UNITTMP_,550)
	    write(UNITTMP_,553) LZ(I),MLT(J),fnorm,XNE(I,J),' omega=',omega, &
     	' Bhiss=',Bw
	    DO 402 L=1,NPA
	     MUBOUN=MU(L)+0.5*WMU(L)
	     DO IZ=1,NCF
	      DO KN=1,ENG 
		DUMP(KN,IZ)=LOG10(NDAAJ(I,KN,L,IZ))
	      ENDDO
	     ENDDO
	     DO 402 K=2,NE
	       ER1=LOG10(EKEV(K))
	       CALL LINTP2(ALENOR,fpofc,DUMP,ENG,NCF, &
     		ER1,xfrl,Y,IER)
	        DWAVE(L)=10.**Y*fnorm/GREL(K)**2* &
                (1.-MUBOUN*MUBOUN)/MUBOUN		! denorm pa [1/s]
	        ATAW(I,J,K,L)=DWAVE(L)           ! call WPADIF twice, implicit
	     taudaa=dwave(l)/MUBOUN/BOUNHS(I,J,L) 	! <Dmu>
	 if (l.eq.2.or.l.eq.6.or.l.eq.10.or.l.eq.15.or.l.eq.20.or. &
           l.eq.25.or.l.eq.30.or.l.eq.40.or.l.eq.50) then
             write(UNITTMP_,552) ekev(k),PAbn(L),DWAVE(L),1/DWAVE(L)
	 endif
402	 CONTINUE
	ENDIF
401	CONTINUE
	CLOSE(UNITTMP_)
	ENDIF				! end diff coeff loop 

550	format(1X,36hbounce-aver diff coeff & time scales,/, &
     	 1X,4Hekev,4X,4Heqpa,4X,9Hdaa [1/s],4X,8Htdaa [s])
553	format(3H L=,F5.2,5H MLT=,F5.2,7H fnorm=,E10.3, &
     	 4H Ne=,F9.2,A7,E10.3,A7,E10.3)
552	format(2F9.3,6(1PE12.3))

100	RETURN
	END

