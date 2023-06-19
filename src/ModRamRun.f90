!============================================================================
!    Copyright (c) 2016, Los Alamos National Laboratory
!    All rights reserved.
!============================================================================

MODULE ModRamRun
! Contains the subroutine for running RAM and calculating pressure
! Calls subroutines from ModRamDrift, ModRamWPI, and ModRamLoss
! Calculates new F2 and parallel/perpendicular pressures

  implicit none

  contains

!==============================================================================
  SUBROUTINE ram_run

    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: electric, DoUseWPI, DoUseBASdiff, DoUseKpDiff, &
                               DoUsePlasmasphere, DoUseCoulomb, verbose, &
                               DoUseEMIC, DoUseFLC
    use ModRamConst,     ONLY: RE
    use ModRamGrids,     ONLY: nS, NR, NE, NT, NPA
    use ModRamTiming,    ONLY: UTs, T, DtEfi, DtsMin, DtsNext, TimeRamNow, Dts, TimeRamElapsed
    use ModRamVariables, ONLY: Kp, VT, VTOL, VTN, TOLV, LZ, PHI, PHIOFS, FFACTOR, &
                               FLUX, FNHS, DtDriftR, DtDriftP, DtDriftE, DtDriftMu, &
                               NECR, F107, LSDR, ELORC, LSCHA, LSWAE, LSATM, &
                               F2, outsideMGNP, Upa, wMu, Mu, species, LSCOE, LSCSC
    !!!! Module Subroutines/Functions
    use ModRamPlasmasphere, ONLY: plasmasphere
    use ModRamDrift, ONLY: DRIFTPARA, DRIFTR, DRIFTP, DRIFTE, DRIFTMU, DRIFTEND
    use ModRamLoss,  ONLY: CEPARA, CHAREXCHANGE, ATMOL, FLC_Radius, PARA_FLC, FLCScatter
    use ModRamWPI,   ONLY: WAPARA_KP, WPADIF, WAVELO
    !!!! Share Modules
    use ModTimeConvert, ONLY: TimeType
    use ModRamCoul,  ONLY: COULPARA, COULEN, COULMU

    implicit none

    real(DP) :: AVS
    integer  :: i, j, k, l
    integer  :: IS ! Species to advect.

    AVS=7.05E-6/(1.-0.159*KP+0.0093*KP**2)**3/RE  ! Voll-Stern parameter
    DO I=1,NR+1
       DO J=1,NT
          if (electric .ne. 'VOLS') then
             VT(I,J)=VTOL(I,J)+(VTN(I,J)-VTOL(I,J))*(T-TOLV)/DtEfi
          else
             VT(I,J)=AVS*(LZ(I)*RE)**2*SIN(PHI(J)-PHIOFS) ! [V]
          endif
       ENDDO
    ENDDO

    T = UTs

    ! Added capability to couple to plasmaspheric density model
    IF (DoUsePlasmasphere) call plasmasphere(2._dp*dts)

    ! Get field line curvature radius
    IF(DoUseFLC) call FLC_Radius
    
  !$OMP PARALLEL DO
    do iS=1,nS
       !   calls for given species S:
       CALL CEPARA(iS)
       CALL DRIFTPARA(iS)
       CALL COULPARA(iS)

       ! Call routines to calculate the changes of distribution function
       ! considering drifts, charge exchange and atmospheric losses
       CALL DRIFTR(iS)
       CALL DRIFTP(iS)
       CALL DRIFTE(iS)
       CALL DRIFTMU(iS)
       CALL SUMRC(iS)
       LSDR(iS)=LSDR(iS)+ELORC(iS)

       ! energy loss via Coulomb collisions
       if (DoUseCoulomb) then
         CALL COULEN(iS)
         CALL SUMRC(iS)
         LSCOE(iS)=LSCOE(iS)+ELORC(iS)
         ! Coulomb collisions scattering -- unstable?
         CALL COULMU(iS)
         CALL SUMRC(iS)
         LSCSC(iS)=LSCSC(iS)+ELORC(iS)
       endif

       if (species(iS)%WPI) then
          IF (DoUseWPI) THEN
             CALL WPADIF(iS) ! pitch angle diffusion
          ELSE
             CALL WAVELO(iS) ! electron lifetimes
          ENDIF
          CALL SUMRC(iS)
          LSWAE(iS)=LSWAE(iS)+ELORC(iS)
       endif

       ! loss via EMIC wave scattering
       if (species(iS)%EMIC .and. DoUseEMIC) then
          call WPADIF(iS)
       end if

       if (species(iS)%CEX) then
          CALL CHAREXCHANGE(iS)
          CALL SUMRC(iS)
          LSCHA(iS)=LSCHA(iS)+ELORC(iS)
       endif
       
       ! ion loss via FLC scattering
       if (species(iS)%FLC .and. DoUseFLC) then
             call PARA_FLC(iS)
             call FLCScatter(iS)
       end if

       ! loss via atmosphere loss cone
       CALL ATMOL(iS)
       CALL SUMRC(iS)
       LSATM(iS)=LSATM(iS)+ELORC(iS)

       ! time splitting, reverse order
       CALL ATMOL(iS)
       CALL SUMRC(iS)
       LSATM(iS)=LSATM(iS)+ELORC(iS)
       
       ! ion loss via FLC scattering
       if (species(iS)%FLC .and. DoUseFLC) then
          call PARA_FLC(iS)
          call FLCScatter(iS)
       end if
       
       if (species(iS)%CEX) then
          CALL CHAREXCHANGE(iS)
          CALL SUMRC(iS)
          LSCHA(iS)=LSCHA(iS)+ELORC(iS)
       endif

       ! ion loss via EMIC wave scattering
       if (species(iS)%EMIC .and. DoUseEMIC) then
          call WPADIF(iS)
       end if

       if (species(iS)%WPI) then
          IF (DoUseWPI) THEN
             CALL WPADIF(iS) ! pitch angle diffusion
          ELSE
             CALL WAVELO(iS) ! electron lifetimes
          ENDIF
          CALL SUMRC(iS)
          LSWAE(iS)=LSWAE(iS)+ELORC(iS)
       endif

       if (DoUseCoulomb) then
         CALL COULMU(iS)
         CALL SUMRC(iS)
         LSCSC(iS)=LSCSC(iS)+ELORC(iS)
         CALL COULEN(iS)
         CALL SUMRC(iS)
         LSCOE(iS)=LSCOE(iS)+ELORC(iS)
       endif

       CALL DRIFTMU(iS)
       CALL DRIFTE(iS)
       CALL DRIFTP(iS)
       CALL DRIFTR(iS)
       CALL SUMRC(iS)
       LSDR(iS)=LSDR(iS)+ELORC(iS)

       ! Interpolate diffusivity as funtion of Kp, if necessary
       if (DoUseWPI.and.DoUseBASdiff.and.DoUseKpDiff) then
          call WAPARA_Kp(iS)
       end if

       ! Routines needed to clear allocated variables for use with OpenMP
       CALL DRIFTEND
    end do
  !$OMP END PARALLEL DO
    F2(:,:,nT,:,:) = F2(:,:,1,:,:)

    if (verbose) then
       write(*,'(1x,a,2E10.2)') 'Max, Min FluxQ = ', &
             maxval(F2(:,:,:,:,:)), minval(F2(:,:,:,:,:), MASK=F2(:,:,:,:,:).GT.0)
    endif

    !!!! For now set the flux to 0 outside magnetopause.
    ! Setting to 0 was causing issues, instead keep the flux but track the
    ! magnetopause in all output routines. This is a bad fix, but works for now.
    ! We will want to change this later. -ME
    DO I = 1, nR
       DO J = 1, nT
          if (outsideMGNP(I,J) == 1) F2(:,I,J,:,:) = 1.d-31
       ENDDO
    ENDDO
    !!!!

    DtsNext = min(minval(DtDriftR), minval(DtDriftP), minval(DtDriftE), minval(DtDriftMu))
    DtsNext = max(DtsNext, DtsMin)

    ! Update flux and pressure totals
    DO iS = 1,nS
       CALL ANISCH(iS)
       DO I = 2, NR
          DO K = 2, NE
             DO L = 2, NPA
                DO J = 1, NT-1
                   FLUX(iS,I,J,K,L) = F2(iS,I,J,K,L)/FFACTOR(iS,I,K,L)/FNHS(I,J,L)
                if (abs(flux(iS,i,j,k,l)) > huge(flux(iS,i,j,k,l))) then
                 write(*,*) 'in ModRamRun flux=inf ', iS,i,j,k,l, flux(iS,i,j,k,l)
                endif
               ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE Ram_Run

!************************************************************************
!                       SUMRC     
!               Calculate the total energy & energy loss 
!************************************************************************
  SUBROUTINE SUMRC(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NR, NE, NPA, NT
    use ModRamVariables, ONLY: EKEV, WE, WMU, F2, SETRC, ELORC

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(DP) :: enold, WEIGHT

    ELORC(S)=0.
    ENOLD=SETRC(S)
    SETRC(S)=0.
    DO I=2,NR
      DO K=2,NE
        DO L=2,NPA
          DO J=1,NT-1
            WEIGHT=F2(S,I,J,K,L)*WE(S,K)*WMU(L)
            SETRC(S)=SETRC(S)+EKEV(S,K)*WEIGHT
          END DO
        END DO
      END DO
    END DO
    ELORC(S)=ENOLD-SETRC(S)

    RETURN
  END SUBROUTINE SUMRC

!*************************************************************************
!                               ANISCH
!               Calculate pressure in equatorial plane
!*************************************************************************
  SUBROUTINE ANISCH(S)
    ! Module Variables
    use ModRamMain,      ONLY: DP, PathRamOut
    use ModRamConst,     ONLY: CS, PI, Q
    use ModRamParams,    ONLY: DoUseWPI, DoUsePlasmasphere, DoUseBASdiff, DoUseEMIC
    use ModRamGrids,     ONLY: NR, NT, NPA, ENG, SLEN, NCO, NCF, Ny, NE, &
                               ENG_emic, NCF_emic
    use ModRamTiming,    ONLY: Dt_bc, T, TimeRamNow
    use ModRamVariables, ONLY: IP1, UPA, WMU, FFACTOR, MU, EKEV, LZ, MLT, PHI,  &
                               PAbn, RMAS, GREL, IR1, EPP, ERNH, CDAAR, BDAAR,  &
                               ENOR, NDAAJ, fpofc, Kp, BOUNHS, FNHS, BNES, XNE, &
                               NECR, PPerT, PParT, F2, ATAW, ATAC, species, AE, &
                               ATAW_emic_h, ATAW_emic_he, EKEV_emic, fp2c_emic, &
                               Daa_emic_h, Daa_emic_he
    ! Module Subroutines/Functions
    use ModRamGSL,       ONLY: GSL_Interpolation_1D, GSL_Interpolation_2D
    use ModRamFunctions, ONLY: ACOSD, RamFileName
    use ModRamWPI,       ONLY: I_emic
    
    implicit none

    integer, intent(in) :: S
    integer :: i, iwa, j, k, klo, l, iz, kn, i1, j1, GSLerr, u
    real(DP) :: cv, rfac, rnhtt, edent, pper, ppar, rnht, Y, &
                eden,sume,suma,sumn,ernm,epma,epme,anis,epar,taudaa, &
                gausgam,anist,epart,fnorm,xfrl,Bw,esu,omega,er1,taudaa_h, taudaa_he
    real(DP) :: MUBOUN, IBw_Hband, IBw_Heband, fnorm_h, fnorm_he
    real(DP), ALLOCATABLE :: DWAVE(:),CMRA(:),BWAVE(:,:),AVDAA(:),TAVDAA(:),&
                             DAA(:,:,:),DUMP(:,:),XFR(:,:),XFRe(:),ALENOR(:),&
                             DUME(:,:),DVV(:,:,:),AVDVV(:),TAVDVV(:),WCDT(:,:),&
                             XFRT(:,:),PA(:),DAMR(:,:),DAMR1(:),GREL_new(:), &
                             DUMP2_h(:,:), DUMP2_he(:,:), logEkeV_emic(:)
    INTEGER :: KHI(5)
    character(len=2)  :: ST2
    character(len=214) :: NameFileOut
    character(len=2), dimension(4) :: speciesString = (/'_h','_o','he','_e'/)
    DATA khi/6, 10, 25, 30, 35/ ! ELB=0.1 keV -> 0.4,1,39,129,325 keV
!    DATA khi/2, 19, 28, 32, 35/ ! ELB=0.1 keV -> 0.1,10,100,200,427 keV

    ALLOCATE(DWAVE(NPA),CMRA(SLEN),BWAVE(NR,NT),AVDAA(NPA),TAVDAA(NPA), &
             DAA(NE,NPA,Slen),DUMP(ENG,NCF),XFR(NR,NT),XFRe(NCF),ALENOR(ENG), &
             DUME(ENG,NCF),DVV(NE,NPA,Slen),AVDVV(NPA),TAVDVV(NPA),WCDT(NR,NT), &
             XFRT(NR,NT),PA(NPA),DAMR(NPA,NCO),DAMR1(NPA),GREL_new(Ny), &
             DUMP2_h(ENG_emic, NCF_emic), DUMP2_he(ENG_emic,NCF_emic), &
             logEkeV_emic(ENG_emic))
    DWAVE = 0.0; CMRA = 0.0; BWAVE = 0.0; AVDAA = 0.0; TAVDAA = 0.0; DAA = 0.0
    DUMP = 0.0; XFR = 0.0; XFRe = 0.0; ALENOR = 0.0; DUME = 0.0; DVV = 0.0
    AVDVV = 0.0; TAVDVV = 0.0; WCDT = 0.0; XFRT = 0.0; PA = 0.0; DAMR = 0.0
    DAMR1 = 0.0; GREL_new = 0.0

    ST2 = speciesString(S)
    cv=CS*100 ! speed of light in [cm/s]
    esu=Q*3E9 ! elementary charge in [esu]
    RFAC=4*pi/cv
    gausgam=1.E-5
    khi(5)=NE

    ! calculate ring current parameters
    DO I=2,NR
      DO J=1,NT
        IF (DoUsePlasmasphere) THEN
!          XNE(I,J)=NECR(I1,J1)          ! CRasm model
          XNE(I,J)=NECR(I,J)
            if (T.gt.0.and.XNE(I,J).le.0) then
              write(*,*) 'in ANISCH T,S,i,j, XNE= ', T/3600,S,i,j, XNE(I,J)
              stop
            endif
        ELSE
          !if (S.EQ.1.and.I.eq.2.and.J.eq.1) write (*,*) " Need to specify Ne if using WPI"
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
            u = int(UPA(I)-1,kind=4)
            DO L=1,u
              ERNM=WMU(L)/FFACTOR(S,I,K,L)/FNHS(I,J,L)
              EPMA=ERNM*MU(L)*MU(L)
              EPME=ERNM-EPMA
              SUME=SUME+F2(S,I,J,K,L)*EPME
              SUMA=SUMA+F2(S,I,J,K,L)*EPMA
              SUMN=SUMN+F2(S,I,J,K,L)*ERNM
            ENDDO
            PPER=PPER+EPP(S,K)*SUME
            PPAR=PPAR+EPP(S,K)*SUMA
            RNHT=RNHT+ERNH(S,K)*SUMN
            EDEN=EDEN+ERNH(S,K)*EKEV(S,K)*SUMN
          ENDDO
          ANIS=PPER/2./PPAR-1.
          EPAR=2*PPAR/RNHT
          RNHT=RNHT*RFAC
          EDEN=EDEN*RFAC
          PPAR=2*RFAC*PPAR
          PPER=RFAC*PPER
          ! write anisotropy, kT parallel, RC density
          klo = khi(iwa)+1
          EDENT=EDENT+EDEN
          PPERT(S,I,J)=PPERT(S,I,J)+PPER
          PPART(S,I,J)=PPART(S,I,J)+PPAR
          RNHTT=RNHTT+RNHT
        enddo
        ANIST=PPERT(S,I,J)/PPART(S,I,J)-1.
        EPART=PPART(S,I,J)/RNHTT
!          if (T.gt.0.and.PPERT(S,I,J).lt.0) then
          if (abs(PPERT(S,I,J))>huge(PPERT(S,I,J))) then
              write(*,*) 'in ANISCH: S,i,j, PPERT= ', S,i,j, PPERT(S,I,J)
              stop
          endif
      ENDDO
    ENDDO

    IF (MOD(INT(T),INT(Dt_bc)).EQ.0.and.DoUseWPI.and.species(S)%s_name.eq.'Electron') THEN
         ! zero PA diffusion coefficients
       DO I=1,NR
          DO J=1,NT
             DO K=1,NE
                DO L=1,NPA
                   ATAW(I,J,K,L)=0.0 ! hiss
                   ATAC(I,J,K,L)=0.0 ! chorus
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ! linear interpolation for CHORUS b-aver diff coefficients

       DO L=1,NPA
          PA(L)=ACOSD(MU(nPa-L+1))
          DO IZ=1,NCO
             DAMR(L,IZ)=0.
          ENDDO
          DAMR1(L)=0.
       ENDDO
       DO I=2,NR
          DO J=1,NT
             IF (XNE(I,J).LE.50.) THEN ! outside pp
                xfrl=CS*SQRT(XNE(I,J)*RMAS(S)*40*PI)/10./BNES(I,J) ! Fpe/Fcyc
                fnorm=1 ! b-av, no Bw-dep
                ! interpolate PA diffusion coefficients for implicit scheme
                DO K=2,NE
                   DO L=1,NPA
                      IF (DoUseBASdiff) THEN
                         DAMR1(L)=LOG10(CDAAR(I,J,K,nPa-L+1))
                      ELSE
                         DAMR1(L)=LOG10(BDAAR(I,J,K,L))
                      ENDIF
                   ENDDO
                   DO L=1,NPA
                      MUBOUN=MU(L)+0.5*WMU(L)
                      CALL GSL_Interpolation_1D(PA,DAMR1,PAbn(L),Y,GSLerr)
                      taudaa=10.**Y*fnorm ! <Daa/p2> [1/s]
                      if (taudaa.gt.1e0) then
                         print*,'taudaa=',taudaa,' L=',LZ(I),' MLT=',MLT(J)
                         taudaa=1e-1
                      endif
                      if (taudaa.lt.1e-30) taudaa=1e-30
                      DWAVE(L)=taudaa*(1.-MUBOUN*MUBOUN)*MUBOUN*BOUNHS(I,J,L)
                      ATAC(I,J,K,L)=DWAVE(L)      ! call WPADIF twice, implicit
                   END DO
                END DO
             ENDIF
          ENDDO
       ENDDO

      ! Wave amplitude for plasmaspheric hiss
      Bw=30.
      IF(KP.GE.4.0) Bw=100.                                 ! pT
      ! linear interpolation for HISS normalized diff coefficients
      DO KN=1,ENG
         ALENOR(KN)=LOG10(ENOR(KN))
         DO IZ=1,NCF
           DUMP(KN,IZ)=0.
         ENDDO
      ENDDO
      DO I=2,NR
         DO J=1,NT
            IF (XNE(I,J).GT.50.) THEN ! inside pp
               omega=esu*10*BNES(I,J)/(RMAS(S)*cv) ! Fcyc, cgs
               xfrl=CS*SQRT(XNE(I,J)*RMAS(S)*40*PI)/10./BNES(I,J) ! Fpe/Fcyc
               if (xfrl.gt.18) xfrl=18. ! upper interp bound
               if (xfrl.lt.2) xfrl=2. ! lower interp bound
               fnorm=omega*(Bw*1e-3)**2*gausgam**2/1e8/BNES(I,J)/BNES(I,J)
               DO L=1,NPA
                  MUBOUN=MU(L)+0.5*WMU(L)
                  DO IZ=1,NCF
                     DO KN=1,ENG
                        DUMP(KN,IZ)=LOG10(NDAAJ(I,KN,L,IZ))
                     ENDDO
                  ENDDO
                  DO K=2,NE
                     ER1=LOG10(EKEV(S, K))
                     CALL GSL_Interpolation_2D(ALENOR,fpofc,DUMP,ER1,xfrl,Y,GSLerr)
                     DWAVE(L)=10.**Y*fnorm/GREL(S,K)**2*(1.-MUBOUN*MUBOUN)/MUBOUN ! denorm pa [1/s]
                     ATAW(I,J,K,L)=DWAVE(L)           ! call WPADIF twice, implicit
                     taudaa=dwave(l)/MUBOUN/BOUNHS(I,J,L) ! <Dmu>
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
    ENDIF ! end diff coeff loop 


    IF (MOD(INT(T),INT(Dt_bc)).EQ.0.and.DoUseEMIC.and.species(S)%s_name.eq.'Hydrogen')THEN
       ! only need to calculate once to obtain the ATAW_emic coefficient
       !.......zero PA diffusion coefficients
       DO  I=1,NR
          DO J=1,NT
             DO K=1,NE
                DO L=1,NPA
                   ATAW_emic_h(I,J,K,L)=0.0                        ! EMIC 
                   ATAW_emic_he(I,J,K,L)=0.0                       ! EMIC 
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       

!       write(NameFileOut1,'(a,a,a,i6.6,a)')&
!            PathRamOut//"ram_",ST2,"_t",nint(TimeRamElapsed/300._Real8_),".wde"
!       open(unit=24,file=trim(namefileout1),status='unknown')
!       write(24,555) T/3600, Kp
!555    Format(2x,3HT =, F8.3, 2X, 4HKp =, F6.2, '  EMIC diff coeff')

       DO KN=1, ENG_emic
          logEkeV_emic(KN) = log10(EkeV_emic(KN))
          DO IZ=1,NCF_emic
             DUMP2_h(KN,IZ)=0.
             DUMP2_he(KN,IZ)=0.
          END DO
       END DO

       DO I=2, NR
          DO J=1, NT
             xfrl = CS*sqrt(XNE(I,J)*RMAS(4)*40*PI)/10./BNES(I,J) !fpe/fce electrons,RMAS(4): electron
             if (xfrl .gt.20)xfrl = 20.
             if (xfrl .lt.2) xfrl = 2.

             call I_emic(AE, I, J, IBw_Hband, IBw_Heband)              ! Bw intensity (nT^2)
             
             fnorm_h = IBw_Hband  ! normalized by wave amplitude
             fnorm_he= IBw_Heband ! normalized by 1nT 

!             write(24,554)LZ(I), MLT(J), xfrl, XNE(I,J), &
!                  'IBw Hband=',IBw_Hband, 'IBw Heband', IBw_Heband,' Fpe/Fcyc=',xfrl

             DO L=1,NPA
                MUBOUN=MU(L)+0.5*WMU(L)
                
                DO IZ=1,NCF_emic
                   DO KN=1,ENG_emic
                      DUMP2_h(KN,IZ) = log10(Daa_emic_h(I,KN,L,IZ))
                      DUMP2_he(KN,IZ) = log10(Daa_emic_he(I,KN,L,IZ))
                   END DO
                END DO
                DO K=2,NE
                   ER1 = log10(EKEV(S,K))
                   ! EkeV_emic: 0.1 keV to 1000 keV 
                   call GSL_Interpolation_2D(logEkeV_emic, fp2c_emic, &
                        DUMP2_h, ER1, xfrl, Y, GSLerr)
                   DWAVE(L) = 10.**Y*fnorm_h * (1.-MUBOUN*MUBOUN) &
                        * MUBOUN*BOUNHS(I,J,L) ! Daa*(1-mu^2)*mu*h = Duu*mu*h = Dwave 
                   ATAW_emic_h(I,J,K,L) = DWAVE(L)
                   
                   call GSL_Interpolation_2D(logEkeV_emic, fp2c_emic, &
                        DUMP2_he, ER1, xfrl, Y, GSLerr)
                   DWAVE(L) = 10.**Y*fnorm_he * (1.-MUBOUN*MUBOUN) &
                        * MUBOUN*BOUNHS(I,J,L) ! Daa*(1-mu^2)*mu*h = Duu*mu*h = Dwave 
                   ATAW_emic_he(I,J,K,L) = DWAVE(L)

                   if (ATAW_emic_h (I,J,K,L) .le. 1.0e-20)ATAW_emic_h (I,J,K,L) = 1.0e-31
                   if (ATAW_emic_he(I,J,K,L) .le. 1.0e-20)ATAW_emic_he(I,J,K,L) = 1.0e-31

                   taudaa_h  = ATAW_emic_h(I,J,K,L)/MUBOUN/BOUNHS(I,J,L)/(1.-MUBOUN*MUBOUN)
                   taudaa_he = ATAW_emic_he(I,J,K,L)/MUBOUN/BOUNHS(I,J,L)/(1.-MUBOUN*MUBOUN)
                   
!                   write(24,556)ekev(k), PABn(L), taudaa_h, taudaa_he, &
!                        1./taudaa_h, 1./taudaa_he
                END DO
             END DO
          END DO
       END DO

!554    format(3H L=,F5.2,5H MLT=,F5.2,9H fpe/fce=,E10.3, &
!            4H Ne=,F9.2,A10,E10.3,A11, E10.3, A10, E10.3)
!556    format(2F9.3,6(1PE12.3))                  
!       close(24)
       
    END IF
    
    DEALLOCATE(DWAVE,CMRA,BWAVE,AVDAA,TAVDAA, &
               DAA,DUMP,XFR,XFRe,ALENOR, &
               DUME,DVV,AVDVV,TAVDVV,WCDT, &
               XFRT,PA,DAMR,DAMR1,GREL_new)
  RETURN
  end subroutine ANISCH

END MODULE ModRamRun
