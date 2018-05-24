!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamRun
! Contains the subroutine for running ram and calculating pressure
! Calls subroutines from ModRamDrift, ModRamWPI, and ModRamLoss
! Calculates new F2 and parallel/perpendicular pressures

  use ModRamVariables, ONLY: SETRC, ELORC, LSDR, LSCHA, LSATM, LSCOE, LSCSC, &
                             LSWAE, XNN, XND, LNCN, LNCD, LECN, LECD, ENERN, &
                             ENERD, ATEW, ATAW, ATAC, ATEC, ATMC, ATAW_EMIC, &
                             NECR, PParH, PPerH, PParO, PPerO, PParE, PPerE, &
                             PParHe, PPerHe, PParT, PPerT, F2

  implicit none
  save

  contains

!==============================================================================
  SUBROUTINE ram_run

    !!!! Module Variables
    use ModRamMain,      ONLY: Real8_
    use ModRamParams,    ONLY: electric, DoUseWPI, DoUseBASdiff, DoUseKpDiff, &
                               DoUsePlane_SCB
    use ModRamConst,     ONLY: RE
    use ModRamGrids,     ONLY: NR, NE, NT, NPA
    use ModRamTiming,    ONLY: UTs, T, DtEfi, Dt_hI, DtsMin, DtsNext, TimeRamNow, &
                               Dts
    use ModRamVariables, ONLY: Kp, VT, VTOL, VTN, TOLV, LZ, PHI, PHIOFS, MU, &
                               WMU, FFACTOR, FLUX, FNHS, DtDriftR, DtDriftP, &
                               DtDriftE, DtDriftMu, NECR, F107, RZ, IG, rsn, &
                               IAPO
    !!!! Module Subroutines/Functions
    use ModRamDrift, ONLY: DRIFTPARA, DRIFTR, DRIFTP, DRIFTE, DRIFTMU, DRIFTEND, &
                           CDriftR, CDriftP, CDriftE, CDriftMu
    use ModRamLoss,  ONLY: CEPARA, CHAREXCHANGE, ATMOL
    use ModRamWPI,   ONLY: WAPARA_KP, WPADIF, WAVELO
    !!!! Share Modules
    use ModTimeConvert, ONLY: TimeType

    implicit none

    real(kind=Real8_)  :: AVS, VT_kV(NR+1,NT), lambda
    integer :: i, j, k, l, jw ! Iterators
    integer :: IS ! Species to advect.
    integer :: DAY, nmonth

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
    IF (DoUsePlane_SCB) THEN
       VT_kV = VT/1e3 ! potential in kV for call to plane
       ! Need total number of days from TimeRamNow
       DAY = TimeRamNow%iMonth*30 + TimeRamNow%iDay
       CALL APFMSIS(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,TimeRamNow%iHour,IAPO)
       CALL TCON(TimeRamNow%iYear,TimeRamNow%iMonth,TimeRamNow%iDay,DAY,RZ,IG,rsn,nmonth)
       CALL PLANE(TimeRamNow%iYear,DAY,T,KP,IAPO(2),RZ(3),F107,2.*DTs,NECR,VT_kV)
    END IF

    do iS=1,4
       !   calls for given species S:
       CALL CEPARA(iS)
       CALL DRIFTPARA(iS)

       ! Call routines to calculate the changes of distribution function
       ! considering drifts, charge exchange and atmospheric losses
       CALL DRIFTR(iS)
       CALL DRIFTP(iS)
       CALL DRIFTE(iS)
       CALL DRIFTMU(iS)

       CALL SUMRC(iS)
       LSDR(iS)=LSDR(iS)+ELORC(iS)

       if (iS.GT.1) then
          CALL CHAREXCHANGE(iS)
          CALL SUMRC(iS)
          LSCHA(iS)=LSCHA(iS)+ELORC(iS)
       else
          IF (DoUseWPI) THEN
             CALL WPADIF(iS) ! pitch angle diffusion
          ELSE
             CALL WAVELO(iS) ! electron lifetimes
          ENDIF
          CALL SUMRC(iS)
          LSWAE(iS)=LSWAE(iS)+ELORC(iS)
       endif

       CALL ATMOL(iS)
       CALL SUMRC(iS)
       LSATM(iS)=LSATM(iS)+ELORC(iS)
       ! time splitting
       CALL ATMOL(iS)
       CALL SUMRC(iS)
       LSATM(iS)=LSATM(iS)+ELORC(iS)

       if (iS.GT.1) then
          CALL CHAREXCHANGE(iS)
          CALL SUMRC(iS)
          LSCHA(iS)=LSCHA(iS)+ELORC(iS)
       else
          IF (DoUseWPI) THEN
             CALL WPADIF(iS)
          ELSE
             CALL WAVELO(iS)
          ENDIF
          CALL SUMRC(iS)
          LSWAE(iS)=LSWAE(iS)+ELORC(iS)
       endif

       CALL DRIFTMU(iS)
       CALL DRIFTE(iS)
       CALL DRIFTP(iS)
       CALL DRIFTR(iS)

       CALL SUMRC(iS)
       LSDR(iS)=LSDR(iS)+ELORC(iS)

       ! Interpolate diffusivity as funtion of Kp, if necessary
       if (DoUseWPI.and.DoUseBASdiff.and.DoUseKpDiff) then
          call WAPARA_Kp()
       end if

       ! Routines needed to clear allocated variables for use with OpenMP
       CALL DRIFTEND
    end do

    DtsNext = min(minval(DtDriftR), minval(DtDriftP), minval(DtDriftE), minval(DtDriftMu))
    DtsNext = max(DtsNext, DtsMin)

    ! Update flux and pressure totals
    DO iS = 1,4
       CALL ANISCH(iS)
       DO I = 2, NR
          DO K = 2, NE
             DO L = 2, NPA
                DO J = 1, NT-1
                   FLUX(iS,I,J,K,L) = F2(iS,I,J,K,L)/FFACTOR(iS,I,K,L)/FNHS(I,J,L)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Update species pressures
    PPerO  = PPerT(4,:,:); PParO  = PParT(4,:,:)
    PPerHe = PPerT(3,:,:); PParHe = PParT(3,:,:)
    PPerE  = PPerT(1,:,:); PParE  = PParT(1,:,:)
    PPerH  = PPerT(2,:,:); PParH  = PParT(2,:,:)

    RETURN
  END SUBROUTINE Ram_Run

!************************************************************************
!                       SUMRC     
!               Calculate the total energy & energy loss 
!************************************************************************
  SUBROUTINE SUMRC(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: Real8_
    use ModRamGrids,     ONLY: NR, NE, NPA, NT
    use ModRamVariables, ONLY: EKEV, WE, WMU, F2

    implicit none
    save

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(kind=Real8_):: enold
    real(kind=Real8_):: WEIGHT
    !$OMP THREADPRIVATE(ENOLD, WEIGHT, I, J, K, L)

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
  END SUBROUTINE SUMRC

!*************************************************************************
!                               ANISCH
!               Calculate pressure in equatorial plane
!*************************************************************************
  SUBROUTINE ANISCH(S)
    ! Module Variables
    use ModRamMain,      ONLY: Real8_
    use ModRamConst,     ONLY: CS, PI, Q
    use ModRamParams,    ONLY: DoUseWPI, DoUsePlane_SCB, DoUseBASdiff
    use ModRamGrids,     ONLY: NR, NT, NPA, ENG, SLEN, NCO, NCF, Nx, Ny, NE
    use ModRamTiming,    ONLY: Dt_bc, T
    use ModRamVariables, ONLY: IP1, UPA, WMU, FFACTOR, MU, EKEV, LZ, PHI, MLT,  &
                               PAbn, RMAS, GREL, IR1, EPP, ERNH, CDAAR, BDAAR,  &
                               ENOR, NDAAJ, fpofc, Kp, BOUNHS, FNHS, BNES, XNE, &
                               NECR
    ! Module Subroutines/Functions
    use ModRamGSL,       ONLY: GSL_Interpolation_1D, GSL_Interpolation_2D
    use ModRamFunctions, ONLY: ACOSD

    implicit none
    save

    integer, intent(in) :: S
    integer :: i, iwa, j, k, klo, l, iz, ier, kn, i1, j1, GSLerr, u
    real(kind=Real8_) :: cv, rfac, rnhtt, edent, pper, ppar, rnht, Y, &
                         eden,sume,suma,sumn,ernm,epma,epme,anis,epar,taudaa,taudea,taudee, &
                         gausgam,anist,epart,fnorm,xfrl,Bw,esu,omega,er1,dx
    real(kind=Real8_) :: MUBOUN,KEVERG,DN(2,25),EPO(2,25),AII(2,25)
    real(kind=Real8_), ALLOCATABLE :: DWAVE(:),CMRA(:),BWAVE(:,:),AVDAA(:),TAVDAA(:),&
                                      DAA(:,:,:),DUMP(:,:),XFR(:,:),XFRe(:),ALENOR(:),&
                                      DUME(:,:),DVV(:,:,:),AVDVV(:),TAVDVV(:),WCDT(:,:),&
                                      XFRT(:,:),PA(:),DAMR(:,:),DAMR1(:),GREL_new(:),BOUNHS_(:)
    INTEGER :: KHI(5)
    INTEGER, ALLOCATABLE :: MINP(:), MAXP(:)
    character(len=80) HEADER
    DATA khi/6, 10, 25, 30, 35/ ! ELB=0.1 keV -> 0.4,1,39,129,325 keV 

    ALLOCATE(DWAVE(NPA),CMRA(SLEN),BWAVE(NR,NT),AVDAA(NPA),TAVDAA(NPA), &
             DAA(NE,NPA,Slen),DUMP(ENG,NCF),XFR(NR,NT),XFRe(NCF),ALENOR(ENG), &
             DUME(ENG,NCF),DVV(NE,NPA,Slen),AVDVV(NPA),TAVDVV(NPA),WCDT(NR,NT), &
             XFRT(NR,NT),PA(NPA),DAMR(NPA,NCO),DAMR1(NPA),GREL_new(Ny),BOUNHS_(Nx))
    ALLOCATE(MINP(nT),MAXP(nT))

    cv=CS*100 ! speed of light in [cm/s]
    esu=Q*3E9 ! elementary charge in [esu]
    RFAC=4*pi/cv
    gausgam=1.E-5
    khi(5)=NE
    ! calculate ring current parameters
    DO I=2,NR
      I1=(I-2)*IR1+3
      DO J=1,NT
        J1=(J-1)*IP1
        IF (DoUsePlane_SCB) THEN
          XNE(I,J)=NECR(I1,J1)
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
            u = UPA(I)-1
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
            EDEN=EDEN+ERNH(S,K)*EKEV(K)*SUMN
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
     ENDDO
  ENDDO

    IF (MOD(INT(T),INT(Dt_bc)).EQ.0.and.DoUseWPI.and.S.eq.1) THEN
         ! zero PA diffusion coefficients
       DO I=1,NR
          DO J=1,NT
             DO K=1,NE
                DO L=1,NPA
                   ATAW(I,J,K,L)=0.0 ! hiss
                   ATAW_emic(I,J,K,L)=0.0 ! EMIC
                   ATAC(I,J,K,L)=0.0 ! chorus
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       ! linear interpolation for CHORUS b-aver diff coefficients

       DO L=1,NPA
          PA(L)=ACOSD(MU(L))
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
                         DAMR1(L)=LOG10(CDAAR(I,J,K,L))
                      ELSE
                         DAMR1(L)=LOG10(BDAAR(I,J,K,L))
                      ENDIF
                   ENDDO
                   DO L=1,NPA
                      MUBOUN=MU(L)+0.5*WMU(L)
                      CALL GSL_Interpolation_1D('Steffen',PA,DAMR1,PAbn(L),Y,GSLerr)
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
                     ER1=LOG10(EKEV(K))
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

    DEALLOCATE(DWAVE,CMRA,BWAVE,AVDAA,TAVDAA, &
               DAA,DUMP,XFR,XFRe,ALENOR, &
               DUME,DVV,AVDVV,TAVDVV,WCDT, &
               XFRT,PA,DAMR,DAMR1,GREL_new,BOUNHS_)
    DEALLOCATE(MINP,MAXP)
  RETURN
  end subroutine ANISCH

END MODULE ModRamRun
