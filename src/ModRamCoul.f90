!============================================================================
!    Copyright (c) 2016, Los Alamos National Laboratory
!    All rights reserved.
!============================================================================

MODULE ModRamCoul
! Contains subroutines for calculating energy and PA losses due to Coulomb collisions

  implicit none

  contains

!**************************************************************************
!                              COULPARA
!   Calculates the time-independent parameters used in the Coulomb terms
!**************************************************************************
  SUBROUTINE COULPARA(S)
    !!!! Module Variables
    use ModRamConst,     ONLY: Q, PI, CS, RE, MP
    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NE, NR, NT, NPA, NS
    use ModRamTiming,    ONLY: DTs
    use ModRamVariables, ONLY: RMAS, species, VBND, V, GREL, COULE, COULI, ATA, &
                               GTA, GTAE, GTAI, CEDR, CIDR, MU, WMU, DMU, EKEV, GRBND
    use ModRamSpecies,   ONLY: RAMSpecies, nSpecies
    use ModRamFunctions, ONLY: FUNT, FUNI, Gcoul, ERFF

    implicit none

    integer, intent(in) :: S
    integer :: l, k, i, j, is
    real(DP) :: EPS, DLN, QE, GAMA, CCO, CCD, EDRCO, EDRE, EDRI
    real(DP) :: CCE, CDE, CCI, CDI, CCDE, CCDI
    real(DP) :: x, xd, MUBOUN, BADIF, AFIR, AFER, ASIC, ASEC, RA
    real(DP) :: Zt, BANE, VF, Zb

    real(DP), allocatable :: COULDE(:), COULDI(:)

    allocate(COULDE(nPa), COULDI(nPa))
    COULDE = 0.0; COULDI = 0.0
    CCE = 0._dp
    CDE = 0._dp
    EDRE = 0._dp
    CCI = 0._dp
    CDI = 0._dp
    EDRI = 0._dp

! Proportions of species in plasmasphere
!old    DATA RA/1.0,0.77,0.2,0.03/    ! e-, H+, He+, O+
!        DATA RA/0.77,0.03,0.2,1.0/    ! H+, O+, He+, e- 
! These proportions were moved to RAMSpecies and can be edited there


!...Parameters used in calculating Coulomb collisions
    EPS=8.854E-12                   ! [F/m]
    DLN=21.5                        ! Coulomb logarithm
    Zt=species(S)%s_charge          ! charge of RC species
    QE=(Q**2/EPS)
    GAMA=Zt**2*DLN/4./PI*QE*1E6*QE  ! *1E6 to transform NE in [m-3]
    CCO=GAMA/Q*DTs/Q/1E3               ! /Q/1E3 to transform WE in [J]
    CCD=GAMA*DTs/RMAS(S)**2/CS**3      ! DT- time-splitting, implicit, CHECK coeff!!!!
    EDRCO=DLN*QE*RE/RMAS(S)*QE*1E9/Q   ! in [eV/cm2/s], check!

    DO k= 1, nE
!...Calculate the thermal velocity of each plasmasph species
!   We assume: Te=Ti=1eV (kT=1eV)
       do iS = 1, nSpecies
          RA = RAMSpecies(iS)%plasmasphereRatio
          if (RA < 1e-9) cycle
          VF = SQRT(2.*Q/(MP*RAMSpecies(iS)%s_mass))   ! [m/s]
          Zb = RAMSpecies(iS)%s_charge                 ! |charge| of plasmasphere species
          X=VBND(S,K)/VF
          XD=V(S,K)/VF
          if (Zb < 0._dp) then  ! Negative charge (electrons)
             CCE = CCE + RA*Gcoul(X)
             CDE = CDE + RA*(ERFF(XD)-Gcoul(XD))
             EDRE = EDRE + RA*Gcoul(XD)
          else                  ! Non-negative charge (ions)
             CCI = CCI + RA*Zb**2*Gcoul(X)
             CDI = CDI + RA*Zb**2*(ERFF(XD)-Gcoul(XD))
             EDRI = EDRI + RA*Gcoul(XD)
          endif
       enddo

!...Collisions with plasmaspheric electrons (minus: dF/dt-d(cF)/dE=0)
       COULE(S,K,1)=-CCE*VBND(S,K)*CCO*GRBND(S,K)**2
!...Collisions with plasmaspheric ions
       COULI(S,K,1)=-CCI*VBND(S,K)*CCO*GRBND(S,K)**2
!...Energy deposition rate [eV/cm2/s] in 1 hemisphere, check:
       CEDR(S,K,1)=EDRCO*EKEV(K)*EDRE*(GREL(S,K)+1)/GREL(S,K)**2/V(S,K)**2
       CIDR(S,K,1)=EDRCO*EKEV(K)*EDRI*(GREL(S,K)+1)/GREL(S,K)**2/V(S,K)**2
!...Coulomb scattering 
       CCDE=CCD*CDE*GREL(S,K)/(GREL(S,K)**2-1)**(1.5)
       CCDI=CCD*CDI*GREL(S,K)/(GREL(S,K)**2-1)**(1.5)
       DO L=2,NPA-1
          BANE=(1.-FUNI(MU(L))/2./FUNT(MU(L)))/(1.-MU(L)*MU(L))   ! assume Bdip
          COULE(S,K,L)=COULE(S,K,1)
          COULI(S,K,L)=COULI(S,K,1)
          CEDR(S,K,L)=CEDR(S,K,1)*BANE*FUNT(MU(L))*MU(L)*WMU(L)   ! use in WRESULT
          CIDR(S,K,L)=CIDR(S,K,1)*BANE*FUNT(MU(L))*MU(L)*WMU(L)   ! use in WRESULT
          MUBOUN=MU(L)+0.5*WMU(L)
          BADIF=(1.-MUBOUN*MUBOUN)/MUBOUN/2.

!...Collisions with plasmaspheric electrons
          COULDE(L)=CCDE*BADIF
          AFER=COULDE(L)/MU(L)/DMU(L)/WMU(L)
          ASEC=COULDE(L-1)/MU(L)/DMU(L-1)/WMU(L)

!...Collisions with plasmaspheric ions
          COULDI(L)=CCDI*BADIF
          AFIR=COULDI(L)/MU(L)/DMU(L)/WMU(L)
          ASIC=COULDI(L-1)/MU(L)/DMU(L-1)/WMU(L)

!...
          ATA(S,K,L)=AFIR+AFER
          GTA(S,K,L)=ASIC+ASEC
       END DO
       CEDR(S,K,NPA)=CEDR(S,K,NPA-1)
       CIDR(S,K,NPA)=CIDR(S,K,NPA-1)
       ATA(S,K,NPA)=0
    END DO

    DEALLOCATE(COULDE, COULDI)
    RETURN
  END SUBROUTINE COULPARA


!**************************************************************************
!                        COULEN
!        Routine calculates the change of distribution function due to
!        Coulomb energy degradation
!**************************************************************************
    SUBROUTINE COULEN(S)

    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NE, NT, NR, NPA, nS
    use ModRamVariables, ONLY: F2, EKEV, WE, DE, RMAS, NECR, GREL, COULE, COULI, &
                           IR1, IP1, XNE, MU, FNHS, FNIS
    use ModRamTiming,    ONLY: TimeRamElapsed

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l, j1, i1, ISIGN, n
    real(DP) :: T, ezero, x, fup, r, corr, LIMITER
    real(DP), ALLOCATABLE :: FBND(:),F(:),GRZERO(:),CccolE(:),BANE(:)

    ALLOCATE(GRZERO(nS),FBND(nE),F(0:nE+2),CccolE(nE),BANE(nPa))
    GRZERO = 0.0; FBND = 0.0; F = 0.0; CccolE = 0.0; BANE = 0.0

    T=TimeRamElapsed

        EZERO=EKEV(1)-WE(1)
        GRZERO(S)=1.+EZERO*1000.*Q/RMAS(S)/CS/CS
        F(NE+1)=0.
        F(NE+2)=0.
        DO 1 J=1,NT
           J1=(J-1)*IP1          
        DO 1 I=2,NR
           I1=(I-2)*IR1+3
!.......equatorial plasmaspheric density [cm-3]: NECR(I1,J1)
!       NEED bounce-averaged plasmasphere density: XNE(I,J,L)=<NEs>=NECR(I1,J1)*BANE(L)
         DO L=2,NPA-1
           BANE(L)=(1.-FNIS(I,J,L)/2./FNHS(I,J,L))/(1.-MU(L)*MU(L))
         END DO
         DO L=NPA-10,NPA
           BANE(L)=BANE(L-1)
         END DO
        DO 1 L=2,NPA
         XNE(I,J)=NECR(I,J)*BANE(L)   ! assume n/B=const
          if (XNE(I,J) < 0._dp) then
              write(*,*) 'in COULEN XNE<0 ', T/3600, S,i,j,l, &
                   XNE(I,J),NECR(I,J),BANE(L)
          endif
         DO K=2,NE
          F(K)=F2(S,I,J,K,L)
         END DO

!.......fix low b.c.
         F(1)=F(2)*GREL(S,1)/GREL(S,2)*SQRT((GREL(S,1)**2-1)/(GREL(S,2)**2-1))
         F(0)=F(1)*GRZERO(S)/GREL(S,1)*SQRT((GRZERO(S)**2-1)/(GREL(S,1)**2-1))

         DO K=1,NE
          CccolE(K)=(COULE(S,K,L)+COULI(S,K,L))*XNE(I,J)
          ISIGN=1
          IF(CccolE(K) < 0._dp) ISIGN=-1
!          if (ABS(CccolE(K)/DE(K)).GT.1) then
!            open(20,file='CccolE.dat',status='unknown',position='append')
!              write(20,*) ' CccolE=', CccolE(K), S,i,j,k,l, COULE(S,K,L), COULI(S,K,L), BANE(L)
!              write(20,*) 'Thr= ', T/3600
!            close(20)
!          endif
          X=F(K+1)-F(K)
          FUP=0.5*(F(K)+F(K+1)-ISIGN*X)
          IF (ABS(X).LE.1.E-27) FBND(K)=FUP
          IF (ABS(X).GT.1.E-27) THEN
            N=K+1-ISIGN
             R=(F(N)-F(N-1))/X
            IF (R.LE.0) FBND(K)=FUP
            IF (R.GT.0) THEN
               LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
               CORR=-0.5*(CccolE(K)/DE(K)-ISIGN)*X   
               FBND(K)=FUP+LIMITER*CORR
            END IF
          END IF  
         END DO

         DO K=2,NE
           F2(S,I,J,K,L)=F2(S,I,J,K,L)-CccolE(K)/WE(K)*FBND(K)+CccolE(K-1)/WE(K)*FBND(K-1)
           if (f2(s,i,j,k,l) < 0._dp) then
!              write(*,*) 'in COULEN f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
              f2(S,i,j,k,l)=1E-15
           endif
         END DO
1        CONTINUE

    DEALLOCATE(GRZERO,FBND,F,CccolE,BANE)
    RETURN
  END SUBROUTINE COULEN


!************************************************************************
!                       COULMU
!        Routine calculates the change of the distribution function
!        due to Coulomb pitch angle scattering
!************************************************************************
  SUBROUTINE COULMU(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NE, NT, NR, NPA, nS
    use ModRamVariables, ONLY: F2, EKEV, WE, DE, RMAS, NECR, ATA, GTA, &
                           IR1, IP1, XNE, MU, WMU, DMU, FNHS, BOUNHS, FNIS, BOUNIS
    use ModRamTiming,    ONLY: TimeRamElapsed

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l, j1, i1
    real(DP) :: T, AN, BN, GN, RP, DENOM
    real(DP), ALLOCATABLE :: RK(:),RL(:),BASCNE(:,:,:)

    ALLOCATE(RK(nPa),RL(nPa),BASCNE(NR,NT,NPA))
    RK = 0.0; RL = 0.0

    T=TimeRamElapsed

        DO 1 J=1,NT
           J1=(J-1)*IP1          
        DO 1 I=2,NR
           I1=(I-2)*IR1+3  
!.......plasmaspheric density [cm-3]
!        XNE(I,J)=NECR(I1,J1)                           ! CRasm model
           XNE(I,J)=NECR(I,J)
        DO 1 K=2,NE
          RK(1)=0.                                   ! lower boundary condition
          RL(1)=-1.                                   ! "
          BASCNE(I,J,1) = XNE(I,J)*BOUNIS(I,J,1)/2./BOUNHS(I,J,1)
          DO L=2,NPA-1
! NEED b-aver: BASCNE=<NEs*Bo/Bs*(1-Bs/Bm)>
             BASCNE(I,J,L)=XNE(I,J)*BOUNIS(I,J,L)/2./BOUNHS(I,J,L)          ! assuming  n/B=const
             if ((T > 0._dp).and.(BASCNE(I,J,L) <= 0._dp)) then
                write(*,*) 'in COULMU BASCNE<0 ', T/3600,S,i,j,l, &
                   BASCNE(I,J,L),XNE(I,J),BOUNIS(I,J,L),BOUNHS(I,J,L)
             endif

             AN=ATA(S,K,L)*BASCNE(I,J,L)/FNHS(I,J,L)*BOUNHS(I,J,L)      ! only Coul scatt
             GN=GTA(S,K,L)*BASCNE(I,J,L-1)/FNHS(I,J,L)*BOUNHS(I,J,L-1)  !        "
             BN=AN+GN                                                   !        "
             RP=F2(S,I,J,K,L)/FNHS(I,J,L)/MU(L)                         ! implicit
             DENOM=BN+GN*RL(L-1)+1
             RK(L)=(RP+GN*RK(L-1))/DENOM
             RL(L)=-AN/DENOM
          END DO           

          F2(S,I,J,K,NPA-1)=RK(NPA-1)/(1+RL(NPA-1))        ! upper b.c.
          DO L=NPA-2,1,-1
             F2(S,I,J,K,L)=RK(L)-RL(L)*F2(S,I,J,K,L+1)
          END DO
          F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)
          DO L=1,NPA
             F2(S,I,J,K,L)=F2(S,I,J,K,L)*FNHS(I,J,L)*MU(L)
             if ((T > 0._dp).and.(f2(s,i,j,k,l) < 0._dp)) then
!                write(*,*) 'in COULMU f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
                f2(S,i,j,k,l)=1E-15
             endif
          ENDDO

1        CONTINUE

    DEALLOCATE(RK,RL,BASCNE)
    RETURN
  END SUBROUTINE COULMU

!==============================================================================
END MODULE ModRamCoul
