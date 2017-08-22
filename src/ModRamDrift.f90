MODULE ModRamDrift
! Contains subroutines responsible for calculating flux changes due to
! different drifts

  use ModRamVariables, ONLY: P1, P2, VR, EDOT, DtDriftR, DtDriftP, DtDriftE, &
                             DtDriftMu, FracCFL, MUDOT

  implicit none
  save

contains

!**************************************************************************
!                               OTHERPARA
!                       Calculate drift parameters
!**************************************************************************
  SUBROUTINE DRIFTPARA

    use ModRamMain,      ONLY: Real8_, S
    use ModRamConst,     ONLY: PI
    use ModRamGrids,     ONLY: NS, NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: RLZ, MDR, DPHI, EKEV, GREL, WMU, EBND, GRBND, &
                               PHIOFS, MU

    implicit none

    real(kind=Real8_) :: RA(NS), MUBOUN, QS
    integer :: i, j, k, l

    DATA RA/1.,.77,.2,.03/    ! Proportions of species in the plasmasph
    ! Electric field offset in radians and particle charge
    PHIOFS=0*PI/12.
    QS=1.
    IF (S.EQ.1) QS=-1.
    ! Parameters used in calculating radial and azimuthal drifts at boundaries
    DO I=1,NR
      DO J=1,NT
        VR(I)=DTs/MDR/(RLZ(I)+0.5*MDR)/2/DPHI
        ! Kp dependent part of azimuthal drift 
        P1(I)=DTs/DPHI/2/MDR/RLZ(I)
      END DO
      DO K=1,NE
        P2(I,K)=DTs*EKEV(K)*1000*(GREL(S,K)+1)/GREL(S,K)/RLZ(I)**2/DPHI/QS
      END DO
    END DO

    ! Pitch angle and energy time derivatives at boundaries of grids
    DO I=1,NR
      DO L=1,NPA-1
        MUBOUN=MU(L)+0.5*WMU(L)
        MUDOT(I,L)=(1.-MUBOUN**2)*DTs/2/MUBOUN/RLZ(I)
      END DO
      MUDOT(I,NPA)=0.
      DO K=1,NE
        EDOT(I,K)=EBND(K)*DTs/RLZ(I)*(GRBND(S,K)+1)/GRBND(S,K)/2.
      END DO
    END DO

    RETURN
  END

!************************************************************************
!                            DRIFTR
!               Calculate changes due to radial drift
!************************************************************************
  SUBROUTINE DRIFTR

    use ModRamMain,      ONLY: Real8_, S
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts, DtsNext
    use ModRamParams,    ONLY: BetaLim
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, MDR, EKEV, GREL, DPHI, &
                               RLZ, CONF1, CONF2, FGEOS, VT, EIP
    implicit none

    integer :: UR, i, j, j0, j1, k, l, n, ISGN(NR,NT)
    real(kind=Real8_) :: p4, qs, x, fup, r, corr, cgr1, cgr2, cgr3
    real(kind=Real8_) :: CGR(NR,NT,NE,NPA),RGR(NR,NT,NE,NPA)
    real(kind=Real8_) :: F(NR+2),FBND(NR),C(NR,NT),LIMITER,CR(NR,NT)

    DTDriftR = 100000.0
    QS=1.
    IF (S.EQ.1) QS=-1.
    DO I=1,NR
      DO J=1,NT
        J0=J-1
        IF (J.EQ.1) J0=NT-1
        J1=J+1
        IF (J.EQ.NT) J1=2
        CR(I,J)=VR(I)*(VT(I,J0)+VT(I+1,J0)-VT(I,J1)-VT(I+1,J1))/(BNES(I,J)+BNES(I+1,J)) &
               +(EIP(I,J)+EIP(I+1,J))/(BNES(I,J)+BNES(I+1,J))*DTs/MDR
      ENDDO
    ENDDO

    DO 1 K=2,NE
      P4=DTs*EKEV(K)*1000.0*(GREL(S,K)+1)/GREL(S,K)/DPHI/MDR/QS
      DO 1 L=2,NPA
        DO 1 J=1,NT !   < FBND(I) >
          J0=J-1 !  |____.____|____.____|____.____|
          IF (J.EQ.1) J0=NT-1 !  <   F(I)  >
          J1=J+1 !      F - average in cell(i,j,k,l)
          IF (J.EQ.NT) J1=2
          DO I=1,NR
            F(I)=F2(S,I,J,K,L)
            CGR1=FNIS(I+1,J1,L)+FNIS(I,J1,L)-FNIS(I+1,J0,L)-FNIS(I,J0,L)
            CGR2=BNES(I+1,J1)+BNES(I,J1)-BNES(I+1,J0)-BNES(I,J0)
            CGR3=CGR1+(FNIS(I+1,J,L)+FNIS(I,J,L)-2*FNHS(I+1,J,L) &
                -2*FNHS(I,J,L))*CGR2/2./(BNES(I+1,J)+BNES(I,J))
            CGR(I,J,K,L)=CGR3/(FNHS(I,J,L)+FNHS(I+1,J,L))*P4/2./(BNES(I,J)+BNES(I+1,J))/(RLZ(I)+0.5*MDR)
            C(I,J)=CR(I,J)+CGR(I,J,K,L)
            RGR(I,J,K,L)=C(I,J)*MDR/DTs
            DTDriftR = min( DTDriftR, FracCFL*DTs/abs(C(I,J)))
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
          ! update the solution for next time step
          DO I=2,NR
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(I,J)*FBND(I)+C(I-1,J)*FBND(I-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
1   CONTINUE
    DtsNext = min(DtsNext, DTDriftR)

    RETURN
  END

!************************************************************************
!                            DRIFTP
!               Calculate changes due to azimuthal drift
!************************************************************************
  SUBROUTINE DRIFTP

    use ModRamMain,      ONLY: Real8_, S
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts
    use ModRamVariables, ONLY: F2, FNIS, FNHS, BNES, VT, EIR, RLZ, MDR, DPHI

    implicit none
    save
    integer :: i, iSign, j, j1, k, l, n
    real(kind=Real8_) :: x, fup, r, corr, ome
    real(kind=Real8_) :: C(NR,NT,NE,NPA),VPA(NR,NT,NE,NPA), &
                         AGR(NR,NT,NE,NPA),GPA1,GPA2,GPA
    real(kind=Real8_) :: FBND(NT),F(NT),LIMITER

    DtDriftP = 100000.0
    OME=7.3E-5 ! Earth's angular velocity [rad/s]
    DO 1 L=2,NPA
      DO 1 K=2,NE
        DO 1 I=2,NR
          DO J=1,NT
            F(J)=F2(S,I,J,K,L)
          END DO
          DO J=2,NT
            J1=J+1
            IF (J.EQ.NT) J1=2
            GPA1=FNIS(I,J,L)+FNIS(I,J1,L)+(FNIS(I+1,J1,L)+FNIS(I+1,J,L)-FNIS(I-1,J,L)-FNIS(I-1,J1,L))*RLZ(I)/2./MDR
            GPA2=RLZ(I)/4./MDR*(FNIS(I,J,L)+FNIS(I,J1,L)-2*FNHS(I,J,L)-2*FNHS(I,J1,L)) &
                *(BNES(I+1,J1)+BNES(I+1,J)-BNES(I-1,J)-BNES(I-1,J1))/(BNES(I,J)+BNES(I,J1))
            GPA=GPA1+GPA2
            VPA(I,J,K,L)=P2(I,K)*GPA/(FNHS(I,J,L)+FNHS(I,J1,L))/(BNES(I,J)+BNES(I,J1))*DPHI/DTs*RLZ(I)
            C(I,J,K,L)=((VT(I+1,J)+VT(I+1,J1)-VT(I-1,J)-VT(I-1,J1))*P1(I)-P2(I,K)*GPA/(FNHS(I,J,L)+FNHS(I,J1,L)) &
                      -(EIR(I,J1)+EIR(I,J))/RLZ(I)*DTs/DPHI)/(BNES(I,J)+BNES(I,J1))+OME*DTs/DPHI
            AGR(I,J,K,L)=C(I,J,K,L)*DPHI/DTs
            DtDriftP = min(DtDriftP, FracCFL*DTs/abs(C(I,J,K,L)))
            ISIGN=1
            IF (C(I,J,K,L).NE.ABS(C(I,J,K,L))) ISIGN=-1
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
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(I,J,K,L)*FBND(J)+C(I,J-1,K,L)*FBND(J-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
          F2(S,I,1,K,L)=F2(S,I,NT,K,L)
1   CONTINUE
    DtsNext = min(DtsNext, DtDriftP)

    RETURN
  END

!**************************************************************************
!                       DRIFTE
!               Calculate energization along the drift path 
!**************************************************************************
  SUBROUTINE DRIFTE

    use ModRamMain,      ONLY: Real8_, S
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, dBdt, dIdt, EKEV, WE, RMAS, &
                               DPHI, RLZ, MDR, EBND, GREL, GRBND, DE, VT, EIR, EIP

    implicit none

    integer :: i, isign, j, j0, j2, k, l, n
    real(kind=Real8_) :: ezero,gpa,gpr1,gpr2,gpr3,gpp1,gpp2,edt1,qs, &
                         drdt, dpdt, dbdt1, didt1, x, fup, r, corr, ome
    real(kind=Real8_) :: GRZERO, DRD1(NR,NT),DPD1(NR,NT),GPR(NR,NT,NPA), &
                         GPP(NR,NT,NPA),DRD2(NR,NT,NPA),DPD2(NR,NT,NPA), &
                         EGR(NR,NT,NE,NPA)
    real(kind=Real8_) :: FBND(NE),F(0:NE+2),C(NE),LIMITER

    DtDriftE=10000.0
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
        DRD1(I,J)=(EIP(I,J)*RLZ(I)-(VT(I,J2)-VT(I,J0))/2./DPHI)/BNES(I,J)
        DPD1(I,J)=OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
        DO 1 L=2,NPA
          GPA=(1.-FNIS(I,J,L)/2./FNHS(I,J,L))/BNES(I,J)
          GPR1=GPA*(BNES(I+1,J)-BNES(I-1,J))/2./MDR
          GPR2=-FNIS(I,J,L)/FNHS(I,J,L)/RLZ(I)
          GPR3=-(FNIS(I+1,J,L)-FNIS(I-1,J,L))/2./MDR/FNHS(I,J,L)
          GPR(I,J,L)=GPR1+GPR2+GPR3
          GPP1=GPA*(BNES(I,J2)-BNES(I,J0))/2./DPHI
          GPP2=-(FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI/FNHS(I,J,L)
          GPP(I,J,L)=GPP1+GPP2
          DRD2(I,J,L)=(FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI &
                     +(FNIS(I,J,L)-2*FNHS(I,J,L))*(BNES(I,J2)-BNES(I,J0))/4/BNES(I,J)/DPHI
          DPD2(I,J,L)=FNIS(I,J,L)+(FNIS(I+1,J,L)-FNIS(I-1,J,L))*RLZ(I)/2/MDR &
                     +RLZ(I)*(FNIS(I,J,L)-2*FNHS(I,J,L))/4/MDR*(BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)
          DO K=2,NE
            F(K)=F2(S,I,J,K,L)
          END DO
          F(1)=F(2)*GREL(S,1)/GREL(S,2)*SQRT((GREL(S,2)**2-1)/(GREL(S,1)**2-1))
          F(0)=F(1)*GRZERO/GREL(S,1)*SQRT((GREL(S,1)**2-1)/(GRZERO**2-1))
          DO K=1,NE
            EDT1=EBND(K)*1e3*(GRBND(S,K)+1)/2/GRBND(S,K)/FNHS(I,J,L)/RLZ(I)/BNES(I,J)/QS
            DRDT=DRD1(I,J)+EDT1*DRD2(I,J,L)*RLZ(I)
            DPDT=DPD1(I,J)-EDT1*DPD2(I,J,L)
            dBdt1=dBdt(I,J)*(1.-FNIS(I,J,L)/2./FNHS(I,J,L))*RLZ(I)/BNES(I,J)
            dIdt1=-dIdt(I,J,L)*RLZ(I)/FNHS(I,J,L)
            C(K)=EDOT(I,K)*(GPR(I,J,L)*DRDT+GPP(I,J,L)*DPDT+dBdt1+dIdt1)
            EGR(I,J,K,L)=C(K)/DTs
            DtDriftE = min(DtDriftE, FracCFL*DTs*DE(K)/abs(C(K)))
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
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(K)/WE(K)*FBND(K)+C(K-1)/WE(K)*FBND(K-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
1   CONTINUE
    DtsNext = min(DtsNext, DtDriftE)

    RETURN
  END

!**************************************************************************
!                       DRIFTMU
!       Calculate pitch angle changes along the drift path
!**************************************************************************
  SUBROUTINE DRIFTMU

    use ModRamMain,      ONLY: Real8_, S
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts
    use ModRamVariables, ONLY: F2, BNES, BOUNIS, BOUNHS, FNHS, dBdt, dIbndt, &
                               RLZ, DPHI, MDR, GREL, EKEV, DMU, WMU, MU, &
                               VT, EIP, EIR

    implicit none

    integer :: i, j, j0, j1, k, l, n
    real(kind=Real8_) :: gmr1, gmr2, gmr3, gmp1, gmp2, qs, &
                         drdm, dpdm, dbdt2, dibndt2, x, fup, r, corr, ome
    real(kind=Real8_) :: CMUDOT(NR,NT,NPA),EDT(NPA),GMR(NR,NT,NPA),GMP(NR,NT,NPA),&
                         DRM2(NR,NT,NPA),DPM2(NR,NT,NPA),UGR(NR,NT,NE,NPA)
    real(kind=Real8_) :: FBND(NPA),F(NPA+1),C(NPA),LIMITER,DRM1(NR,NT), &
                         DPM1(NR,NT)
    integer :: URP(NPA),ISGM(NPA)

    DtDriftMu = 10000.0
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
          DRM1(I,J)=(EIP(I,J)*RLZ(I)-(VT(I,J1)-VT(I,J0))/2/DPHI)/BNES(I,J)
          DPM1(I,J)=OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
          C(1)=0.
          FBND(1)=0.
          UGR(I,J,K,1)=0.
          C(NPA)=0.
          FBND(NPA)=0.
          DO L=2,NPA
            F(L)=F2(S,I,J,K,L)
            CMUDOT(I,J,L)=MUDOT(I,L)*BOUNIS(I,J,L)/BOUNHS(I,J,L)
            GMR1=(BNES(I+1,J)-BNES(I-1,J))/4/MDR/BNES(I,J)
            GMR2=1/RLZ(I)
            GMR3=(BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))/2/MDR/BOUNIS(I,J,L)
            GMR(I,J,L)=GMR1+GMR2+GMR3
            GMP1=(BNES(I,J1)-BNES(I,J0))/4/DPHI/BNES(I,J)
            GMP2=(BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI/BOUNIS(I,J,L)
            GMP(I,J,L)=GMP1+GMP2
            EDT(L)=EKEV(K)*1e3*(GREL(S,K)+1)/2/GREL(S,K)/BOUNHS(I,J,L)/RLZ(I)/BNES(I,J)/QS
            DRM2(I,J,L)=(BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI+(BOUNIS(I,J,L)-2*BOUNHS(I,J,L)) &
                       *(BNES(I,J1)-BNES(I,J0))/4/BNES(I,J)/DPHI
            DPM2(I,J,L)=BOUNIS(I,J,L)+(BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))*RLZ(I)/2/MDR &
                       +(BOUNIS(I,J,L)-2*BOUNHS(I,J,L))*RLZ(I)/4/MDR*(BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)
            DRDM=DRM1(I,J)+EDT(L)*DRM2(I,J,L)*RLZ(I)
            DPDM=DPM1(I,J)-EDT(L)*DPM2(I,J,L)
            dBdt2=dBdt(I,J)/2./BNES(I,J)*RLZ(I)
            dIbndt2=dIbndt(I,J,L)*RLZ(I)/BOUNIS(I,J,L)
            C(L)=-CMUDOT(I,J,L)*(GMR(I,J,L)*DRDM+GMP(I,J,L)*DPDM+dBdt2+dIbndt2)
            UGR(I,J,K,L)=C(L)/DTs
            DtDriftMu=min(DtDriftMu,FracCFL*DTs*DMU(L)/max(1e-32,abs(C(L))))
            ISGM(L)=1
            IF(C(L).NE.ABS(C(L))) ISGM(L)=-1
            IF (ISGM(L).EQ.1) THEN
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
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-C(L)/WMU(L)*FBND(L)+C(L-1)/WMU(L)*FBND(L-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
          F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)*FNHS(I,J,NPA)*MU(NPA)/FNHS(I,J,NPA-1)/MU(NPA-1)
1    CONTINUE
    DtsNext = min(DtsNext, DtDriftMu)

    RETURN
  END

!==============================================================================
END MODULE ModRamDrift
