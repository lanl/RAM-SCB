MODULE ModRamDrift
! Contains subroutines responsible for calculating flux changes due to
! different drifts

  use ModRamVariables, ONLY: P1, P2, VR, EDOT, DtDriftR, DtDriftP, DtDriftE, &
                             DtDriftMu, FracCFL, MUDOT

  implicit none

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

    real(kind=Real8_) :: MUBOUN, QS!, RA(NS)
    integer :: i, j, k, l

    !DATA RA/1.,.77,.2,.03/    ! Proportions of species in the plasmasph
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
  END SUBROUTINE DriftPara

!==============================================================================
!************************************************************************
!                            DRIFTR
!               Calculate changes due to radial drift
!************************************************************************
  SUBROUTINE DRIFTR(ReCalculate)

    use ModRamMain,      ONLY: Real8_, S
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts, DtsNext, DtsMin
    use ModRamParams,    ONLY: BetaLim
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, MDR, EKEV, GREL, DPHI, &
                               RLZ, CONF1, CONF2, FGEOS, VT, EIP, CDriftR, &
                               sgnDriftR
    implicit none

    logical, intent(in) :: ReCalculate
    integer :: UR, i, j, j0, j1, k, l, n
    real(kind=Real8_) :: p4, qs, x, fup, r, corr, cgr1, cgr2, cgr3, ctemp
    real(kind=Real8_) :: CGR,CR(NR,NT),LIMITER
    real(kind=Real8_) :: F(NR+2),FBND(NR)

    If (ReCalculate) DTDriftR = 100000.0
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
        DO 1 J=1,NT
          F(1:NR) = F2(S,:,J,K,L)
          J0=J-1
          IF (J.EQ.1) J0=NT-1
          J1=J+1
          IF (J.EQ.NT) J1=2
          if (ReCalculate) then
            DO I=1,NR
              CGR1=FNIS(I+1,J1,L)+FNIS(I,J1,L)-FNIS(I+1,J0,L)-FNIS(I,J0,L)
              CGR2=BNES(I+1,J1)+BNES(I,J1)-BNES(I+1,J0)-BNES(I,J0)
              CGR3=CGR1+(FNIS(I+1,J,L)+FNIS(I,J,L)-2*FNHS(I+1,J,L) &
                  -2*FNHS(I,J,L))*CGR2/2./(BNES(I+1,J)+BNES(I,J))
              CGR=CGR3/(FNHS(I,J,L)+FNHS(I+1,J,L))*P4/2./(BNES(I,J)+BNES(I+1,J))/(RLZ(I)+0.5*MDR)
              CDriftR(I,J,K,L)=CR(I,J)+CGR
              ctemp = max(abs(CDriftR(I,J,K,L)),1E-10)
              DTDriftR = min( DTDriftR, FracCFL*DTs/ctemp)
              sgnDriftR(I,J,K,L)=1
              IF (CDriftR(I,J,K,L).NE.ABS(CDriftR(I,J,K,L))) sgnDriftR(I,J,K,L)=-1
            END DO
          endif
          IF (sgnDriftR(NR,J,K,L).EQ.1) THEN
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
            FUP=0.5*(F(I)+F(I+1)-sgnDriftR(I,J,K,L)*X)
            IF (ABS(X).LE.1.E-27) FBND(I)=FUP
            IF (ABS(X).GT.1.E-27) THEN
              N=I+1-sgnDriftR(I,J,K,L)
              R=(F(N)-F(N-1))/X
              IF (R.LE.0) FBND(I)=FUP
              IF (R.GT.0) THEN
                LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                CORR=-0.5*(CDriftR(I,J,K,L)-sgnDriftR(I,J,K,L))*X
                FBND(I)=FUP+LIMITER*CORR
              END IF
            END IF
          END DO
          ! update the solution for next time step
          DO I=2,NR
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftR(I,J,K,L)*FBND(I)+CDriftR(I-1,J,K,L)*FBND(I-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
1   CONTINUE
    DtsNext = min(DtsNext, max(DTDriftR,DtsMin))

    RETURN
  END

!************************************************************************
!                            DRIFTP
!               Calculate changes due to azimuthal drift
!************************************************************************
  SUBROUTINE DRIFTP(ReCalculate)

    use ModRamMain,      ONLY: Real8_, S
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, FNIS, FNHS, BNES, VT, EIR, RLZ, MDR, DPHI, CDriftP

    implicit none
    save

    logical, intent(in) :: ReCalculate
    integer :: i, sgn, j, j1, k, l, n
    real(kind=Real8_) :: x, fup, r, corr, ome, ctemp
    real(kind=Real8_) :: GPA1,GPA2
    real(kind=Real8_) :: FBND(NT),F(NT),LIMITER

    if (ReCalculate) DtDriftP = 100000.0
    OME=7.3E-5 ! Earth's angular velocity [rad/s]
    DO 1 L=2,NPA
      DO 1 K=2,NE
        DO 1 I=2,NR
          F(:)=F2(S,I,:,K,L)
          DO J=2,NT
            J1=J+1
            IF (J.EQ.NT) J1=2
            if (ReCalculate) then
              GPA1 = FNIS(I,J,L)+FNIS(I,J1,L) &
                    +(FNIS(I+1,J1,L)+FNIS(I+1,J,L)-FNIS(I-1,J,L)-FNIS(I-1,J1,L))*RLZ(I)/2./MDR
              GPA2 = RLZ(I)/4./MDR*(FNIS(I,J,L)+FNIS(I,J1,L)-2*FNHS(I,J,L)-2*FNHS(I,J1,L)) &
                    *(BNES(I+1,J1)+BNES(I+1,J)-BNES(I-1,J)-BNES(I-1,J1))/(BNES(I,J)+BNES(I,J1))
              CDriftP(I,J,K,L) = ((VT(I+1,J)+VT(I+1,J1)-VT(I-1,J)-VT(I-1,J1))*P1(I) &
                                  -P2(I,K)*(GPA1+GPA2)/(FNHS(I,J,L)+FNHS(I,J1,L)) &
                                  -(EIR(I,J1)+EIR(I,J))/RLZ(I)*DTs/DPHI)/(BNES(I,J) &
                                +BNES(I,J1))+OME*DTs/DPHI
              ctemp = max(abs(CDriftP(I,J,K,L)),1E-10)
              DtDriftP = min(DtDriftP, FracCFL*DTs/ctemp)
            endif
            sgn=1
            IF (CDriftP(I,J,K,L).NE.ABS(CDriftP(I,J,K,L))) sgn=-1
            X=F(J1)-F(J)
            FUP=0.5*(F(J)+F(J1)-sgn*X)
            IF (ABS(X).LE.1.E-27) FBND(J)=FUP
            IF (ABS(X).GT.1.E-27) THEN
              N=J+1-sgn
              IF (N.GT.NT) N=N-NT+1
              R=(F(N)-F(N-1))/X
              IF (R.LE.0) FBND(J)=FUP
              IF (R.GT.0) THEN
                LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                CORR=-0.5*(CDriftP(I,J,K,L)-sgn)*X
                FBND(J)=FUP+LIMITER*CORR
              END IF
            END IF
          END DO
          CDriftP(I,1,K,L)=CDriftP(I,NT,K,L)
          FBND(1)=FBND(NT)
          DO J=2,NT
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftP(I,J,K,L)*FBND(J)+CDriftP(I,J-1,K,L)*FBND(J-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
          F2(S,I,1,K,L)=F2(S,I,NT,K,L)
1   CONTINUE
    DtsNext = min(DtsNext, max(DtDriftP,DtsMin))

    RETURN
  END

!**************************************************************************
!                       DRIFTE
!               Calculate energization along the drift path 
!**************************************************************************
  SUBROUTINE DRIFTE(ReCalculate)

    use ModRamMain,      ONLY: Real8_, S
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, dBdt, dIdt, EKEV, WE, RMAS, &
                               DPHI, RLZ, MDR, EBND, GREL, GRBND, DE, VT, EIR, EIP, &
                               CDriftE

    implicit none

    logical, intent(in) :: ReCalculate
    integer :: i, sgn, j, j0, j2, k, l, n
    real(kind=Real8_) :: ezero,gpa,gpr1,gpr2,gpr3,gpp1,gpp2,edt1,qs, &
                         drdt, dpdt, dbdt1, didt1, x, fup, r, corr, ome
    real(kind=Real8_) :: GRZERO, DRD1,DPD1,DRD2,DPD2, ctemp
    real(kind=Real8_) :: FBND(NE),F(0:NE+2),LIMITER

    if (ReCalculate) DtDriftE=10000.0
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
        DRD1=(EIP(I,J)*RLZ(I)-(VT(I,J2)-VT(I,J0))/2./DPHI)/BNES(I,J)
        DPD1=OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
        DO 1 L=2,NPA
          if (ReCalculate) then
            GPA  = (1.-FNIS(I,J,L)/2./FNHS(I,J,L))/BNES(I,J)
            GPR1 = GPA*(BNES(I+1,J)-BNES(I-1,J))/2./MDR
            GPR2 = -FNIS(I,J,L)/FNHS(I,J,L)/RLZ(I)
            GPR3 = -(FNIS(I+1,J,L)-FNIS(I-1,J,L))/2./MDR/FNHS(I,J,L)
            GPP1 = GPA*(BNES(I,J2)-BNES(I,J0))/2./DPHI
            GPP2 = -(FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI/FNHS(I,J,L)
            DRD2 = (FNIS(I,J2,L)-FNIS(I,J0,L))/2./DPHI &
                  +(FNIS(I,J,L)-2*FNHS(I,J,L))*(BNES(I,J2)-BNES(I,J0))/4/BNES(I,J)/DPHI
            DPD2 = FNIS(I,J,L)+(FNIS(I+1,J,L)-FNIS(I-1,J,L))*RLZ(I)/2/MDR &
                  +RLZ(I)*(FNIS(I,J,L)-2*FNHS(I,J,L))/4/MDR*(BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)
          endif
          F(1:NE) = F2(S,I,J,:,L)
          F(1) = F(2)*GREL(S,1)/GREL(S,2)*SQRT((GREL(S,2)**2-1)/(GREL(S,1)**2-1))
          F(0) = F(1)*GRZERO/GREL(S,1)*SQRT((GREL(S,1)**2-1)/(GRZERO**2-1))
          DO K=1,NE
            if (ReCalculate) then
              EDT1  = EBND(K)*1e3*(GRBND(S,K)+1)/2/GRBND(S,K)/FNHS(I,J,L)/RLZ(I)/BNES(I,J)/QS
              DRDT  = DRD1+EDT1*DRD2*RLZ(I)
              DPDT  = DPD1-EDT1*DPD2
              dBdt1 = dBdt(I,J)*(1.-FNIS(I,J,L)/2./FNHS(I,J,L))*RLZ(I)/BNES(I,J)
              dIdt1 = -dIdt(I,J,L)*RLZ(I)/FNHS(I,J,L)
              CDriftE(I,J,K,L) = EDOT(I,K)*((GPR1+GPR2+GPR3)*DRDT+(GPP1+GPP2)*DPDT+dBdt1+dIdt1)
              ctemp = max(abs(CDriftE(I,J,K,L)),1E-10)
              DtDriftE = min(DtDriftE, FracCFL*DTs*DE(K)/ctemp)
            endif
            sgn=1
            IF(CDriftE(I,J,K,L).NE.ABS(CDriftE(I,J,K,L))) sgn=-1
            X=F(K+1)-F(K)
            FUP=0.5*(F(K)+F(K+1)-sgn*X)
            IF (ABS(X).LE.1.E-27) FBND(K)=FUP
            IF (ABS(X).GT.1.E-27) THEN
              N=K+1-sgn
              R=(F(N)-F(N-1))/X
              IF (R.LE.0) FBND(K)=FUP
              IF (R.GT.0) THEN
                LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                CORR=-0.5*(CDriftE(I,J,K,L)/DE(K)-sgn)*X
                FBND(K)=FUP+LIMITER*CORR
              END IF
            END IF
          END DO
          DO K=2,NE
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftE(I,J,K,L)/WE(K)*FBND(K)+CDriftE(I,J,K-1,L)/WE(K)*FBND(K-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
1   CONTINUE
    DtsNext = min(DtsNext, max(DtDriftE,DtsMin))

    RETURN
  END

!**************************************************************************
!                       DRIFTMU
!       Calculate pitch angle changes along the drift path
!**************************************************************************
  SUBROUTINE DRIFTMU(ReCalculate)

    use ModRamMain,      ONLY: Real8_, S
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, BNES, BOUNIS, BOUNHS, FNHS, dBdt, dIbndt, &
                               RLZ, DPHI, MDR, GREL, EKEV, DMU, WMU, MU, &
                               VT, EIP, EIR, CDriftMu

    implicit none

    logical, intent(in) :: ReCalculate
    integer :: i, j, j0, j1, k, l, n
    real(kind=Real8_) :: gmr1, gmr2, gmr3, gmp1, gmp2, qs, &
                         drdm, dpdm, dbdt2, dibndt2, x, fup, r, corr, ome
    real(kind=Real8_) :: CMUDOT,EDT,DRM2,DPM2,DRM1,DPM1, ctemp
    real(kind=Real8_) :: FBND(NPA),F(NPA),LIMITER
    integer :: ISGM

    if (ReCalculate) DtDriftMu = 10000.0
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
          F(:) = F2(S,I,J,K,:)
          F(1) = F(2)
          DRM1 = (EIP(I,J)*RLZ(I)-(VT(I,J1)-VT(I,J0))/2/DPHI)/BNES(I,J)
          DPM1 = OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
          DO L=2,NPA
            if (ReCalculate) then
              CMUDOT = MUDOT(I,L)*BOUNIS(I,J,L)/BOUNHS(I,J,L)
              GMR1 = (BNES(I+1,J)-BNES(I-1,J))/4/MDR/BNES(I,J)
              GMR2 = 1/RLZ(I)
              GMR3 = (BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))/2/MDR/BOUNIS(I,J,L)
              GMP1 = (BNES(I,J1)-BNES(I,J0))/4/DPHI/BNES(I,J)
              GMP2 = (BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI/BOUNIS(I,J,L)
              EDT  = EKEV(K)*1e3*(GREL(S,K)+1)/2/GREL(S,K)/BOUNHS(I,J,L)/RLZ(I)/BNES(I,J)/QS
              DRM2 = (BOUNIS(I,J1,L)-BOUNIS(I,J0,L))/2/DPHI+(BOUNIS(I,J,L)-2*BOUNHS(I,J,L)) &
                    *(BNES(I,J1)-BNES(I,J0))/4/BNES(I,J)/DPHI
              DPM2 = BOUNIS(I,J,L)+(BOUNIS(I+1,J,L)-BOUNIS(I-1,J,L))*RLZ(I)/2/MDR &
                    +(BOUNIS(I,J,L)-2*BOUNHS(I,J,L))*RLZ(I)/4/MDR*(BNES(I+1,J)-BNES(I-1,J))/BNES(I,J)
              DRDM = DRM1+EDT*DRM2*RLZ(I)
              DPDM = DPM1-EDT*DPM2
              dBdt2   = dBdt(I,J)/2./BNES(I,J)*RLZ(I)
              dIbndt2 = dIbndt(I,J,L)*RLZ(I)/BOUNIS(I,J,L)
              CDriftMu(I,J,K,L) = -CMUDOT*((GMR1+GMR2+GMR3)*DRDM+(GMP1+GMP2)*DPDM+dBdt2+dIbndt2)
              ctemp = max(abs(CDriftMu(I,J,K,L)),1E-32)
              DtDriftMu=min(DtDriftMu,FracCFL*DTs*DMU(L)/ctemp)
            endif
            ISGM=1
            IF(CDriftMu(I,J,K,L).NE.ABS(CDriftMu(I,J,K,L))) ISGM=-1
            if (L.LE.NPA-2) then
              X=F(L+1)-F(L)
              FUP=0.5*(F(L)+F(L+1)-ISGM*X)
              IF (ABS(X).LE.1.E-27) FBND(L)=FUP
              IF (ABS(X).GT.1.E-27) THEN
                N=L+1-ISGM
                R=(F(N)-F(N-1))/X
                IF (R.LE.0) FBND(L)=FUP
                IF (R.GT.0) THEN
                  LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                  CORR=-0.5*(CDriftMu(I,J,K,L)/DMU(L)-ISGM)*X
                  FBND(L)=FUP+LIMITER*CORR
                END IF
              END IF
            endif
          END DO
          CDriftMu(I,J,K,1)=0.
          FBND(1)     = 0.
          FBND(NPA-1) = F(NPA)
          DO L=2,NPA-1
            F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftMu(I,J,K,L)/WMU(L)*FBND(L)+CDriftMu(I,J,K,L-1)/WMU(L)*FBND(L-1)
            if (f2(s,i,j,k,l).lt.0) then
              f2(S,i,j,k,l)=1E-15
            endif
          END DO
          F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)*FNHS(I,J,NPA)*MU(NPA)/FNHS(I,J,NPA-1)/MU(NPA-1)
1    CONTINUE
    DtsNext = min(DtsNext, max(DtDriftMu,DtsMin))

    RETURN
  END

!==============================================================================
END MODULE ModRamDrift
