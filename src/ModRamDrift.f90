!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamDrift
! Contains subroutines responsible for calculating flux changes due to
! different drifts

  use ModRamMain,      ONLY: Real8_
  use ModRamGrids,     ONLY: nR, nE, nPA
  use ModRamVariables, ONLY: FracCFL, DtDriftR, DtDriftP, DtDriftE, DtDriftMu

  implicit none
  save

  integer :: QS
  real(kind=Real8_), ALLOCATABLE :: VR(:), P1(:), P2(:,:), MUDOT(:,:), EDOT(:,:), &
                                    CDriftR(:,:,:,:), CDriftP(:,:,:,:), CDriftE(:,:,:,:), &
                                    CDriftMu(:,:,:,:)
  !$OMP THREADPRIVATE(QS, VR, P1, P2, MUDOT, EDOT, CDriftR, CDriftP, CDriftE, CDriftMu)

contains
!==============================================================================
  SUBROUTINE DRIFTEND

    implicit none
    
    DEALLOCATE(VR, P1, P2, EDOT, MUDOT, CDriftR, CDriftP, CDriftE, CDriftMu)

  END SUBROUTINE DRIFTEND

!**************************************************************************
!                               OTHERPARA
!                       Calculate drift parameters
!**************************************************************************
  SUBROUTINE DRIFTPARA(S)

    use ModRamMain,      ONLY: Real8_
    use ModRamConst,     ONLY: PI
    use ModRamGrids,     ONLY: NS, NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: RLZ, MDR, DPHI, EKEV, GREL, WMU, EBND, GRBND, &
                               PHIOFS, MU

    implicit none
    save

    integer, intent(in) :: S
    real(kind=Real8_) :: MUBOUN!, RA(NS)
    integer :: i, j, k, l
    !$OMP THREADPRIVATE(MUBOUN)

    ALLOCATE(VR(nR), P1(nR), P2(nR,nE), EDOT(nR,nE), MUDOT(nR,nPA))
    ALLOCATE(CDriftR(nR,nT,nE,nPa), CDriftP(nR,nT,nE,nPa), &
             CDriftE(nR,nT,nE,nPa), CDriftMu(nR,nT,nE,nPa))
    
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
  SUBROUTINE DRIFTR(S)

    use ModRamMain,      ONLY: Real8_
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts, DtsNext, DtsMin
    use ModRamParams,    ONLY: BetaLim
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, MDR, EKEV, GREL, DPHI, &
                               RLZ, CONF1, CONF2, FGEOS, VT, EIP
    implicit none
    save

    integer, intent(in) :: S
    integer :: UR, i, j, j0, j1, k, l, n, sgn
    real(kind=Real8_) :: p4, x, fup, r, corr, cgr1, cgr2, cgr3, ctemp
    real(kind=Real8_) :: CGR,CR(NR,NT),LIMITER, DtTemp
    real(kind=Real8_) :: F(NR+2),FBND(NR)
    !$OMP THREADPRIVATE(UR, J0, J1, P4, X, FUP, R, CORR, CGR1, CGR2, CGR3, CTEMP)
    !$OMP THREADPRIVATE(CGR, CR, LIMITER, F, FBND, sgn)

    DTDriftR(S) = 100000.0
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

    DO K=2,NE
       P4=DTs*EKEV(K)*1000.0*(GREL(S,K)+1)/GREL(S,K)/DPHI/MDR/QS
       DO L=2,NPA
          DO J=1,NT
             F(1:NR) = F2(S,:,J,K,L)
             J0=J-1
             IF (J.EQ.1) J0=NT-1
             J1=J+1
             IF (J.EQ.NT) J1=2
             DO I=1,NR
                CGR1=FNIS(I+1,J1,L)+FNIS(I,J1,L)-FNIS(I+1,J0,L)-FNIS(I,J0,L)
                CGR2=BNES(I+1,J1)+BNES(I,J1)-BNES(I+1,J0)-BNES(I,J0)
                CGR3=CGR1+(FNIS(I+1,J,L)+FNIS(I,J,L)-2*FNHS(I+1,J,L) &
                     -2*FNHS(I,J,L))*CGR2/2./(BNES(I+1,J)+BNES(I,J))
                CGR=CGR3/(FNHS(I,J,L)+FNHS(I+1,J,L))*P4/2./(BNES(I,J)+BNES(I+1,J))/(RLZ(I)+0.5*MDR)
                CDriftR(I,J,K,L)=CR(I,J)+CGR
                ctemp = max(abs(CDriftR(I,J,K,L)),1E-10)
                DTDriftR(S) = min( DTDriftR(S), FracCFL*DTs/ctemp)
             END DO
             sgn = 1
             IF (CDriftR(nR,J,K,L).NE.ABS(CDriftR(nR,J,K,L))) sgn=-1
             IF (sgn.EQ.1) THEN
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
                FUP=0.5*(F(I)+F(I+1)-sgn*X)
                IF (ABS(X).LE.1.E-27) FBND(I)=FUP
                IF (ABS(X).GT.1.E-27) THEN
                   N=I+1-sgn
                   R=(F(N)-F(N-1))/X
                   IF (R.LE.0) FBND(I)=FUP
                   IF (R.GT.0) THEN
                      LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                      CORR=-0.5*(CDriftR(I,J,K,L)-sgn)*X
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
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE DriftR

!************************************************************************
!                            DRIFTP
!               Calculate changes due to azimuthal drift
!************************************************************************
  SUBROUTINE DRIFTP(S)

    use ModRamMain,      ONLY: Real8_
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, FNIS, FNHS, BNES, VT, EIR, RLZ, MDR, DPHI

    implicit none
    save

    integer, intent(in) :: S
    integer :: i, sgn, j, j1, k, l, n
    real(kind=Real8_) :: x, fup, r, corr, ome, ctemp
    real(kind=Real8_) :: GPA1,GPA2
    real(kind=Real8_) :: FBND(NT),F(NT),LIMITER
    !$OMP THREADPRIVATE(SGN, J1, X, FUP, R, CORR, OME, CTEMP, GPA1, GPA2, FBND, F, LIMITER)

    DtDriftP(S) = 100000.0
    OME=7.3E-5 ! Earth's angular velocity [rad/s]
    DO L=2,NPA
       DO K=2,NE
          DO I=2,NR
             F(:)=F2(S,I,:,K,L)
             DO J=2,NT
                J1=J+1
                IF (J.EQ.NT) J1=2
                GPA1 = FNIS(I,J,L)+FNIS(I,J1,L)+(FNIS(I+1,J1,L)+FNIS(I+1,J,L)-FNIS(I-1,J,L) &
                      -FNIS(I-1,J1,L))*RLZ(I)/2./MDR
                GPA2 = RLZ(I)/4./MDR*(FNIS(I,J,L)+FNIS(I,J1,L)-2*FNHS(I,J,L)-2*FNHS(I,J1,L)) &
                      *(BNES(I+1,J1)+BNES(I+1,J)-BNES(I-1,J)-BNES(I-1,J1))/(BNES(I,J)+BNES(I,J1))
                CDriftP(I,J,K,L) = ((VT(I+1,J)+VT(I+1,J1)-VT(I-1,J)-VT(I-1,J1))*P1(I) &
                                  -P2(I,K)*(GPA1+GPA2)/(FNHS(I,J,L)+FNHS(I,J1,L)) &
                                  -(EIR(I,J1)+EIR(I,J))/RLZ(I)*DTs/DPHI)/(BNES(I,J) &
                                  +BNES(I,J1))+OME*DTs/DPHI
                ctemp = max(abs(CDriftP(I,J,K,L)),1E-10)
                DtDriftP(S) = min(DtDriftP(S), FracCFL*DTs/ctemp)
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
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE DriftP

!**************************************************************************
!                       DRIFTE
!               Calculate energization along the drift path 
!**************************************************************************
  SUBROUTINE DRIFTE(S)

    use ModRamMain,      ONLY: Real8_
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA, nS
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, dBdt, dIdt, EKEV, WE, RMAS, &
                               DPHI, RLZ, MDR, EBND, GREL, GRBND, DE, VT, EIR, EIP

    implicit none
    save

    integer, intent(in) :: S
    integer :: i, sgn, j, j0, j2, k, l, n
    real(kind=Real8_) :: ezero,gpa,gpr1,gpr2,gpr3,gpp1,gpp2,edt1, &
                         drdt, dpdt, dbdt1, didt1, x, fup, r, corr, ome
    real(kind=Real8_) :: GRZERO(nS), DRD1,DPD1,DRD2,DPD2, ctemp
    real(kind=Real8_) :: FBND(NE),F(0:NE+2),LIMITER
    !$OMP THREADPRIVATE(SGN, J0, J2, GPA, GPR1, GPR2, GPR3, GPP1, GPP2, EDT1, DRDT)
    !$OMP THREADPRIVATE(DPDT, DBDT1, DIDT1, X, FUP, R, CORR, DRD1, DPD1, DRD2, DPD2)
    !$OMP THREADPRIVATE(CTEMP, FBND, F, LIMITER)

    DtDriftE(S)=10000.0
    QS=1.
    IF (S.EQ.1) QS=-1.
    OME=7.3E-5
    EZERO=EKEV(1)-WE(1)
    GRZERO(S)=1.+EZERO*1000.*Q/RMAS(S)/CS/CS
    F(NE+1)=0.
    F(NE+2)=0.
    DO J=1,NT
       J0=J-1
       IF (J.EQ.1) J0=NT-1
       J2=J+1
       IF (J.EQ.NT) J2=2
       DO I=2,NR
          DRD1=(EIP(I,J)*RLZ(I)-(VT(I,J2)-VT(I,J0))/2./DPHI)/BNES(I,J)
          DPD1=OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
          DO L=2,NPA
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
             F(1:NE) = F2(S,I,J,:,L)
             F(1) = F(2)*GREL(S,1)/GREL(S,2)*SQRT((GREL(S,2)**2-1)/(GREL(S,1)**2-1))
             F(0) = F(1)*GRZERO(S)/GREL(S,1)*SQRT((GREL(S,1)**2-1)/(GRZERO(S)**2-1))
             DO K=1,NE
                EDT1  = EBND(K)*1e3*(GRBND(S,K)+1)/2/GRBND(S,K)/FNHS(I,J,L)/RLZ(I)/BNES(I,J)/QS
                DRDT  = DRD1+EDT1*DRD2*RLZ(I)
                DPDT  = DPD1-EDT1*DPD2
                dBdt1 = dBdt(I,J)*(1.-FNIS(I,J,L)/2./FNHS(I,J,L))*RLZ(I)/BNES(I,J)
                dIdt1 = -dIdt(I,J,L)*RLZ(I)/FNHS(I,J,L)
                CDriftE(I,J,K,L) = EDOT(I,K)*((GPR1+GPR2+GPR3)*DRDT+(GPP1+GPP2)*DPDT+dBdt1+dIdt1)
                ctemp = max(abs(CDriftE(I,J,K,L)),1E-10)
                DtDriftE(S) = min(DtDriftE(S), FracCFL*DTs*DE(K)/ctemp)
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
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE DriftE

!**************************************************************************
!                       DRIFTMU
!       Calculate pitch angle changes along the drift path
!**************************************************************************
  SUBROUTINE DRIFTMU(S)

    use ModRamMain,      ONLY: Real8_
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: DtsNext, Dts, DtsMin
    use ModRamVariables, ONLY: F2, BNES, BOUNIS, BOUNHS, FNHS, dBdt, dIbndt, &
                               RLZ, DPHI, MDR, GREL, EKEV, DMU, WMU, MU, &
                               VT, EIP, EIR

    implicit none
    save

    integer, intent(in) :: S
    integer :: i, j, j0, j1, k, l, n
    real(kind=Real8_) :: gmr1, gmr2, gmr3, gmp1, gmp2, &
                         drdm, dpdm, dbdt2, dibndt2, x, fup, r, corr, ome
    real(kind=Real8_) :: CMUDOT,EDT,DRM2,DPM2,DRM1,DPM1, ctemp
    real(kind=Real8_) :: FBND(NPA),F(NPA),LIMITER
    integer :: ISGM
    !$OMP THREADPRIVATE(J0, J1, GMR1, GMR2, GMR3, GMP1, GMP2, DRDM, DPDM, DBDT2, DIBNDT2)
    !$OMP THREADPRIVATE(X, FUP, R, CORR, CMUDOT, EDT, DRM2, DPM2, DRM1, DPM1, CTEMP)
    !$OMP THREADPRIVATE(FBND, F, LIMITER, ISGM)

    DtDriftMu(S) = 10000.0
    QS=1.
    IF (S.EQ.1) QS=-1.
    OME=7.3E-5
    DO K=2,NE
       DO J=1,NT
          J0=J-1
          IF (J.EQ.1) J0=NT-1
          J1=J+1
          IF (J.EQ.NT) J1=2
          DO I=2,NR
             F(:) = F2(S,I,J,K,:)
             F(1) = F(2)
             DRM1 = (EIP(I,J)*RLZ(I)-(VT(I,J1)-VT(I,J0))/2/DPHI)/BNES(I,J)
             DPM1 = OME*RLZ(I)+((VT(I+1,J)-VT(I-1,J))/2/MDR-EIR(I,J))/BNES(I,J)
             DO L=2,NPA
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
                DtDriftMu(S)=min(DtDriftMu(S),FracCFL*DTs*DMU(L)/ctemp)
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
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE DriftMu

!!==============================================================================
! SUBROUTINE DriftCalculations(order) 
!! An attempt to merge all four drift subroutines into one, turns out to be
!! slower since large arrays are needed to store C and FBND variables, but
!! we might get a large speed up if we MPI the species and OpenMP the drifts
!
!    use ModRamMain,      ONLY: Real8_, S
!    use ModRamParams,    ONLY: betalim
!    use ModRamConst,     ONLY: Q, CS
!    use ModRamTiming,    ONLY: Dts, DtsNext, UTs
!    use ModRamGrids,     ONLY: NR, NE, NT, NPA
!    use ModRamVariables, ONLY: BNES, FNIS, FNHS, BOUNHS, BOUNIS, VT, EIP, EIR, &
!                               RLZ, EKEV, WE, RMAS, DPHI, MDR, DE, dIdt, dBdt, &
!                               dIbndt, P1, P2, GREL, CONF1, CONF2, FGEOS, F2,  &
!                               DE, DMU, WMU, MU, GRBND, EBND
!    implicit none
!
!    integer, intent(in) :: order
!    integer :: Io, Im, Ip, Jo, Jm, Jp, K, L, sgn, N
!    real(kind=Real8_) :: QS, OME, EZERO, GRZERO
!    real(kind=Real8_) :: X, R, FUP, CORR, LIMITER
!    real(kind=Real8_) :: CDriftR(NR,NT,NE,NPA), CDriftP(NR,NT,NE,NPA), &
!                         CDriftE(NR,NT,NE,NPA), CDriftMu(NR,NT,NE,NPA)
!
!    real(kind=Real8_) :: CR, P4, CGR1, CGR2, CGR3, FRp1, FRp2, FBNDR(NR,NT,NE,NPA)
!
!    real(kind=Real8_) :: GPA1, GPA2, GPA3, FBNDP(NR,NT,NE,NPA)
!
!    real(kind=Real8_) :: gpa, gpr1, gpr2, gpr3, gpp1, gpp2, edt1, drdt, dpdt,  &
!                         dbdt1, didt1, DRD1, DPD1, GPR, GPP, DRD2, DPD2, FEm1, &
!                         FEp1, FEp2, FBNDE(NR,NT,NE,NPA)
!
!    real(kind=Real8_) :: gmr1, gmr2, gmr3, gmp1, gmp2, drm1, dpm1, drdm, dpdm, &
!                         dbdt2, dibndt2, CMUDOT, EDT, GMR, GMP, DRM2, DPM2, FMm1, &
!                         FBNDMU(NR,NT,NE,NPA)
!
!DtDriftR  = 10000.0
!DtDriftP  = 10000.0
!DtDriftE  = 10000.0
!DtDriftMu = 10000.0
!QS = 1.
!IF (S.EQ.1) QS = -1.
!OME = 7.3E-5
!EZERO = EKEV(1) - WE(1)
!GRZERO = 1. + EZERO*1000.*Q/RMAS(S)/CS/CS
!Do Io=1,NR
!  Im = Io-1
!  IF (Io.EQ.1) Im = 1
!  Ip = Io+1
!!  IF (Io.EQ.NR) Ip = 2
!  DO Jo=1,NT
!    Jm = Jo-1
!    IF (Jo.EQ.1) Jm = NT-1
!    Jp = Jo+1
!    if (Jo.EQ.NT) Jp = 2
!! DriftR
!    CR=VR(Io)*(VT(Io,Jm)+VT(Ip,Jm)-VT(Io,Jp)-VT(Ip,Jp))/(BNES(Io,Jo)+BNES(Ip,Jo)) &
!      +(EIP(Io,Jo)+EIP(Ip,Jo))/(BNES(Io,Jo)+BNES(Ip,Jo))*DTs/MDR
!! DriftE
!    DRD1=(EIP(Io,Jo)*RLZ(Io)-(VT(Io,Jp)-VT(Io,Jm))/2./DPHI)/BNES(Io,Jo)
!    DPD1=OME*RLZ(Io)+((VT(Ip,Jo)-VT(Im,Jo))/2/MDR-EIR(Io,Jo))/BNES(Io,Jo)
!! DriftMu (actually the same values as the DriftE values)
!    DRM1 = DRD1
!    DPM1 = DPD1
!    DO L=2,NPA
!! DriftR
!      CGR1=FNIS(Ip,Jp,L)+FNIS(Io,Jp,L)-FNIS(Ip,Jm,L)-FNIS(Io,Jm,L)
!      CGR2=BNES(Ip,Jp)+BNES(Io,Jp)-BNES(Ip,Jm)-BNES(Io,Jm)
!      CGR3=CGR1+(FNIS(Ip,Jo,L)+FNIS(Io,Jo,L)-2.*FNHS(Ip,Jo,L)-2.*FNHS(Io,Jo,L))*CGR2/2./(BNES(Ip,Jo)+BNES(Io,Jo))
!! DriftP
!      GPA1=FNIS(Io,Jo,L)+FNIS(Io,Jp,L)+(FNIS(Ip,Jp,L)+FNIS(Ip,Jo,L)-FNIS(Im,Jo,L)-FNIS(Im,Jp,L))*RLZ(Io)/2./MDR
!      GPA2=RLZ(Io)/4./MDR*(FNIS(Io,Jo,L)+FNIS(Io,Jp,L)-2.*FNHS(Io,Jo,L)-2.*FNHS(Io,Jp,L)) &
!          *(BNES(Ip,Jp)+BNES(Ip,Jo)-BNES(Im,Jo)-BNES(Im,Jp))/(BNES(Io,Jo)+BNES(Io,Jp))
!      GPA3=GPA1+GPA2
!! DriftE
!      GPA=(1.-FNIS(Io,Jo,L)/2./FNHS(Io,Jo,L))/BNES(Io,Jo)
!      GPR1=GPA*(BNES(Ip,Jo)-BNES(Im,Jo))/2./MDR
!      GPR2=-FNIS(Io,Jo,L)/FNHS(Io,Jo,L)/RLZ(Io)
!      GPR3=-(FNIS(Ip,Jo,L)-FNIS(Im,Jo,L))/2./MDR/FNHS(Io,Jo,L)
!      GPR=GPR1+GPR2+GPR3
!      GPP1=GPA*(BNES(Io,Jp)-BNES(Io,Jm))/2./DPHI
!      GPP2=-(FNIS(Io,Jp,L)-FNIS(Io,Jm,L))/2./DPHI/FNHS(Io,Jo,L)
!      GPP=GPP1+GPP2
!      DRD2=(FNIS(Io,Jp,L)-FNIS(Io,Jm,L))/2./DPHI+(FNIS(Io,Jo,L)-2.*FNHS(Io,Jo,L))
!&
!          *(BNES(Io,Jp)-BNES(Io,Jm))/4./BNES(Io,Jo)/DPHI
!      DPD2=FNIS(Io,Jo,L)+(FNIS(Ip,Jo,L)-FNIS(Im,Jo,L))*RLZ(Io)/2./MDR+RLZ(Io) &
!          *(FNIS(Io,Jo,L)-2.*FNHS(Io,Jo,L))/4./MDR*(BNES(Ip,Jo)-BNES(Im,Jo))/BNES(Io,Jo)
!      dBdt1=dBdt(Io,Jo)*(1.-FNIS(Io,Jo,L)/2./FNHS(Io,Jo,L))*RLZ(Io)/BNES(Io,Jo)
!      dIdt1=-dIdt(Io,Jo,L)*RLZ(Io)/FNHS(Io,Jo,L)
!! DriftMu
!      CMUDOT=MUDOT(Io,L)*BOUNIS(Io,Jo,L)/BOUNHS(Io,Jo,L)
!      GMR1=(BNES(Ip,Jo)-BNES(Im,Jo))/4./MDR/BNES(Io,Jo)
!      GMR2=1/RLZ(Io)
!      GMR3=(BOUNIS(Ip,Jo,L)-BOUNIS(Im,Jo,L))/2./MDR/BOUNIS(Io,Jo,L)
!      GMR=GMR1+GMR2+GMR3
!      GMP1=(BNES(Io,Jp)-BNES(Io,Jm))/4./DPHI/BNES(Io,Jo)
!      GMP2=(BOUNIS(Io,Jp,L)-BOUNIS(Io,Jm,L))/2./DPHI/BOUNIS(Io,Jo,L)
!      GMP=GMP1+GMP2
!      DRM2=(BOUNIS(Io,Jp,L)-BOUNIS(Io,Jm,L))/2./DPHI+(BOUNIS(Io,Jo,L)-2.*BOUNHS(Io,Jo,L)) &
!          *(BNES(Io,Jp)-BNES(Io,Jm))/4./BNES(Io,Jo)/DPHI
!      DPM2=BOUNIS(Io,Jo,L)+(BOUNIS(Ip,Jo,L)-BOUNIS(Im,Jo,L))*RLZ(Io)/2./MDR &
!          +(BOUNIS(Io,Jo,L)-2.*BOUNHS(Io,Jo,L))*RLZ(Io)/4./MDR*(BNES(Ip,Jo)-BNES(Im,Jo))/BNES(Io,Jo)
!      dBdt2=dBdt(Io,Jo)/2./BNES(Io,Jo)*RLZ(Io)
!      dIbndt2=dIbndt(Io,Jo,L)*RLZ(Io)/BOUNIS(Io,Jo,L)
!      DO K=1,NE
!!!!!!!! Calculate Drift Stuff
!! DriftR
!        IF (K.GT.1) THEN
!          P4=DTs*EKEV(K)*1000.0*(GREL(S,K)+1)/GREL(S,K)/DPHI/MDR/QS
!          CDriftR(Io,Jo,K,L)=CR + CGR3/(FNHS(Io,Jo,L)+FNHS(Ip,Jo,L))*P4/2./(BNES(Io,Jo)+BNES(Ip,Jo))/(RLZ(Io)+0.5*MDR)
!          DTDriftR = min( DTDriftR, FracCFL*DTs/abs(CDriftR(Io,Jo,K,L)))
!        ENDIF
!! DriftP
!        IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.GT.1)) THEN
!          CDriftP(Io,Jo,K,L)=((VT(Ip,Jo)+VT(Ip,Jp)-VT(Im,Jo)-VT(Im,Jp))*P1(Io)-P2(Io,K)*GPA3/(FNHS(Io,Jo,L)+FNHS(Io,Jp,L)) &
!                            -(EIR(Io,Jp)+EIR(Io,Jo))/RLZ(Io)*DTs/DPHI)/(BNES(Io,Jo)+BNES(Io,Jp))+OME*DTs/DPHI
!          DtDriftP = min(DtDriftP, FracCFL*DTs/abs(CDriftP(Io,Jo,K,L)))
!        ENDIF
!! DriftE
!        IF (Io.GT.1) THEN
!          EDT1=EBND(K)*1e3*(GRBND(S,K)+1)/2/GRBND(S,K)/FNHS(Io,Jo,L)/RLZ(Io)/BNES(Io,Jo)/QS
!          DRDT=DRD1+EDT1*DRD2*RLZ(Io)
!          DPDT=DPD1-EDT1*DPD2
!          CDriftE(Io,Jo,K,L)=EDOT(Io,K)*(GPR*DRDT+GPP*DPDT+dBdt1+dIdt1)
!          DtDriftE = min(DtDriftE, FracCFL*DTs*DE(K)/abs(CDriftE(Io,Jo,K,L)))
!        ENDIF
!! DriftMu
!        IF ((K.GT.1).AND.(Io.GT.1)) THEN
!          EDT=EKEV(K)*1e3*(GREL(S,K)+1)/2/GREL(S,K)/BOUNHS(Io,Jo,L)/RLZ(Io)/BNES(Io,Jo)/QS
!          DRDM=DRM1+EDT*DRM2*RLZ(Io)
!          DPDM=DPM1-EDT*DPM2
!          CDriftMu(Io,Jo,K,L)=-CMUDOT*(GMR*DRDM+GMP*DPDM+dBdt2+dIbndt2)
!          DtDriftMu=min(DtDriftMu,FracCFL*DTs*DMU(L)/max(1e-32,abs(CDriftMu(Io,Jo,K,L))))
!        ENDIF
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO
!
!DO Io = 2,NR
!  Im = Io-1
!  IF (Io.EQ.1) Im = 1
!  Ip = Io+1
!!  IF (Io.EQ.NR) Ip = 2
!  DO Jo=2,NT
!    Jm = Jo-1
!    IF (Jo.EQ.1) Jm = NT-1
!    Jp = Jo+1
!    if (Jo.EQ.NT) Jp = 2
!    DO K = 1,NE
!      FBNDMu(Io,Jo,K,1) = 0.
!      FBNDMu(Io,Jo,K,NPA-1) = F2(S,Io,Jo,K,NPA)
!      CDriftMu(Io,Jo,K,1) = 0.
!      DO L = 2,NPA
!!!!!!!! Calculate Drift Stuff
!! DriftR
!        IF (K.GT.1) THEN
!          IF (Io.GT.1) THEN
!            if (Io.eq.2) FBNDR(1,Jo,K,L) = F2(S,2,Jo,K,L)
!            sgn = 1
!            if (CDriftR(NR,Jo,K,L).NE.ABS(CDriftR(NR,Jo,K,L))) sgn = -1
!            N = Io+1-sgn
!            IF (Io.EQ.NR) THEN
!              X = FGEOS(S,Jo,K,L)*CONF1*FNHS(NR,Jo,L)-F2(S,Io,Jo,K,L)
!              FUP = 0.5*(FGEOS(S,Jo,K,L)*CONF1*FNHS(NR,Jo,L)+F2(S,Io,Jo,K,L)-sgn*X)
!            ELSE
!              X   = F2(S,Ip,Jo,K,L)-F2(S,Io,Jo,K,L)
!              FUP = 0.5*(F2(S,Ip,Jo,K,L)+F2(S,Io,Jo,K,L)-sgn*X)
!            ENDIF
!            IF (ABS(X).LE.1.E-27) THEN
!              FBNDR(Io,Jo,K,L) = FUP
!            ELSE
!              IF (N.EQ.NR+2) then
!                R = (FGEOS(S,Jo,K,L)*CONF2*FNHS(NR,Jo,L)-FGEOS(S,Jo,K,L)*CONF1*FNHS(NR,Jo,L))/X
!              ELSEIF (N.EQ.NR+1) then
!                R = (FGEOS(S,Jo,K,L)*CONF1*FNHS(NR,Jo,L)-F2(S,NR,Jo,K,L))/X
!              ELSE
!                R = (F2(S,N,Jo,K,L)-F2(S,N-1,Jo,K,L))/X
!              ENDIF
!              if (R.le.0) FBNDR(Io,Jo,K,L) = FUP
!              If (R.gt.0) then
!                LIMITER = MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
!                CORR = -0.5*(CDriftR(Io,Jo,K,L)-sgn)*X
!                FBNDR(Io,Jo,K,L) = FUP+LIMITER*CORR
!              endif
!            ENDIF
!            if ((sgn.eq.1).and.(Io.eq.NR)) FBNDR(NR,Jo,K,L) = F2(S,NR,Jo,K,L)
!            if ((sgn.eq.1).and.(Io.eq.2))  FBNDR(1,Jo,K,L)  = 0.
!          ENDIF
!        ENDIF
!! DriftP
!        IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.GT.1)) THEN
!          sgn = 1
!          if (CDriftP(Io,Jo,K,L).NE.ABS(CDriftP(Io,Jo,K,L))) sgn = -1
!          X   = F2(S,Io,Jp,K,L)-F2(S,Io,Jo,K,L)
!          FUP = 0.5*(F2(S,Io,Jp,K,L)+F2(S,Io,Jo,K,L)-sgn*X)
!          IF (ABS(X).LE.1.E-27) THEN
!            FBNDP(Io,Jo,K,L) = FUP
!          ELSE
!            N   = Jo+1-sgn
!            if (N.GT.NT) N=N-NT+1
!            R = (F2(S,Io,N,K,L)-F2(S,Io,N-1,K,L))/X
!            IF (R.LE.0) FBNDP(Io,Jo,K,L)=FUP
!            IF (R.GT.0) THEN
!              LIMITER = MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
!              CORR    = -0.5*(CDriftP(Io,Jo,K,L)-sgn)*X
!              FBNDP(Io,Jo,K,L)=FUP+LIMITER*CORR
!            ENDIF
!          ENDIF
!        ENDIF
!        IF (Jo.EQ.NT) THEN
!          CDriftP(Io,1,K,L) = CDriftP(Io,NT,K,L)
!          FBNDP(Io,1,K,L)   = FBNDP(Io,NT,K,L)
!        ENDIF
!! DriftE
!        IF (Io.GT.1) THEN
!          sgn = 1
!          if (CDriftE(Io,Jo,K,L).NE.ABS(CDriftE(Io,Jo,K,L))) sgn = -1
!          N = K+1-sgn
!          IF (K.EQ.NE) THEN
!            X   = 0.-F2(S,Io,Jo,K,L)
!            FUP = 0.5*(0.+F2(S,Io,Jo,K,L)-sgn*X)
!          ELSE
!            X   = F2(S,Io,Jo,K+1,L)-F2(S,Io,Jo,K,L)
!            FUP = 0.5*(F2(S,Io,Jo,K+1,L)+F2(S,Io,Jo,K,L)-sgn*X)
!          ENDIF
!          IF (ABS(X).LE.1.E-27) THEN
!            FBNDE(Io,Jo,K,L) = FUP
!          ELSE
!            IF (N.EQ.NE+2) THEN
!              R = (0.-0.)/X
!            ELSEIF (N.EQ.NE+1) THEN
!              R = (0.-F2(S,Io,Jo,N-1,L))/X
!            ELSEIF (N.EQ.2) THEN
!              R = (F2(S,Io,Jo,2,L) &
!                 - F2(S,Io,Jo,2,L)*GREL(S,1)/GREL(S,2) &
!                 *SQRT((GREL(S,2)**2-1)/(GREL(S,1)**2-1)))/X
!            ELSEIF (N.EQ.1) THEN
!              R = (F2(S,Io,Jo,2,L)*GREL(S,1)/GREL(S,2)*SQRT((GREL(S,2)**2-1)/(GREL(S,1)**2-1)) &
!                 - F2(S,Io,Jo,1,L)*GRZERO/GREL(S,1)*SQRT((GREL(S,1)**2-1)/(GRZERO**2-1)))/X
!            ELSE
!              R = (F2(S,Io,Jo,N,L)-F2(S,Io,Jo,N-1,L))/X
!            ENDIF
!            IF (R.LE.0) FBNDE(Io,Jo,K,L)=FUP
!            IF (R.GT.0) THEN
!              LIMITER = MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
!              CORR    = -0.5*(CDriftE(Io,Jo,K,L)-sgn)*X
!              FBNDE(Io,Jo,K,L)=FUP+LIMITER*CORR
!            ENDIF
!          ENDIF
!        ENDIF
!! DriftMu
!        IF ((K.GT.1).AND.(Io.GT.1)) THEN
!          IF (L.LE.NPA-2) THEN
!            sgn = 1
!            if (CDriftMu(Io,Jo,K,L).NE.ABS(CDriftMu(Io,Jo,K,L))) sgn = -1
!            X   = F2(S,Io,Jo,K,L+1)-F2(S,Io,Jo,K,L)
!            FUP = 0.5*(F2(S,Io,Jo,K,L+1)+F2(S,Io,Jo,K,L)-sgn*X)
!            IF (ABS(X).LE.1.E-27) THEN
!              FBNDMu(Io,Jo,K,L) = FUP
!            ELSE
!              N = L+1-sgn
!              R = (F2(S,Io,Jo,K,N)-F2(S,Io,Jo,K,N-1))/X
!              IF (R.LE.0) FBNDMu(Io,Jo,K,L)=FUP
!              IF (R.GT.0) THEN
!                LIMITER = MAX(MIN(BetaLim*R,1._8),MIN(R,BetaLim))
!                CORR    = -0.5*(CDriftMu(Io,Jo,K,L)-sgn)*X
!                FBNDMu(Io,Jo,K,L)=FUP+LIMITER*CORR
!              ENDIF
!            ENDIF
!          ENDIF
!        ENDIF
!!!!!!! Calculate New Flux
!       IF (ORDER.EQ.1) THEN
!! DriftR
!         IF ((K.GT.1).AND.(Io.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftR(Io,Jo,K,L)*FBNDR(Io,Jo,K,L)
!&
!                                            +CDriftR(Im,Jo,K,L)*FBNDR(Im,Jo,K,L)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!! DriftP
!         IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftP(Io,Jo,K,L)*FBNDP(Io,Jo,K,L)
!&
!                                            +CDriftP(Io,Jm,K,L)*FBNDP(Io,Jm,K,L)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!         IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.EQ.NT)) F2(S,Io,1,K,L) = F2(S,Io,NT,K,L)
!! DriftE
!         IF ((K.GT.1).AND.(Io.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftE(Io,Jo,K,L)/WE(K)*FBNDE(Io,Jo,K,L) &
!                                            +CDriftE(Io,Jo,K-1,L)/WE(K)*FBNDE(Io,Jo,K-1,L)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!! DriftMu
!         IF ((L.LT.NPA).AND.(Io.GT.1).AND.(K.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftMu(Io,Jo,K,L)/WMU(L)*FBNDMu(Io,Jo,K,L) &
!                                            +CDriftMu(Io,Jo,K,L-1)/WMU(L)*FBNDMu(Io,Jo,K,L-1)
!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!         IF ((L.EQ.NPA).AND.(Io.GT.1).AND.(K.GT.1)) THEN 
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L-1)*FNHS(Io,Jo,L)*MU(L)/FNHS(Io,Jo,L-1)/MU(L-1)
!         ENDIF
!       ELSEIF (ORDER.EQ.-1) THEN
!! DriftMu
!         IF ((L.LT.NPA).AND.(Io.GT.1).AND.(K.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftMu(Io,Jo,K,L)/WMU(L)*FBNDMu(Io,Jo,K,L) &
!                                            +CDriftMu(Io,Jo,K,L-1)/WMU(L)*FBNDMu(Io,Jo,K,L-1)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!         IF ((L.EQ.NPA).AND.(Io.GT.1).AND.(K.GT.1)) THEN 
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L-1)*FNHS(Io,Jo,L)*MU(L)/FNHS(Io,Jo,L-1)/MU(L-1)
!         ENDIF
!! DriftE
!         IF ((K.GT.1).AND.(Io.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftE(Io,Jo,K,L)/WE(K)*FBNDE(Io,Jo,K,L) &
!                                            +CDriftE(Io,Jo,K-1,L)/WE(K)*FBNDE(Io,Jo,K-1,L)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!! DriftP
!         IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftP(Io,Jo,K,L)*FBNDP(Io,Jo,K,L) &
!                                            +CDriftP(Io,Jm,K,L)*FBNDP(Io,Jm,K,L)
!!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!         IF ((K.GT.1).AND.(Io.GT.1).AND.(Jo.EQ.NT)) F2(S,Io,1,K,L) = F2(S,Io,NT,K,L)
!! DriftR
!         IF ((K.GT.1).AND.(Io.GT.1)) THEN
!           F2(S,Io,Jo,K,L) = F2(S,Io,Jo,K,L)-CDriftR(Io,Jo,K,L)*FBNDR(Io,Jo,K,L) &
!                                            +CDriftR(Im,Jo,K,L)*FBNDR(Im,Jo,K,L)
!           if (F2(S,Io,Jo,K,L).lt.0) F2(S,Io,Jo,K,L)=1E-15
!         ENDIF
!       ENDIF
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO
!
!END SUBROUTINE DriftCalculations
!!==============================================================================

END MODULE ModRamDrift
