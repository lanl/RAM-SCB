!============================================================================
!    Copyright (c) 2016, Los Alamos National Laboratory`
!    All rights reserved.
!============================================================================

MODULE ModRamDrift
! Contains subroutines responsible for calculating flux changes due to different drifts

  use ModRamMain,      ONLY: Real8_
  use ModRamGrids,     ONLY: nR, nE, nPA
  use ModRamVariables, ONLY: FracCFL, DtDriftR, DtDriftP, DtDriftE, DtDriftMu

  implicit none

  integer, save :: QS
  real(kind=Real8_), save, ALLOCATABLE :: VR(:), P1(:), P2(:,:), MUDOT(:,:), EDOT(:,:), &
                                          CDriftR(:,:,:,:), CDriftP(:,:,:,:), &
                                          CDriftE(:,:,:,:), CDriftMu(:,:,:,:)
  !$OMP THREADPRIVATE(QS,VR,P1,P2,MUDOT,EDOT,CDriftR,CDriftP,CDriftE,CDriftMu)

  contains
!==============================================================================
  SUBROUTINE DRIFTEND
    ! Deallocates arrays needed for drift equations

    implicit none
    
    DEALLOCATE(VR, P1, P2, EDOT, MUDOT, CDriftR, CDriftP, CDriftE, CDriftMu)

  END SUBROUTINE DRIFTEND

!**************************************************************************
!                               DRIFTPARA
!                       Calculate drift parameters
!**************************************************************************
  SUBROUTINE DRIFTPARA(S)
    ! Initializes the parameters needed for the drift equations

    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: PI
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: RLZ, MDR, DPHI, EKEV, GREL, WMU, EBND, GRBND, &
                               PHIOFS, MU, species

    implicit none

    integer, intent(in) :: S
    real(DP) :: MUBOUN
    integer :: i, j, k, l

    if (.not.ALLOCATED(VR)) then
       ALLOCATE(VR(nR), P1(nR), P2(nR,nE), EDOT(nR,nE), MUDOT(nR,nPA))
       VR = 0.0; P1 = 0.0; P2 = 0.0; EDOT = 0.0; MUDOT = 0.0
       ALLOCATE(CDriftR(nR,nT,nE,nPa), CDriftP(nR,nT,nE,nPa), &
                CDriftE(nR,nT,nE,nPa), CDriftMu(nR,nT,nE,nPa))
       CDriftR = 0.0; CDriftP = 0.0; CDriftE = 0.0; CDriftMu = 0.0
    endif

    ! Electric field offset in radians and particle charge
    PHIOFS=0*PI/12.
    QS = species(S)%s_charge

    ! Parameters used in calculating radial and azimuthal drifts at boundaries
    DO I=1,NR
      DO J=1,NT
        VR(I)=DTs/MDR/(RLZ(I)+0.5*MDR)/2/DPHI
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

    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamParams,    ONLY: BetaLim
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, MDR, EKEV, GREL, DPHI, &
                               RLZ, CONF1, CONF2, FGEOS, VT, EIP, outsideMGNP

    implicit none

    integer, intent(in) :: S
    integer :: UR, i, j, j0, j1, k, l, n
    integer, ALLOCATABLE :: sgn(:,:)
    real(DP) :: p4, x, fup, r, corr, cgr1, cgr2, cgr3, ctemp, CGR,LIMITER
    real(DP), ALLOCATABLE :: F(:),FBND(:), CR(:,:)

    ALLOCATE(sgn(nR,nT),CR(nR,nT),F(NR+2),FBND(nR))
    sgn = 1; CR = 0.0; F = 0.0; FBND = 0.0

    DTDriftR(S) = 100000.0

    ! ExB Radial Drift
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

    ! Gradient Curvature Radial Drift
    DO K=1,NE
       P4=DTs*EKEV(K)*1000.0*(GREL(S,K)+1)/GREL(S,K)/DPHI/MDR/QS
       DO L=1,NPA
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
                ! ExB Radial Drift + Gradient Curvature Radial Drift
                CDriftR(I,J,K,L)=CR(I,J)+CGR
                if (outsideMGNP(i,j) == 0) then
                   ctemp = max(abs(CDriftR(I,J,K,L)),1E-10)
                   DTDriftR(S) = min( DTDriftR(S), FracCFL*DTs/ctemp)
                endif
                sgn(i,j) = 1
                if (CDriftR(I,J,K,L).lt.0) sgn(i,j) = -1
             END DO
             IF (sgn(nR,j).EQ.1) THEN
                FBND(1)=0.
                FBND(NR)=F(NR)
                UR=NR-1
             ELSE
                FBND(1)=F(2)
                UR=NR
                if (outsideMGNP(nR,j) == 1) then
                   F(nR+1) = 0._dp
                   F(nR+2) = 0._dp
                else
                   F(NR+1)=FGEOS(S,J,K,L)*CONF1*FNHS(NR,J,L)
                   F(NR+2)=FGEOS(S,J,K,L)*CONF2*FNHS(NR,J,L)
                endif
             END IF
             DO I=2,UR
                X=F(I+1)-F(I)
                FUP=0.5*(F(I)+F(I+1)-sgn(i,j)*X)
                IF (ABS(X).LE.1.E-27) FBND(I)=FUP
                IF (ABS(X).GT.1.E-27) THEN
                   N=I+1-sgn(i,j)
                   R=(F(N)-F(N-1))/X
                   IF (R.LE.0) FBND(I)=FUP
                   IF (R.GT.0) THEN
                      LIMITER=MAX(MIN(BetaLim*R,1.),MIN(R,BetaLim))
                      CORR=-0.5*(CDriftR(I,J,K,L)-sgn(i,j))*X
                      FBND(I)=FUP+LIMITER*CORR
                   END IF
                END IF
             END DO
             ! update the solution for next time step
             DO I=2,NR
                F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftR(I,J,K,L)*FBND(I)+CDriftR(I-1,J,K,L)*FBND(I-1)
                if (f2(s,i,j,k,l).lt.0) then
!                   write(*,*) 'in DRIFTR f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
                   f2(S,i,j,k,l)=1E-15
                endif
             END DO
          END DO
       END DO
    END DO
    DEALLOCATE(sgn,CR,F,FBND)

    RETURN
  END SUBROUTINE DriftR

!************************************************************************
!                            DRIFTP
!               Calculate changes due to azimuthal drift
!************************************************************************
  SUBROUTINE DRIFTP(S)

    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, FNIS, FNHS, BNES, VT, EIR, RLZ, MDR, DPHI, &
                               outsideMGNP

    implicit none

    integer, intent(in) :: S
    integer :: i, sgn, j, j1, k, l, n
    real(DP) :: x, fup, r, corr, ome, ctemp, GPA1,GPA2,LIMITER
    real(DP), ALLOCATABLE :: FBND(:),F(:)

    ALLOCATE(FBND(nT),F(nT))
    FBND = 0.0; F = 0.0

    DtDriftP(S) = 100000.0
    OME=7.3E-5 ! Earth's angular velocity [rad/s]
    DO L=1,NPA
       DO K=1,NE
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
                if (outsideMGNP(i,j) == 0) then
                   ctemp = max(abs(CDriftP(I,J,K,L)),1E-10)
                   DtDriftP(S) = min(DtDriftP(S), FracCFL*DTs/ctemp)
                endif
                sgn = 1
                if (CDriftP(I,J,K,L).lt.0) sgn = -1
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
                J1=J+1
                IF (J.EQ.NT) J1=2
                F2(S,I,J,K,L)=F2(S,I,J,K,L)-CDriftP(I,J,K,L)*FBND(J)+CDriftP(I,J-1,K,L)*FBND(J-1)
                if (f2(s,i,j,k,l).lt.0) then
 !                  write(*,*) 'in DRIFTP f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
                   f2(S,i,j,k,l)=1E-15
                endif
             END DO
             F2(S,I,1,K,L)=F2(S,I,NT,K,L)
          END DO
       END DO
    END DO
    
    DEALLOCATE(FBND,F)
    RETURN
  END SUBROUTINE DriftP

!**************************************************************************
!                       DRIFTE
!               Calculate energization along the drift path 
!**************************************************************************
  SUBROUTINE DRIFTE(S)

    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: CS, Q
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA, nS
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, BNES, FNIS, FNHS, dBdt, dIdt, EKEV, WE, RMAS, &
                               DPHI, RLZ, MDR, EBND, GREL, GRBND, DE, VT, EIR, &
                               EIP, outsideMGNP

    implicit none

    integer, intent(in) :: S
    integer :: i, sgn, j, j0, j2, k, l, n
    real(DP) :: ezero,gpa,gpr1,gpr2,gpr3,gpp1,gpp2,edt1, &
                drdt, dpdt, dbdt1, didt1, x, fup, r, corr, ome, &
                DRD1,DPD1,DRD2,DPD2, ctemp, LIMITER
    real(DP), ALLOCATABLE :: FBND(:),F(:),GRZERO(:)

    ALLOCATE(GRZERO(nS),FBND(nE),F(0:nE+2))
    GRZERO = 0.0; FBND = 0.0; F = 0.0

    DtDriftE(S)=10000.0
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
          DO L=1,NPA
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
                if (outsideMGNP(i,j) == 0) then
                   ctemp = max(abs(CDriftE(I,J,K,L)),1E-10)
                   DtDriftE(S) = min(DtDriftE(S), FracCFL*DTs*DE(K)/ctemp)
                endif
                sgn = 1
                if (CDriftE(I,J,K,L).lt.0) sgn = -1
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
!                   write(*,*) 'in DRIFTE f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
                   f2(S,i,j,k,l)=1E-15
                endif
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(GRZERO,FBND,F)
    RETURN
  END SUBROUTINE DriftE

!**************************************************************************
!                       DRIFTMU
!       Calculate pitch angle changes along the drift path
!**************************************************************************
  SUBROUTINE DRIFTMU(S)

    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: BetaLim
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, BNES, BOUNIS, BOUNHS, FNHS, dBdt, dIbndt, &
                               RLZ, DPHI, MDR, GREL, EKEV, DMU, WMU, MU, &
                               VT, EIP, EIR, outsideMGNP

    implicit none

    integer, intent(in) :: S
    integer :: i, j, j0, j1, k, l, n, ISGM
    real(DP) :: gmr1, gmr2, gmr3, gmp1, gmp2, drdm, dpdm, dbdt2, dibndt2, &
                x, fup, r, corr, ome, CMUDOT,EDT,DRM2,DPM2,DRM1,DPM1, &
                ctemp, LIMITER
    real(DP), ALLOCATABLE :: FBND(:),F(:)

    ALLOCATE(FBND(nPa),F(nPa))
    FBND = 0.0; F = 0.0

    DtDriftMu(S) = 10000.0
    OME=7.3E-5
    DO K=1,NE
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
                if (outsideMGNP(i,j) == 0) then
                   ctemp = max(abs(CDriftMu(I,J,K,L)),1E-32)
                   DtDriftMu(S)=min(DtDriftMu(S),FracCFL*DTs*DMU(L)/ctemp)
                endif
                ISGM = 1
                if (CDriftMu(I,J,K,L).lt.0.0) ISGM = -1
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
!                   write(*,*) 'in DRIFTMU f2<0 ', S,i,j,k,l, f2(S,i,j,k,l)
                   f2(S,i,j,k,l)=1E-15
                endif
             END DO
             F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)*FNHS(I,J,NPA)*MU(NPA)/FNHS(I,J,NPA-1)/MU(NPA-1)
          END DO
       END DO
    END DO

    DEALLOCATE(FBND,F)
    RETURN
  END SUBROUTINE DriftMu

!!==============================================================================

END MODULE ModRamDrift
