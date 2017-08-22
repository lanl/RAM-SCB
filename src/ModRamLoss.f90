MODULE ModRamLoss
! Contains subroutines for calculating particle loss and loss rates
! currently contains charge exchange and atmospheric loss

  use ModRamVariables, ONLY: ACHAR, ATLOS

  implicit none
  save

  contains
!**************************************************************************
!                              CEPARA
!       Routine calculates the charge exchange and atmosphere loss rates
!**************************************************************************
  SUBROUTINE CEPARA
    !!!! Module Variables
    use ModRamMain,      ONLY: Real8_, S
    use ModRamGrids,     ONLY: NE, NR, NT, NPA
    use ModRamTiming,    ONLY: DTs
    use ModRamVariables, ONLY: EKEV, V, RLZ, HDNS

    implicit none

    integer :: l, k, i, j
    real(kind=Real8_) :: x, y, alpha, taub

    ! Calculate charge exchange crosss-section of ring species S with H
    ! and then the charge exchange decay rate ACHAR
    IF (S.EQ.2) THEN
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF (X.LT.-2.) X=-2.
                   Y=-18.767-0.11017*X-3.8173e-2*X**2-0.1232*X**3-5.0488e-2*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross section of H+ in m2
                   ACHAR(I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSEIF (S.EQ.3) THEN
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF(X.LT.-2.) X=-2.
                   Y=-20.789+0.92316*X-0.68017*X**2+0.66153*X**3-0.20998*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross sect of He+ in m2
                   ACHAR(I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSEIF (S.EQ.4) THEN
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF(X.LT.-2.) X=-2.
                   Y=-18.987-0.10613*X-5.4841E-3*X**2-1.6262E-2*X**3-7.0554E-3*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross sect of O+ in m2
                   ACHAR(I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! Calculate the losses due to the collis with atmosphere
    DO K=2,NE
       DO I=2,NR
          TAUB=2*RLZ(I)/V(S,K)
          ATLOS(I,K)=EXP(-DTs/TAUB)
       ENDDO
    ENDDO

    RETURN
  END

!************************************************************************
!                       CHAREXCHANGE 
!               Calculate the decay due to charge exchange
!************************************************************************
  SUBROUTINE CHAREXCHANGE
    !!!! Module Variables
    use ModRamMain,      ONLY: S
    use ModRamGrids,     ONLY: NE, NT, NR, NPA
    use ModRamVariables, ONLY: F2

    implicit none

    integer :: i, j, k, l

    DO K=2,NE
       DO J=1,NT
          DO I=2,NR
             DO L=2,NPA
                F2(S,I,J,K,L)=F2(S,I,J,K,L)*ACHAR(I,J,K,L)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END


!************************************************************************
!                       ATMOSPHERIC LOSSES
!               Calculate loss to the atmosphere due to collisions
!************************************************************************
  SUBROUTINE ATMOL
    !!!! Module Variables
    use ModRamMain,      ONLY: S
    use ModRamGrids,     ONLY: NE, NT, NR, NPA
    use ModRamVariables, ONLY: F2, FNHS, UPA

    implicit none

    integer :: i, j, k, l

    DO K=2,NE
       DO J=1,NT
          DO I=2,NR
             DO L=UPA(I),NPA
                F2(S,I,J,K,L)=F2(S,I,J,K,L)*ATLOS(I,K)**(1/FNHS(I,J,L))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END

!==============================================================================
END MODULE ModRamLoss
