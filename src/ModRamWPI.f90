!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamWPI
! Contains subroutines related to wave particle interactions
! and electron lifetimes

  implicit none

  contains

!**************************************************************************
!                              WAVEPARA1
!               Life time due to WPI inside plasmasphere
!**************************************************************************
  SUBROUTINE WAVEPARA1(S)

    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NE, NR
    use ModRamVariables, ONLY: EKEV, LZ, WALOS1

    implicit none
    integer, intent(in) :: S

    integer :: i, ii, j, k
    real(DP):: TAU_WAVE,xE,xL,xlife, rEa(5),rL(8),rlife(5,8),clife(5)
    DATA rEa/0.2,0.5,1.0,1.5,2.0/
    DATA rL/5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.65/

    ! life time (as a function of energy and L)
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

    ! Calculates the losses due to the w-p interaction
    DO K=2,NE
      DO II=2,NR
        xE=EKEV(K)/1000.
        xL=LZ(II)
        if (xL.ge.1.65.and.xL.le.5.0) then
          do i=8,2,-1
            if (xL.ge.rL(i).and.xL.lt.rL(i-1)) then
              do j=1,5
                clife(j)=(log10(rlife(j,i-1))-log10(rlife(j,i)))/(rL(i-1)-rL(i))*(xL-rL(i))+log10(rlife(j,i))
                clife(j)=10.**clife(j)
              enddo
              EXIT
            endif
          enddo
        elseif (xL.gt.5.0) then
          do j=1,5
            clife(j)=(log10(rlife(j,1))-log10(rlife(j,2)))/(rL(1)-rL(2))*(xL-rL(1))+log10(rlife(j,1))
            clife(j)=10.**clife(j)
          enddo
        elseif (xL.lt.1.65) then
          do j=1,5
            clife(j)=(log10(rlife(j,7))-log10(rlife(j,8)))/(rL(7)-rL(8))*(xL-rL(8))+log10(rlife(j,8))
            clife(j)=10.**clife(j)
          enddo
        endif

        if (xE.ge.0.2.and.xE.lt.2.0) then
          do i=1,4
            if (xE.ge.rEa(i).and.xE.lt.rEa(i+1)) then
              xlife=(log10(clife(i+1))-log10(clife(i)))/(log10(rEa(i+1)) &
                    -log10(rEa(i)))*(log10(xE)-log10(rEa(i)))+log10(clife(i))
              xlife=10.**xlife
              EXIT
            endif
          enddo
        elseif (xE.lt.0.2) then
          xlife=(log10(clife(2))-log10(clife(1)))/(log10(rEa(2)) &
                -log10(rEa(1)))*(log10(xE)-log10(rEa(1)))+log10(clife(1))
          xlife=10.**xlife
        elseif (xE.ge.2.0) then
          xlife=(log10(clife(5))-log10(clife(4)))/(log10(rEa(5)) &
                -log10(rEa(4)))*(log10(xE)-log10(rEa(5)))+log10(clife(5))
          xlife=10.**xlife
        endif

        tau_wave=xlife*60.*60.*24.  ! day -> sec
        WALOS1(II,K)=tau_wave
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE WAVEPARA1

!**************************************************************************
!                              WAVEPARA2
!       Another life time due to diffusion not everywhere strong
!**************************************************************************
  SUBROUTINE WAVEPARA2(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: HMIN, RE, PI
    use ModRamGrids,     ONLY: NE, NR
    use ModRamVariables, ONLY: EKEV, LZ, RLZ, V, WALOS2, WALOS3
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: asind

    implicit none
    integer, intent(in) :: S

    integer :: i, k
    real(DP):: TAU_WAVE,EMEV,R1,R2,CLC
    real(DP), ALLOCATABLE :: CONE(:)

    ALLOCATe(CONE(nR+4))
    CONE = 0.0

    DO K=2,NE
      DO I=2,NR
        EMEV=EKEV(K)*0.001
        R1=0.08*EMEV**(-1.32)
        R2=0.4*10.**(2.*LZ(I)-6.+0.4*log10(29.*EMEV))
        tau_wave=min(R1,R2)
        tau_wave=1.0/tau_wave
        tau_wave=tau_wave*60.*60.*24.  ! day -> sec
        WALOS2(I,K)=tau_wave
      ENDDO
    ENDDO

    ! CONE - in degree
    DO I=1,NR
      CLC=(RE+HMIN)/RLZ(I)
      CONE(I)=ASIND(SQRT(CLC**3/SQRT(4.-3.*CLC)))
    END DO
    CONE(NR+1)=2.5
    CONE(NR+2)=1.5
    CONE(NR+3)=1.
    CONE(NR+4)=0.

    DO I=2,NR
      DO K=2,NE
        WALOS3(I,K)=64.0*LZ(I)*RE/35./(1-0.25)/ &
        SIN(CONE(I)*PI/180.)/SIN(CONE(I)*PI/180.)/V(S,K)
      ENDDO
    ENDDO

    DEALLOCATE(CONE)
    RETURN
  END SUBROUTINE WAVEPARA2

!*************************************************************************
!                              WAPARA_EMIC 
!       Routine reading normalized PA diffusion coeff
!          based on B. Ni's diffusion coefficients.
!**************************************************************************
  SUBROUTINE WAPARA_EMIC

    use ModRamVariables, ONLY: Daa_emic_h, Daa_emic_he, EKEV_emic, &
                               fp2c_emic, PAbn, DP, LZ, Ihs_emic, Ihes_emic
    use MOdRamGrids,     ONLY: ENG_emic, NCF_emic, NPA_emic, NLZ_emic, NR, NPA
    use ModIoUnit,       ONLY: UnitTMP_
    use ModRamMain,      ONLY: PathRamIn
    use ModRamGSL,       ONLY: GSL_Interpolation_1D

    implicit none
    integer :: jfpc, iR, L, iE, IER,I, ipa, iu,im,ih
    character(len=2) :: fpetofce
    character :: LL
    character*100 :: filename
    character :: ind1
    real(DP) :: Y, YY1, YY2, Lz_emic(NLZ_emic), PA_emic(NPA_emic)
    real(DP) :: Daahe_emic_tmp(NPA_emic), Daah_emic_tmp(NPA_emic)
    real(DP), allocatable :: Daahe_emic(:,:,:), Daah_emic(:,:,:)


    ! read and store the EMIC related pitch angle diffusion coefficients
    allocate(Daahe_emic(NLZ_emic, ENG_emic, NPA_emic), &
             Daah_emic( NLZ_emic, ENG_emic, NPA_emic))
    
    ! pitch angle in the Daa file
    do L=1, NPA_emic
       PA_emic(L) = L
    end do

    ! Lz in the Daa file
    do iR=1, NLZ_emic
       Lz_emic(iR) = 3+(iR-1)
    end do
    
    ! energies in the Daa file
    do iE=1, ENG_emic
       EkeV_emic(iE) = 10.**(-1.+0.1*(iE-1))
    end do
    
26  format(89(1PE16.7))
    
    ! fp2c in the Daa file
    do jfpc = 1, NCF_emic
       fp2c_emic(jfpc) = jfpc*2
       
       write(fpetofce,'(I2.2)')jfpc*2
       
       do iR=3, 7
          write(LL,'(I1.1)')iR
          filename = trim(PathRamIn)//'/HBand/EMIC_Proton_HBand_L'//LL//&
               '_fpetofce'//fpetofce//'_Daa.txt'
          
          open(unit=unittmp_,file=filename, status='old')
          do iE=1,ENG_emic
             read(unittmp_,26)(DaaH_emic_tmp(ipa), ipa=1,NPA_emic)
             ! daa_emic: interpolate it into ram pitch angle grid
             where(Daah_emic_tmp < 1.0e-31)Daah_emic_tmp = 1.0e-31
             Daah_emic(iR,iE,1) = DaaH_emic_tmp(NPA_emic) ! around 89 degree
             do L=2,NPA-1
                if(PAbn(L) .GE. 1.0 .and. PAbn(L) .LE. 89)then
                   call GSL_Interpolation_1D(PA_emic, DaaH_emic_tmp, &
                        PAbn(L), Daah_emic(iR,iE,L), IER)                  
                else
                   Daah_emic(iR,iE,L) = 1.0e-31
                end if
             end do
          end do
          
          close(unittmp_)
          
          open(unit=unittmp_,file=trim(PathRamIn)//'/HeBand/EMIC_Proton_HeBand_L'//LL//&
               '_fpetofce'//fpetofce//'_Daa.txt', status='old')
          do iE=1,ENG_emic
             read(unittmp_,26)(DaaHe_emic_tmp(ipa), ipa=1,NPA_emic)
             ! daa_emic: interpolate it into ram pitch angle grid
             where(Daahe_emic_tmp < 1.0e-31)Daahe_emic_tmp = 1.0e-31
             Daahe_emic(iR,iE,L) = DaaHe_emic_tmp(NPA_emic)                   
             do L=2,NPA-1
                if(PAbn(L) .GE. PA_emic(1) .and. PAbn(L) .LE. PA_emic(NPA_emic))then
                   call GSL_Interpolation_1D(PA_emic, DaaHe_emic_tmp, &
                        PAbn(L), Daahe_emic(iR,iE,L), IER)
                else
                   Daahe_emic(iR,iE,L) = 1.0e-31
                end if
             end do
          end do
          close(unittmp_)
       end do
       
       do iE=1,ENG_emic
          do L=2,NPA-1
             do I=2,NR
                if(maxval(Daah_emic(:,IE,L)) > 1.0e-31 .and. Lz(I) .GE.Lz_emic(1) .and. LZ(I) .LE. Lz_emic(NLz_emic))then
                   call GSL_Interpolation_1D(Lz_emic, daah_emic(:,iE, L),  Lz(I), YY1, IER)
                   Daa_emic_h(I,iE,L,jfpc) = YY1
                   
                else
                   Daa_emic_h(I,iE,L,jfpc) = 1.0e-31
                end if
                
                if(maxval(Daahe_emic(:,IE,L)) > 1.0e-31 .and. Lz(I) .GE.Lz_emic(1) .and. LZ(I) .LE. Lz_emic(NLz_emic))then
                   call GSL_Interpolation_1D(Lz_emic, daahe_emic(:,iE, L),  Lz(I), YY1, IER)
                   Daa_emic_he(I,iE,L,jfpc) = YY1
                else
                   Daa_emic_he(I,iE,L,jfpc) = 1.0e-31
                end if
                
             end do
          end do
       end do
       
    end do
    
    where(Daa_emic_h  ==0.0)Daa_emic_h = 1.0e-31
    where(Daa_emic_he ==0.0)Daa_emic_he= 1.0e-31

    deallocate(Daahe_emic, Daah_emic)


    !Read and store EMIC WAVE Intensity from Saikin model
    do iu=1,4
       write(ind1,'(I1.1)')iu
       filename =  trim(PathRamIn)//'/EMIC_model/EMIC_H_intensity_AE_'//ind1//'.txt'
       open(unit=unittmp_,file=filename, status='old')
       do im=1,20
          read(unittmp_,*)(Ihs_emic(iu,im,ih), ih=1,25)
       enddo
       close(unittmp_)

       filename = trim(PathRamIn)//'/EMIC_model/EMIC_He_intensity_AE_'//ind1//'.txt'
       open(unit=unittmp_,file=filename, status='old')
       do im=1,20
          read(unittmp_,*)(Ihes_emic(iu,im,ih), ih=1,25)
       enddo
       close(unittmp_)       
    enddo
    
  END SUBROUTINE WAPARA_EMIC

!*************************************************************************
!                              WAPARA_HISS
!       Routine reading normalized Energy & PA hiss diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_HISS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRamIn
    use ModRamGrids,     ONLY: NR, NCF, ENG, NPA
    use ModRamParams,    ONLY: HissFilePath
    use ModRamVariables, ONLY: LZ, fpofc, ndaaj, ENOR, ndvvj
    !!!! Share Modules
    use ModIoUnit,   ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i, ix, kn, l
    character(len=80) HEADER
    character(len=3) ST4
    character(len=2) ST3, ST2

    ST2 = '_e'
    DO I=1,NR
      fpofc(1)=2.
      write(ST4,'(I3.3)') INT(LZ(I)*100)
      DO IX=1,NCF
        write(ST3,'(I2.2)') INT(fpofc(ix))
        OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/whis_L'//ST4//'_'//ST3//ST2//'.aan',STATUS='old')
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ENOR(KN)
          read(UNITTMP_,27)(ndaaj(i,kn,l,ix),l=1,npa)
        ENDDO
        ndvvj = 0
        IF (IX.LT.NCF) fpofc(ix+1)=fpofc(ix)+4.

        DO KN=1,ENG
          DO L=1,NPA
            if (ndaaj(i,kn,l,ix).lt.1e-20) ndaaj(i,kn,l,ix)=1e-20
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

17  FORMAT(E13.4)
20  FORMAT(A80)
27  FORMAT(80(1PE12.3))

    RETURN
  END SUBROUTINE WAPARA_HISS

!*************************************************************************
!                              WAPARA_CHORUS
!       Routine reading bounce-aver PA wave diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_CHORUS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, PathRamIn
    use ModRamGrids,     ONLY: NR, NT, ENG, NPA
    use ModRamVariables, ONLY: KP, ECHOR, BDAAR
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i, j, kn, l, ikp
    real(DP), ALLOCATABLE :: RLDAA(:,:),RUDAA(:,:)
    character(len=1) ST3
    character(len=2) ST2
    character(len=80) HEADER

    ST2 = '_e'

    ALLOCATE(RLDAA(ENG,NPA),RUDAA(ENG,NPA))
    RLDAA = 0.0; RUDAA = 0.0

    ikp=INT(KP)
    IF (ikp.gt.4) ikp=4
    write(ST3,'(I1.1)') ikp

    OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/wlowcho_Kp'//ST3//ST2//'.aan',STATUS='old')
    READ(UNITTMP_,20) HEADER
    print*,'in WAPARA_CHORUS: ',trim(HEADER)

    DO I=1,NR
      DO J=1,NT
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ECHOR(KN)
          read(UNITTMP_,27)(RLDAA(kn,l),l=1,npa)
        ENDDO
        DO KN=1,ENG
          DO L=1,NPA
            if (RLDAA(kn,l).lt.1e-20) RLDAA(kn,l)=1e-20
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

    OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/wuppcho_Kp'//ST3//ST2//'.aan',STATUS='old')
    READ(UNITTMP_,20) HEADER

    DO I=1,NR
      DO J=1,NT
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ECHOR(KN)
          read(UNITTMP_,27)(RUDAA(kn,l),l=1,npa)
        ENDDO
        DO KN=1,ENG
          DO L=1,NPA
            if (RUDAA(kn,l).lt.1e-20) RUDAA(kn,l)=1e-20
            BDAAR(i,j,kn,l)=(RLDAA(kn,l)+RUDAA(kn,l))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

17  FORMAT(E13.4)
20  FORMAT(A80)
27  FORMAT(80(1PE12.3))

    DEALLOCATE(RLDAA,RUDAA)

    RETURN
  END SUBROUTINE WAPARA_CHORUS

! *************************************************************************
!                              WAPARA_Kp
!       Interpolate Kp-dependent diffusion cofficients
!**************************************************************************
  SUBROUTINE WAPARA_Kp(S)

    use ModRamVariables, ONLY: KP, CDAAR, CDAAR_chorus, NKpDiff, Kp_chorus

    implicit none
    integer, intent(in) :: S

    integer :: i1,i2

    if (Kp.gt.maxval(Kp_chorus)) then
      CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,NKpDiff)
    else
      i1 = minloc(abs(Kp-Kp_chorus),dim=1)
      if (Kp.lt.Kp_chorus(i1)) then
        i1 = i1-1
      end if
      i2 = i1+1
      ! Linear interpolation of Kp
      CDAAR(:,:,:,:) = (Kp-Kp_chorus(i1))*CDAAR_chorus(:,:,1:35,:,i2) &
                     + (Kp_chorus(i2)-Kp)*CDAAR_chorus(:,:,1:35,:,i1)
      CDAAR = CDAAR/(Kp_chorus(i2) - Kp_chorus(i1))
    end if
    CDAAR(:,25,:,:) = CDAAR(:,1,:,:)

    return
  end SUBROUTINE WAPARA_Kp

! *************************************************************************
!                              WAPARA_BAS
!       Routine reading normalized Energy & PA wave diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_BAS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: DoUseKpDiff, BASFilePath
    use ModRamGrids,     ONLY: NPA, NT, NE, NR
    use ModRamVariables, ONLY: MU, nR_Dxx, nE_Dxx, nPa_Dxx, RCHOR_Dxx, &
                               TCHOR_Dxx, ECHOR_Dxx, CDAAR_chorus, &
                               PACHOR_Dxx, CDAAR, nKpDiff, nT_Dxx
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: ACOSD
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i,j,k,l,nkp,nloop
    character(len=32) :: H1,H2,H3,nchar
    character(len=200) :: fname
    real(DP), ALLOCATABLE :: Dxx_hold(:,:,:), PA(:)

    ALLOCATE(Dxx_hold(NR_Dxx,NE_Dxx,NPA_Dxx), PA(NPA))
    Dxx_hold = 0.0; Pa = 0.0

    write(*,*) "Starting WAPARA_BAS"

    DO L=1,NPA
      PA(L)=ACOSD(MU(L))
    END DO

    if (DoUseKpDiff) then
      nloop = NKpDiff
    else
      nloop = 1
    endif

    do nkp=1,nloop
      write(nchar,'(i1)') nkp-1
      fname = trim(BASFilePath)//'bav_diffcoef_chorus_rpa_Kp'//trim(nchar)//'.PAonly.dat'
      OPEN(UNIT=UNITTMP_,FILE=trim(fname),STATUS='old')
      ! First skip over header
      do i=1,12,1
        read(UNITTMP_,*)
      end do
      do i=1,NR_Dxx
        do j=1,NT_Dxx
          do k=1,NE_Dxx
            read(UNITTMP_,'(A19,F6.4,A9,F6.4,A12,F8.2)') H1,RCHOR_Dxx(i), H2, &
                           TCHOR_Dxx(j), H3, ECHOR_Dxx(k)
            do l=1,NPA_Dxx
              read(UNITTMP_,'(F15.6,3E18.6)') PACHOR_Dxx(l),CDAAR_chorus(i,j,k,l,nkp)
            end do
            ! skip over blank lines
            do l=1,4
              read(UNITTMP_,*)
            end do
          end do
        end do
      end do
      CLOSE(UNIT=UNITTMP_)
      ! Interpolate onto L-shell, energy and pitch angle, assuming time is
      ! normal
      if ((NT.ne.25).and.(NE.ne.35)) then
        write(*,*) "We have a problem... assuming NT=25 & NE=35"
        stop
      end if

      if ((NR_Dxx.eq.NR).and.(NPA.eq.NPA_Dxx)) then
        write(*,*) "No interpolation of diffusion coeff for", fname
      end if

    end do   ! end NKpDiff loop

    ! Initialization
    CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,1)

    write(*,*) "Finished WAPARA_BAS"

    DEALLOCATE(Dxx_hold, PA)
    RETURN
  END SUBROUTINE WAPARA_BAS

!************************************************************************
!                       WAVELO
!       Calculate loss due to waves everywhere using electron lifetimes
!************************************************************************
  SUBROUTINE WAVELO(S)

    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: DoUsePlasmasphere
    use ModRamGrids,     ONLY: NE, NR, NT, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, KP, LZ, IP1, IR1, EKEV, NECR, WALOS1, &
                               WALOS2, WALOS3

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l, j1, i1
    real(DP) :: TAU_LIF,Bw,Kpmax
    real(DP), ALLOCATABLE :: RLpp(:)

    ALLOCATE(RLpp(nT))
    RLpp = 0.0

    Bw=30.
    IF (KP.GE.4.0) Bw=100.                                 ! pT
    Kpmax=KP  ! need Kpmax from previous 12 hrs, TO FIX!!!
    DO J=1,NT
      RLpp(J)=5.39-0.382*KPmax  ! PP from Moldwin et al. [2002]
      DO I=2,NR
        IF (DoUsePlasmasphere.and.NECR(I,J).gt.50.) RLpp(J)=LZ(I)
      ENDDO
    ENDDO

    DO K=2,NE
      DO I=2,NR
        DO J=1,NT
          DO L=2,NPA
            IF (LZ(I).LE.RLpp(J)) THEN
              TAU_LIF=WALOS1(I,K)*((10./Bw)**2)
            ELSEIF (LZ(I).GT.RLpp(J)) THEN
              IF (EKEV(K).LE.1000.) THEN
                TAU_LIF=WALOS2(I,K)*(1+WALOS3(I,K)/WALOS2(I,K))
                if (ekev(k).le.1.1) then
                  tau_lif=tau_lif*37.5813*exp(-1.81255*ekev(k))
                elseif (ekev(k).gt.1.1.and.ekev(k).le.5.) then
                  tau_lif=tau_lif*(7.5-1.15*ekev(k))
                else
                  tau_lif=tau_lif
                endif
              ELSEIF(EKEV(K).GT.1000.) THEN
                TAU_LIF=5.*3600*24/KP
              ENDIF
            ENDIF
            F2(S,I,J,K,L)=F2(S,I,J,K,L)*EXP(-DTs/TAU_LIF)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(RLpp)
    RETURN
  END SUBROUTINE WAVELO

!*************************************************************************
!                               WPADIF
!     Routine calculates the decay of the distribution function
!        due to WPI pitch angle diffusion using implicit scheme
!*************************************************************************
  SUBROUTINE WPADIF(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, PathRamOut
    use ModRamGrids,     ONLY: NT, NR, NE, NPA
    use ModRamTiming,    ONLY: Dts, T
    use ModRamParams,    ONLY: DoUseEMIC
    use ModRamVariables, ONLY: F2, FNHS, MU, DMU, WMU, ATAC, ATAW, species,&
                               ATAW_emic_h, ATAW_emic_he        
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(DP) :: AN,BN,GN,RP,DENOM
    real(DP), ALLOCATABLE :: F(:), RK(:), RL(:), FACMU(:)

    ALLOCATE(F(NPA),RK(NPA),RL(NPA),FACMU(NPA))
    F = 0.0; RK = 0.0; RL = 0.0; FACMU = 0.0

    DO J=1,NT
      DO I=2,NR
        DO K=2,NE
          DO L=2,NPA
            FACMU(L)=FNHS(I,J,L)*MU(L)
            F(L)=F2(S,I,J,K,L)/FACMU(L)
          END DO
          FACMU(1)=FNHS(I,J,1)*MU(1)
          F(1)=F(2) ! lower b.c.
          RK(1)=0.
          RL(1)=-1.
          DO L=2,NPA-1
            IF(species(S)%s_name.eq.'Electron')then
               AN=(ATAW(I,J,K,L)+ATAC(I,J,K,L))/DMU(L)       ! Hiss & chorus
               GN=(ATAW(I,J,K,L-1)+ATAC(I,J,K,L-1))/DMU(L-1) !  "
            ELSE
               IF(DoUseEMIC)then
                  AN = (ATAW_emic_h(I,J,K,L)+ATAW_emic_he(I,J,K,L))/DMU(L)       ! EMIC waves for ions
                  GN = (ATAW_emic_h(I,J,K,L-1)+ATAW_emic_he(I,J,K,L-1))/DMU(L-1) !  "   
               END IF
            END IF
            AN=AN*DTs/FACMU(L)/WMU(L)
            GN=GN*DTs/FACMU(L)/WMU(L)
            BN=AN+GN
            if (abs(-1-bn).lt.(abs(an)+abs(gn))) then
              open(UNITTMP_,file=trim(PathRamOut)//'diffcf_e.dat',status='unknown',position='append')
              write(UNITTMP_,*) ' T=',T/3600,' hr'
              write(UNITTMP_,*) 'i=',i,' j=',j,' k=',k,' l=',l
              write(UNITTMP_,*) 'an=',AN,' -1-bn=',(-1-BN),' gn=',GN
              close(UNITTMP_)
            endif
            RP=F(L)
            DENOM=BN+GN*RL(L-1)+1
            RK(L)=(RP+GN*RK(L-1))/DENOM
            RL(L)=-AN/DENOM
          END DO
          F2(S,I,J,K,NPA-1)=RK(NPA-1)/(1+RL(NPA-1)) ! upper b.c.
          DO L=NPA-2,1,-1
            F2(S,I,J,K,L)=RK(L)-RL(L)*F2(S,I,J,K,L+1)
          END DO
          F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)
          DO L=1,NPA
            F2(S,I,J,K,L)=F2(S,I,J,K,L)*FACMU(L)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(F,RK,RL,FACMU)
    RETURN
  END SUBROUTINE WPADIF

!*************************************************************************
!                               EMIC wave amplitude
!                           Based on Empirical Model by Saikin et al.
!*************************************************************************
  subroutine I_emic(AEind, xl, nmlt, I_H, I_He)

    use ModKind
    use ModRamVariables, ONLY: Ihs_emic, Ihes_emic
    implicit none

    integer, intent(in) :: xl, nmlt, AEind
    real(kind=real8_), intent(out) :: I_H,I_He

    I_H=0.0
    I_He=0.0

    if(AEind.ge.0.and.AEind.lt.100) then
       I_H = Ihs_emic(1,xl,nmlt)
       I_He= Ihes_emic(1,xl,nmlt)

    elseif(AEind.ge.100.and.AEind.lt.300) then
       I_H = Ihs_emic(2,xl,nmlt)
       I_He= Ihes_emic(2,xl,nmlt)
       
    elseif(AEind.ge.300.and.AEind.lt.400) then
       I_H = Ihs_emic(3,xl,nmlt)
       I_He= Ihes_emic(3,xl,nmlt)
       
    elseif(AEind.ge.400) then
       I_H = Ihs_emic(4,xl,nmlt)
       I_He= Ihes_emic(4,xl,nmlt)
       
    endif

  end subroutine I_emic
  
!==============================================================================  
END MODULE ModRamWPI
