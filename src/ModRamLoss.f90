!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamLoss
! Contains subroutines for calculating particle loss and loss rates
! currently contains charge exchange, atmospheric loss, and field line curvature
! scattering.

  implicit none

  contains

!**************************************************************************
!                              CEPARA
!       Routine calculates the charge exchange and atmosphere loss rates
!**************************************************************************
  SUBROUTINE CEPARA(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NE, NR, NT, NPA
    use ModRamTiming,    ONLY: DTs
    use ModRamVariables, ONLY: EKEV, V, RLZ, HDNS, CHARGE, ATLOS, species, ODNS, NDNS
    use ModRamIO,        ONLY: read_CHEX_file
    use ModRamGSL,       ONLY: GSL_Interpolation_1D
    implicit none

    integer, intent(in) :: S
    integer :: l, k, i, j, GSLErr
    real(DP) :: x, y, alpha, taub
    real(DP), allocatable :: CEXsig(:,:,:,:)

    CHARGE(S,:,:,:,:) = 1._dp

    ! Calculate charge exchange crosss-section of ring species S with H
    ! and then the charge exchange decay rate ACHAR
    select case(species(S)%s_name)
    case("Hydrogen")
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF (X.LT.-2.) X=-2.
                   Y=-18.767-0.11017*X-3.8173e-2*X**2-0.1232*X**3-5.0488e-2*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross section of H+ in m2
                   CHARGE(S,I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    case ("HeliumP1")
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF(X.LT.-2.) X=-2.
                   Y=-20.789+0.92316*X-0.68017*X**2+0.66153*X**3-0.20998*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross sect of He+ in m2
                   CHARGE(S,I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    case("OxygenP1")
       DO L=2,NPA
          DO K=2,NE
             DO I=2,NR
                DO J=1,NT
                   X=LOG10(EKEV(K))
                   IF(X.LT.-2.) X=-2.
                   Y=-18.987-0.10613*X-5.4841E-3*X**2-1.6262E-2*X**3-7.0554E-3*X**4
                   ALPHA=10.**Y*V(S,K)*HDNS(I,J,L)*DTs
                   ! 10**Y is cross sect of O+ in m2
                   CHARGE(S,I,J,K,L)=EXP(-ALPHA)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    case default
       if (trim(species(S)%CEX_file) == 'na') then
          ! No charge exchange loss
       else
          ! Charge exchange loss from cross sections in a file
          allocate(CEXsig(nR,nT,nE,nPa))

          ! Hydrogen Charge Exchange
          if (allocated(species(S)%CEX_nH)) then
             CEXsig = 1._dp
             DO K = 1, nE
                if (V(S,K) > species(S)%CEX_velocities(1)) then
                   call GSL_Interpolation_1D(species(S)%CEX_velocities, species(S)%CEX_nH, &
                                             V(S,K), Y, GSLerr)
                else
                   Y = species(S)%CEX_nH(1)
                endif
                DO L = 1, nPa
                   DO I = 1, nR
                      DO J = 1, nT
                         ALPHA = Y*V(S,K)*HDNS(I,J,L)*DTs
                         CEXsig(I,J,K,L) = EXP(-ALPHA)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             CHARGE(S,:,:,:,:) = CHARGE(S,:,:,:,:)*CEXsig(:,:,:,:)
          endif

          ! Oxygen Charge Exchange
          if (allocated(species(S)%CEX_nO)) then
             CEXsig = 1._dp
                if (V(S,K) > species(S)%CEX_velocities(1)) then
                   call GSL_Interpolation_1D(species(S)%CEX_velocities, species(S)%CEX_nO, &
                                             V(S,K), Y, GSLerr)
                else
                   Y = species(S)%CEX_nH(1)
                endif
             DO K = 1, nE
                DO L = 1, nPa
                   DO I = 1, nR
                      DO J = 1, nT
                         ALPHA = Y*V(S,K)*ODNS(I,J,L)*DTs
                         CEXsig(I,J,K,L) = EXP(-ALPHA)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             CHARGE(S,:,:,:,:) = CHARGE(S,:,:,:,:)*CEXsig(:,:,:,:)
          endif

          ! Nitrogen Charge Exchange
          if (allocated(species(S)%CEX_nN)) then
             CEXsig = 1._dp
                if (V(S,K) > species(S)%CEX_velocities(1)) then
                   call GSL_Interpolation_1D(species(S)%CEX_velocities, species(S)%CEX_nN, &
                                             V(S,K), Y, GSLerr)
                else
                   Y = species(S)%CEX_nH(1)
                endif
             DO K = 1, nE
                DO L = 1, nPa
                   DO I = 1, nR
                      DO J = 1, nT
                         ALPHA = Y*V(S,K)*NDNS(I,J,L)*DTs
                         CEXsig(I,J,K,L) = EXP(-ALPHA)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             CHARGE(S,:,:,:,:) = CHARGE(S,:,:,:,:)*CEXsig(:,:,:,:)
          endif

          deallocate(CEXsig)
       end if
    end select

    ! Calculate the losses due to the collis with atmosphere
    DO K=2,NE
       DO I=2,NR
          TAUB=2*RLZ(I)/V(S,K)
          ATLOS(S,I,K)=EXP(-DTs/TAUB)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE CEPARA

! *************************************************************************
  !                              FLC_Radius
  !       Routine determining the field line curvature radius
!**************************************************************************                    
  subroutine FLC_Radius
    use ModRamGrids,     ONLY: NR, NT, NPA
    use ModRamConst,     ONLY: PI
    use ModRamMain,      ONLY: PathRamOut, DP
    use ModRamTiming,    ONLY: TimeRamNow, TimeRamElapsed, Dt_bc
    use ModRamVariables, ONLY: r_curvEq, zeta1Eq, zeta2Eq, LZ, MLT
    use ModScbMain,      ONLY: REarth
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables

    use ModRamGSL,       ONLY: GSL_Interpolation_2D
    use ModRamIO,        ONLY: RamFileName
    use ModIoUnit,       ONLY: UnitTMP_

    implicit none

    integer :: i,j,k,GSLerr
    real(DP), dimension(nthe,npsi,nzeta) :: r_curv, dRcdS, dBdS, dsqRcdsqS, dsqBdsqS,&
         zeta1, zeta2, r2,ds,bb,  r, vbx,vby,vbz,ax,ay,az
    real(DP) :: dx, dy, dz, ddx_vbx, ddy_vbx, ddz_vbx, ddx_vby, ddy_vby, ddz_vby, &
         ddx_vbz,ddy_vbz,ddz_vbz
    real(DP) ::  epsl, alph, alpha0, MUBOUN,PAb, Nmin,tau_bounce(nPA)
    character(len=200) :: filename
    real(DP) :: coordpoint(2)
    !-------------------------------------------------------------------------------------
    ! The following calculation of the pitch angle diffusion coefficient due 
    ! to field line curvature scattering is based on Ebihara 2011, Young 2002, Young 2008. 
    !-------------------------------------------------------------------------------------

    ! calculate the FLC radius every Dt_bc seconds.
    if (mod(int(TimeRamElapsed), int(Dt_bc)) .gt. 1e-6) return

    do i=1,nthe
       r2(i,:,1:nzeta)   = x(i,:,1:nzeta)*x(i,:,1:nzeta)+y(i,:,1:nzeta)*y(i,:,1:nzeta)+z(i,:,1:nzeta)*z(i,:,1:nzeta)
       r(i,:,1:nzeta)    = sqrt(r2(i,:,1:nzeta))
    end do
    
    ! curvature radius: 1/Rc = (b dot delta)b = (vbx*d/dx + vby*d/dy + vbz*d/dz)(vbx, vby, vbz) 
    bb = sqrt(bx*bx+by*by+bz*bz)
    vbx = bx(:,:,1:nzeta)/bb(:,:,1:nzeta)
    vby = by(:,:,1:nzeta)/bb(:,:,1:nzeta)
    vbz = bz(:,:,1:nzeta)/bb(:,:,1:nzeta)

    do i=2,nthe-1
       do j=1,npsi
          do k=1,nzeta
             dx = x(i+1,j,k)-x(i,j,k)
             dy = y(i+1,j,k)-y(i,j,k)
             dz = z(i+1,j,k)-z(i,j,k)

             if (dx == 0.0)then ! may occur on dawn/dusk                                              
                ddy_vbx = (vbx(i+1,j,k)-vbx(i,j,k))/dy
                ddz_vbx = (vbx(i+1,j,k)-vbx(i,j,k))/dz
                
                ddy_vby = (vby(i+1,j,k)-vby(i,j,k))/dy
                ddz_vby = (vby(i+1,j,k)-vby(i,j,k))/dz
                
                ddy_vbz = (vbz(i+1,j,k)-vbz(i,j,k))/dy
                ddz_vbz = (vbz(i+1,j,k)-vbz(i,j,k))/dz
                
                ax(i,j,k) = vby(i,j,k) * ddy_vbx + vbz(i,j,k) * ddz_vbx
                ay(i,j,k) = vby(i,j,k) * ddy_vby + vbz(i,j,k) * ddz_vby
                az(i,j,k) = vby(i,j,k) * ddy_vbz + vbz(i,j,k) * ddz_vbz
                
             else if (dy == 0.0)then !may occur on noon/midnight                                      
                
                ddx_vbx = (vbx(i+1,j,k)-vbx(i,j,k))/dx
                ddz_vbx = (vbx(i+1,j,k)-vbx(i,j,k))/dz
                
                ddx_vby = (vby(i+1,j,k)-vby(i,j,k))/dx
                ddz_vby = (vby(i+1,j,k)-vby(i,j,k))/dz
                
                ddx_vbz = (vbz(i+1,j,k)-vbz(i,j,k))/dx
                ddz_vbz = (vbz(i+1,j,k)-vbz(i,j,k))/dz
                
                ax(i,j,k) = vbx(i,j,k) * ddx_vbx + vbz(i,j,k) * ddz_vbx
                ay(i,j,k) = vbx(i,j,k) * ddx_vby + vbz(i,j,k) * ddz_vby
                az(i,j,k) = vbx(i,j,k) * ddx_vbz + vbz(i,j,k) * ddz_vbz
                
             else
                ax(i,j,k) = &
                     vbx(i,j,k) * (vbx(i+1,j,k)-vbx(i,j,k))/dx +&
                     vby(i,j,k) * (vbx(i+1,j,k)-vbx(i,j,k))/dy +&
                     vbz(i,j,k) * (vbx(i+1,j,k)-vbx(i,j,k))/dz
                ay(i,j,k) = &
                     vbx(i,j,k) * (vby(i+1,j,k)-vby(i,j,k))/dx +&
                     vby(i,j,k) * (vby(i+1,j,k)-vby(i,j,k))/dy +&
                     vbz(i,j,k) * (vby(i+1,j,k)-vby(i,j,k))/dz
                az(i,j,k) = &
                     vbx(i,j,k) * (vbz(i+1,j,k)-vbz(i,j,k))/dx +&
                     vby(i,j,k) * (vbz(i+1,j,k)-vbz(i,j,k))/dy +&
                     vbz(i,j,k) * (vbz(i+1,j,k)-vbz(i,j,k))/dz
                
             end if
             
             ds(i,j,k) = sqrt(&
                  (x(i+1,j,k)-x(i,j,k))**2 + &
                  (y(i+1,j,k)-y(i,j,k))**2 + &
                  (z(i+1,j,k)-z(i,j,k))**2)*REarth ! across two grids                                 
          end do
       end do
    end do
    
    r_curv = 1./sqrt(ax**2+ay**2+az**2)*REarth
    r_curv(1,:,1:nzeta) = r_curv(2,:,1:nzeta)
    r_curv(nthe,:,1:nzeta) = r_curv(nthe-1,:,1:nzeta)

    do i=2,nthe-1
       dBdS(i,:,1:nzeta)  = bnormal*(bb(i+1,:,1:nzeta) - bb(i,:,1:nzeta))*1.0e-9/ds(i,:,1:nzeta)
       dRcdS(i,:,1:nzeta) = (r_curv(i+1,:,1:nzeta) - r_curv(i,:,1:nzeta))/ds(i,:,1:nzeta)
    end do
    dBdS(1,:,1:nzeta) = dBdS(2,:,1:nzeta)
    dRcdS(1,:,1:nzeta) = dRcdS(2,:,1:nzeta)
    
    ! get the zeta parameters: zeta1 = rc*(dsqRc/dsqS); zeta2 = Rc^2/(dsqB0/dsqS)      
    do i=2, nthe-2
       dsqBdsqS(i,:,1:nzeta)  = (dBdS(i+1,:,1:nzeta) - dBdS(i,:,1:nzeta))/ds(i,:,1:nzeta)
       dsqRcdsqS(i,:,1:nzeta) = (dRcdS(i+1,:,1:nzeta) - dRcdS(i,:,1:nzeta))/ds(i,:,1:nzeta)
    end do
    zeta1 = r_curv*dsqRcdsqS
    zeta2 = r_curv**2/(bnormal*1.0e-9*bb(:,:,1:nzeta)) * dsqBdsqS
  
    !filename = trim(PathRamOut)//RamFileName('r_curv_','dat',TimeRamNow)
    !open(unit=UnitTmp_, file=filename, status='unknown')
    !do i=1,npsi
    !   do j = 2, nzeta
    !      write(UnitTmp_,'(5e11.3)')x(nThetaEquator, i,j), y(nThetaEquator, i,j), r_curv(nThetaEquator,i,j), &
    !           zeta1(nThetaEquator,i,j),zeta2(nThetaEquator,i,j)
    !   end do
    !end do
    !close(UnitTmp_)
    
    do i=2, nR
       do j=1, nT-1
          coordPoint(1) = radRaw(i)* cos(azimRaw(j)*2*pi/24.-pi)
          coordPoint(2) = radRaw(i)* sin(azimRaw(j)*2*pi/24.-pi)
          call GSL_Interpolation_2D(x(nThetaEquator,1:npsi, 2:nzeta), y(nThetaEquator,1:npsi, 2:nzeta), &
               r_curv(nThetaEquator,1:npsi, 2:nzeta), coordPoint(1), coordPoint(2), r_curvEq(i,j), GSLerr)
          call GSL_Interpolation_2D(x(nThetaEquator,1:npsi, 2:nzeta), y(nThetaEquator,1:npsi, 2:nzeta), &
               zeta1(nThetaEquator,1:npsi, 2:nzeta), coordPoint(1), coordPoint(2), zeta1Eq(i,j), GSLerr)
          call GSL_Interpolation_2D(x(nThetaEquator,1:npsi, 2:nzeta), y(nThetaEquator,1:npsi, 2:nzeta), &
               zeta2(nThetaEquator,1:npsi, 2:nzeta), coordPoint(1), coordPoint(2), zeta2Eq(i,j), GSLerr)
       end do
    end do
    r_curvEq(:,nT) = r_curvEq(:,1) ! at MLT=0 & 24 
    r_curvEq(1,:)  = r_curvEq(2,:)    ! at the inner-most circle
    Zeta1Eq(:,nT) = Zeta1Eq(:,1) ! at MLT=0 & 24 
    Zeta1Eq(1,:)  = Zeta1Eq(2,:)    ! at the inner-most circle
    Zeta2Eq(:,nT) = Zeta2Eq(:,1) ! at MLT=0 & 24 
    Zeta2Eq(1,:)  = Zeta2Eq(2,:)    ! at the inner-most circle
    
    !filename = trim(PathRamOut)//RamFileName('r_curvEq_','dat',TimeRamNow)
    !open(unit=unittmp_, file=filename, status='unknown')
    !do i=2,nR
    !   do j = 1, nT-1
    !      write(unittmp_,'(2f16.4, 3e16.4)')Lz(i),MLT(j), r_curvEq(i,j), Zeta1Eq(i,j), Zeta2Eq(i,j)
    !   end do
    !end do
    !close(unittmp_)
    
  end subroutine FLC_Radius

 ! *************************************************************************
 !                              PARA_FLC
 !       Routine determining the field line curvature scattering coefficient
 !**************************************************************************
  subroutine PARA_FLC(S)
    use ModRamVariables, ONLY: MU, WMU, BOUNHS, BNES, EKEV, RMAS, &
                               LZ, MLT, V, r_curvEq, zeta1Eq, zeta2Eq, FLC_coef
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamConst,     ONLY: Q, CS
    use ModRamMain,      ONLY: PathRamOut, DP
    use ModRamTiming,    ONLY: TimeRamNow, T, TimeRamElapsed, Dt_bc

    use ModScbMain,      ONLY: REarth
    use ModRamIO,        ONLY: RamFileName
    use ModIoUnit,       ONLY: UnitTMP_
    use ModRamParams,    ONLY: DoWriteFLCDiffCoeff
    
    implicit none

    integer, intent(in) :: S
    integer :: i,j,k,l,lmin
    real(DP), dimension(nR,NT,nE) :: r_gyro,epsil
    real(DP) :: a1, a2, ba, ca, da, omegaa, Am
    real(DP), dimension(NPA) :: D, Daa, Nfactor
    real(DP) ::  epsl, alph, alpha0, MUBOUN,PAb, Nmin,tau_bounce(nPA)
    character(len=200) :: filename
    character(len=2)  :: ST2
    character(len=2), dimension(4) :: speciesString = (/'_h','_o','he','_e'/)
    !-------------------------------------------------------------------------------------
    ! The following calculation of the pitch angle diffusion coefficient due
    ! to field line curvature scattering is based on Ebihara 2011, Young 2002, Young 2008.
    !-------------------------------------------------------------------------------------

    ! calculate the FLC Daa every Dt_bc seconds.
    if (mod(int(TimeRamElapsed), int(Dt_bc)) .gt. 1e-6) return

    if(DoWriteFLCDiffCoeff)then
       st2 = speciesString(S)
       filename = trim(PathRamOut)//RamFileName('FLC_coeff_'//st2,'dat',TimeRamNow)
       open(unit=unittmp_, file=filename, status='unknown')
       write(unittmp_,'(a, f8.3, a)')'T= ', T/3600., ' field line curvature Daa (E, pitch_angle, Rg, Rc, Epsil, tau_bounce, Daa)'
    end if

    FLC_coef(S,:,:,:,:) = 0.0
    do i=1, NR
       do j=1, NT
          if(DoWriteFLCDiffCoeff)&
               write(unittmp_, '(a, 2x, f8.3, a, f8.3, 2(2x,f8.3))')'L=',LZ(i),'MLT=', MLT(j), zeta1Eq(i,j), zeta2Eq(i,j)
          do k=1,nE
             
             Nmin = 1.0e20
             lmin = 1
             Nfactor = 0.0
             tau_bounce = 0.0
             D = 0.0
             Daa = 0.0
             
             ! BNES in T unit                        
             r_gyro(i,j,k) = RMAS(S)*V(S,k)/abs(BNES(i,j)*Q) ! m*V/Bq 
             epsil(i,j,k) = r_gyro(i,j,k)/r_curvEq(i,j)
             
             if (epsil(i,j,k) .gt. 0.584) epsil(i,j,k) = 0.584 ! upper bound of epsil 
             
             epsl = epsil(i,j,k)
             if (epsl .ge. 0.1)then ! lower bound of epsil 
                a1 = -0.35533865 + 0.12800347 *epsl**(-1) + 0.0017113113*epsl**(-2)
                a2 =  0.23156321 + 0.15561211 *epsl**(-1) - 0.001860433 *epsl**(-2)
                ba = -0.51057275 + 0.93651781 *epsl**(-1) - 0.031690658 *epsl**(-2)
                ca =   1.0663037 - 1.0944973  *epsl**(-1) + 0.016679378 *epsl**(-2) - 0.000499*epsl**(-3)
                da = -0.49667826 -0.0081941799*epsl**(-1) + 0.0013621659*epsl**(-2)
                omegaa =1.0513540+   0.1351358*epsl       - 0.50787555  *epsl**2
                
                Am = exp(ca)*(zeta1Eq(i,j)**a1 * zeta2Eq(i,j)**a2 + da)
                
                
                do l=1,nPA-1
                   MUBOUN=MU(L)+0.5*WMU(L)
                   alph = acos(MUBOUN)
                   Nfactor(l) = (sin(omegaa*alph) * MUBOUN**ba)**(-1)
                   
                   ! bounce period
                   tau_bounce(l) = 4*LZ(i)*REarth*BOUNHS(i,j,l)/V(S,k) ! tau_bounce = 4R0*h/v
                   D(l) = Am**2/(2*tau_bounce(l)) ! tau_b is the bounce period, ~ h(pa)??

                   if(Nfactor(l) .le. Nmin) then
                      Nmin = Nfactor(l)
                      lmin = l
                   end if
                end do
                alpha0 = acos(MU(lmin)+0.5*WMU(lmin))

                do l=1, nPA-1
                   MUBOUN=MU(L)+0.5*WMU(L)
                   alph = acos(MUBOUN)
                   ! pitch angle diffusion coefficients from FLC 
                   Daa(l) = D(l) * Nfactor(lmin)**2 * sin(omegaa * alph)**2 * MUBOUN**(2*ba)/( (1-MUBOUN**2) * MUBOUN**2 )
                   
                   FLC_coef(S,i,j,k,l) =  Daa(l)*(1-MUBOUN**2) * MUBOUN * BOUNHS(i,j,l)
                end do
             end if
             
             if(DoWriteFLCDiffCoeff)then
                do l=1, nPA-1
                   write(unittmp_, '(2(2x, f11.3),5(2x,e11.3))') EkeV(k), acos(MU(l)+0.5*WMU(l)), r_gyro(i,j,k), &
                        r_curvEq(i,j)/REarth, epsil(i,j,k),tau_bounce(l),FLC_coef(S,i,j,k,l)
                end do
             end if

          end do
       end do
    end do
    
    close(unittmp_)
  end subroutine PARA_FLC

!************************************************************************
!                       CHAREXCHANGE 
!               Calculate the decay due to charge exchange
!************************************************************************
  SUBROUTINE CHAREXCHANGE(S)
    !!!! Module Variables
    use ModRamGrids,     ONLY: NE, NT, NR, NPA
    use ModRamVariables, ONLY: F2, CHARGE

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l

    DO K=2,NE
       DO J=1,NT
          DO I=2,NR
             DO L=2,NPA
                F2(S,I,J,K,L)=F2(S,I,J,K,L)*CHARGE(S,I,J,K,L)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE CHAREXCHANGE


!************************************************************************
!                       ATMOSPHERIC LOSSES
!               Calculate loss to the atmosphere due to collisions
!************************************************************************
  SUBROUTINE ATMOL(S)
    !!!! Module Variables
    use ModRamGrids,     ONLY: NE, NT, NR, NPA
    use ModRamVariables, ONLY: F2, FNHS, UPA, ATLOS

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l, u

    DO K=2,NE
       DO J=1,NT
          DO I=2,NR
             u = int(UPA(I),kind=4)
             DO L=u,NPA
                F2(S,I,J,K,L)=F2(S,I,J,K,L)*ATLOS(S,I,K)**(1/FNHS(I,J,L))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE ATMOL

!************************************************************************
!                       Flield line curvature scattering
!    Calculate loss to the atmosphere due to collisionless scattering
!************************************************************************
  subroutine FLCscatter(S)

    use ModRamGrids,     ONLY: NT, NR, NE, NPA
    use ModRamVariables, ONLY: F2, FNHS, MU, DMU, FLC_coef, DP, WMU
    use ModRamMain,      ONLY: PathRamOut
    use ModRamTiming,    ONLY: Dts, TimeRamElapsed, Dt_bc, T
    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(DP):: F(NPA),AN,BN,GN,RP,DENOM,RK(NPA),RL(NPA),FACMU(NPA)


    if (TimeRamElapsed .lt. Dt_bc) return ! let one cycle complete and have B distribution

    DO J=1,NT
       DO I=2,NR
          DO K=2,NE
             DO L=2,NPA
                FACMU(L)=FNHS(I,J,L)*MU(L)
                F(L)=F2(S,I,J,K,L)/FACMU(L)
             END DO
             FACMU(1)=FNHS(I,J,1)*MU(1)

             F(1)=F(2)                            ! lower b.c.   
             RK(1)=0.                             ! " 
             RL(1)=-1.                            ! " 

             DO L=2,NPA-1
                AN=FLC_coef(S,I,J,K,L)/DMU(L)     ! FLC_coef: pitch angle diffusion coefficent = Duu*u*h 
                GN=FLC_coef(S,I,J,K,L-1)/DMU(L-1)
                
                AN=AN*DTs/FACMU(L)/WMU(L)
                GN=GN*DTs/FACMU(L)/WMU(L)

                BN=AN+GN
                if (abs(-1-bn).lt.(abs(an)+abs(gn))) then
                   open(20,file=trim(PathRamOut)//'flc_cf_ion.dat',status='unknown', &
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

             F2(S,I,J,K,NPA-1)=RK(NPA-1)/(1+RL(NPA-1))    ! upper b.c.
             DO L=NPA-2,1,-1
                F2(S,I,J,K,L)=RK(L)-RL(L)*F2(S,I,J,K,L+1)
             END DO
             F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)
             DO L=1,NPA
                F2(S,I,J,K,L)=F2(S,I,J,K,L)*FACMU(L)
             ENDDO
          END DO
       END DO
    END DO

END subroutine FLCscatter

!==============================================================================
END MODULE ModRamLoss
