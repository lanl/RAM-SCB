MODULE ModSceIono

  use nrtype, ONLY: DP

  implicit none

  real(DP), allocatable :: x(:), y(:), rhs(:), b(:), Bnd_I(:), d_I(:), e_I(:), &
                           f_I(:), e1_I(:), f1_I(:)

  contains
!==================================================================================================
  subroutine FACs_to_fluxes_North

    !\
    ! The goal here is to convert the ionospheric FAC pattern into a 
    ! particle precipitation pattern, which can then be turned into
    ! a conductance pattern.
    !/

    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables, ONLY: Hall_to_Ped_Ratio, PolarCapPedConductance, IONO_Min_Ave_E, &
                               IONO_Min_EFlux, iono_north_im_eflux, iono_north_im_avee, &
                               iono_north_im_dis_eflux, iono_north_im_dis_avee, &
                               iono_north_ave_e, iono_north_eflux

    use nrtype, ONLY: DP

    implicit none

    real(DP) :: PolarCapHallConductance, PolarCap_AveE, PolarCap_EFlux 
    real(DP) :: MulFac_Dae, MulFac_Def, MulFac_ae, MulFac_ef
    real(DP), allocatable :: discrete_ef(:,:), discrete_ae(:,:), &
                             diffuse_ef(:,:), diffuse_ae(:,:)
    !---------------------------------------------------------------------------
    allocate(discrete_ef(Iono_nTheta,Iono_nPsi), discrete_ae(Iono_nTheta,Iono_nPsi), &
             diffuse_ef(Iono_nTheta,Iono_nPsi), diffuse_ae(Iono_nTheta,Iono_nPsi))

    Hall_to_Ped_Ratio = 1.5

    if (PolarCapPedConductance > 0.0) then
       PolarCapHallConductance = Hall_to_Ped_Ratio * PolarCapPedConductance
       PolarCap_AveE = (Hall_to_Ped_Ratio/0.45)**(1.0/0.85)
       PolarCap_EFlux = ((PolarCapPedConductance*(16.0 + PolarCap_AveE**2) / &
            (40.0*PolarCap_AveE))**2)/1000.0 ! convert ergs/cm2/s to W/m2
    else
       PolarCap_AveE  = IONO_Min_Ave_E
       PolarCap_EFlux = IONO_Min_EFlux
    endif

    MulFac_Dae = 1.0e22
    MulFac_Def = 1.0e19
    MulFac_ef = 0.3e6
    MulFac_ae = 4.0e-12

    !Yu: recalculate the coefficients (by comparing to J. Raeder's formulas. take F1=4)
    MulFac_Dae = 0.0267
    MulFac_Def = 6.0e17

    !\
    ! pass the diffuse energy from IM integrated over the loss cone due to WPI
    !/

    diffuse_ef = iono_north_im_eflux/1000.0 ! from ergs/cm^2s to W/m^2
    where (diffuse_ef < IONO_Min_EFlux) diffuse_ef = IONO_Min_EFlux

    diffuse_ae = iono_north_im_avee         ! keV
    where (diffuse_ae > 20.0) diffuse_ae = 20.0

    !\
    ! pass the discrete energy from IM using Te, Ne from RAM, following Zhang et al.2015 formulas
    !/
    ! it is only applied in Jr>0 region
    discrete_ef = iono_north_im_dis_eflux/1000.0
    where (discrete_ef < IONO_Min_EFlux) discrete_ef = IONO_Min_EFlux

    discrete_ae = iono_north_im_dis_avee
    where (discrete_ae > 20.0) discrete_ae = 20.0

    !\
    !/
    where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
    where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

    ! Let's weight the average energy by the number flux, which is ef/av
    iono_north_ave_e = (diffuse_ef + discrete_ef) &
                      /(diffuse_ef/diffuse_ae + discrete_ef/discrete_ae)

    ! The energy flux should be weighted by the average energy
    iono_north_eflux = (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae) * iono_north_ave_e
    where (iono_north_ave_e < IONO_Min_Ave_E) iono_north_ave_e = IONO_Min_Ave_E

    deallocate(discrete_ef, discrete_ae, diffuse_ef, diffuse_ae)
    return
  end subroutine FACs_to_Fluxes_North

!==================================================================================================
  subroutine ionosphere_conductance(Sigma0, SigmaH, SigmaP, SigmaThTh, SigmaThPs, &
                                    SigmaPsPs, dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                    dSigmaPsPs_dTheta, dSigmaThTh_dPsi, dSigmaThPs_dPsi, &
                                    dSigmaPsPs_dPsi, Eflux, Ave_E, Theta, Psi, nTheta, &
                                    nPsi, dTheta, dPsi, x, y, z)
  
    !\
    ! This subroutine computes the height-integrated field-aligned and
    ! Hall and Pedersen conductances for the ionosphere at each
    ! location of the discretized solution domain.  The gradients of
    ! these quantities are also computed.
    !/

    use ModRamVariables, ONLY: F107

    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables, ONLY: IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                               StarLightPedConductance, SAVE_NORTH_SigmaH, &
                               SAVE_NORTH_SigmaP, cosThetaTilt, sinThetaTilt

    use nrtype, ONLY: DP, cDegtoRad

    implicit none

    integer, intent(in)     :: nTheta, nPsi
    real(DP), intent(inout) :: Sigma0(:,:), SigmaH(:,:), SigmaP(:,:), SigmaThTh(:,:), &
                               SigmaThPs(:,:), SigmaPsPs(:,:), dSigmaThTh_dTheta(:,:), &
                               dSigmaThPs_dTheta(:,:), dSigmaPsPs_dTheta(:,:), &
                               dSigmaThTh_dPsi(:,:), dSigmaThPs_dPsi(:,:), &
                               dSigmaPsPs_dPsi(:,:), Eflux(:,:), Ave_E(:,:), Theta(:,:), &
                               Psi(:,:), x(:,:), y(:,:), z(:,:), dTheta(:), dPsi(:)
    integer  :: i, j
    logical  :: old
    real(DP) :: f107p53, f107p49, cos_limit, meeting_value_p, meeting_value_h, &
                SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, SigmaH_EUV_2, &
                SigmaP_EUV_2, SigmaH_STAR, SigmaP_STAR, sn, cs, sn2, cs2, cs3, &
                cs4, C, SigmaH_Particles, SigmaP_Particles

    real(DP), allocatable :: cos_SZA(:,:)

    allocate(cos_SZA(nTheta,nPsi))
 
    cos_SZA = (x*cosTHETATilt-z*sinTHETATilt)/sqrt(x**2 + y**2 + z**2)
  
    ! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
    ! So, let's just store them straight away:  
    f107p53 = f107**0.53
    f107p49 = f107**0.49
    cos_limit = cos(70.0*cDegToRad)
    meeting_value_p = f107p49*(0.34*cos_limit+0.93*sqrt(cos_limit))
    meeting_value_h = f107p53*(0.81*cos_limit+0.54*sqrt(cos_limit))
  
    do j = 1, nPsi
       do i = 1, nTheta
  
          Sigma0(i,j) = 1000.00
  
          if (cos_SZA(i,j) > 0) then
             SigmaH_EUV=f107p53*(0.81*cos_SZA(i,j)+0.54*sqrt(cos_SZA(i,j)))
             SigmaP_EUV=f107p49*(0.34*cos_SZA(i,j)+0.93*sqrt(cos_SZA(i,j)))
             SigmaH_SCAT = 1.00
             SigmaP_SCAT = 0.50
             if (cos_SZA(i,j) < cos_limit) then
                SigmaH_EUV_2 = meeting_value_h *   &
                     exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                SigmaP_EUV_2 = meeting_value_p *   &
                     exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
                SigmaH_EUV = (SigmaH_EUV + SigmaH_EUV_2)/2.0
                SigmaP_EUV = (SigmaP_EUV + SigmaP_EUV_2)/2.0
             endif
          else
             SigmaH_EUV = meeting_value_h *   &
                  exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
             SigmaP_EUV = meeting_value_p *   &
                  exp(-((cos_SZA(i,j)-cos_limit)**2.0)*15.0)
             SigmaH_SCAT = 1.00*(10.00**cos_SZA(i,j))
             SigmaP_SCAT = 0.50*(10.00**cos_SZA(i,j))
          end if
  
          SigmaH_STAR = StarLightPedConductance*2.0
          SigmaP_STAR = StarLightPedConductance
  
          !\
          ! Use Robinson's Formula to convert the Ave_E and E_Flux to SigmaP and SigmaH
          !/
  
          SigmaP_Particles = 40.0 * Ave_E(i,j) / (16.0 + Ave_E(i,j)*Ave_E(i,j))  &
                                  * sqrt(EFlux(i,j)*1000.0)
  
          SigmaH_Particles = 0.45 * (Ave_E(i,j)**0.85) * SigmaP_Particles
  
          SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + SigmaH_SCAT*SigmaH_SCAT &
                           + SigmaH_STAR*SigmaH_STAR + SigmaH_Particles*SigmaH_Particles)
  
          SigmaP_EUV = SigmaP_EUV*SigmaP_EUV + SigmaP_SCAT*SigmaP_SCAT + SigmaP_STAR*SigmaP_STAR
  
          SigmaP_Particles = SigmaP_Particles*SigmaP_Particles
  
          SigmaP(i,j) = sqrt(SigmaP_EUV + SigmaP_Particles)
       enddo
    enddo
  
    do j = 1, nPsi
       do i = 1, nTheta
  
          sn = sin(Theta(i,j))
          cs = cos(Theta(i,j))
          sn2= sn*sn
          cs2 = cs*cs
          cs3 = 1.00 + 3.00*cs2
          cs4 = sqrt(cs3)
          C = 4.00*Sigma0(i,j)*cs2 + SigmaP(i,j)*sn2
  
          SigmaThTh(i,j) = Sigma0(i,j)*SigmaP(i,j)*cs3/C
          SigmaThPs(i,j) = 2.00*Sigma0(i,j)*SigmaH(i,j)*cs*cs4/C
          SigmaPsPs(i,j) = SigmaP(i,j)+SigmaH(i,j)*SigmaH(i,j)*sn2/C
  
          dSigmaThTh_dTheta(i,j) = 0.00
          dSigmaThTh_dPsi(i,j) = 0.00
          dSigmaThPs_dTheta(i,j) = 0.00
          dSigmaThPs_dPsi(i,j) = 0.00
          dSigmaPsPs_dTheta(i,j) = 0.00
          dSigmaPsPs_dPsi(i,j) = 0.00
  
       end do
    end do
  
    do j = 1, nPsi
       if (j > 1 .and. j < nPsi ) then
          do i = 2, nTheta-1
             dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/dTheta(i)
             dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,j-1))/dPsi(j)
  
             dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/dTheta(i)
             dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,j-1))/dPsi(j)
  
             dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/dTheta(i)
             dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,j-1))/dPsi(j)
          end do
       else if (j == 1) then
          do i = 2, nTheta-1
             dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/dTheta(i)
             dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,j+1)-SigmaThTh(i,nPsi-1))/dPsi(j)
  
             dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/dTheta(i)
             dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,j+1)-SigmaThPs(i,nPsi-1))/dPsi(j)
  
             dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/dTheta(i)
             dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,j+1)-SigmaPsPs(i,nPsi-1))/dPsi(j)
          end do
       else
          do i = 2, nTheta-1
             dSigmaThTh_dTheta(i,j) = (SigmaThTh(i+1,j)-SigmaThTh(i-1,j))/dTheta(i)
             dSigmaThTh_dPsi(i,j) = (SigmaThTh(i,2)-SigmaThTh(i,j-1))/dPsi(j)
  
             dSigmaThPs_dTheta(i,j) = (SigmaThPs(i+1,j)-SigmaThPs(i-1,j))/dTheta(i)
             dSigmaThPs_dPsi(i,j) = (SigmaThPs(i,2)-SigmaThPs(i,j-1))/dPsi(j)
  
             dSigmaPsPs_dTheta(i,j) = (SigmaPsPs(i+1,j)-SigmaPsPs(i-1,j))/dTheta(i)
             dSigmaPsPs_dPsi(i,j) = (SigmaPsPs(i,2)-SigmaPsPs(i,j-1))/dPsi(j)
          end do
       end if
    end do
  
    deallocate(cos_SZA)
    return

  end subroutine ionosphere_conductance

!==================================================================================================
  subroutine ionosphere_currents(Jx, Jy, Jz, Ex, Ey, Ez, ETh, EPs, Ux, Uy, Uz, &
                                 PHI, SigmaThTh, SigmaThPs, SigmaPsPs, X, Y, Z, &
                                 Theta, Psi, dTheta, dPsi)
  
    !\
    ! For the calculated ionospheric potential solution,
    ! this routine determines the ionospheric currents and
    ! electric fields, as well as convection velocities.
    !/
  
    use ModRamTiming,    ONLY: TimeRamElapsed
    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi 
    use ModSceVariables, ONLY: IONO_TOLER, IONO_NORTH_JTh, IONO_NORTH_JPs, &
                               IONO_Radius, IONO_Height, Radius

    use nrtype, ONLY: DP

    use ModCoordTransform, ONLY: dir_to_xyz, cross_product
    use CON_planet_field,  ONLY: get_planet_field

    implicit none
  
    integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi
  
    real(DP), intent(inout) :: PHI(:,:), SigmaThTh(:,:), SigmaThPs(:,:), &
                               SigmaPsPs(:,:), Jx(:,:), Jy(:,:), Jz(:,:), &
                               Ex(:,:), Ey(:,:), Ez(:,:), ETh(:,:), EPs(:,:), &
                               Ux(:,:), Uy(:,:), Uz(:,:), X(:,:), Y(:,:), &
                               Z(:,:), Theta(:,:), Psi(:,:), dTheta(:), &
                               dPsi(:)
  
    integer  :: i, j
    real(DP) :: cosTheta, sinTheta, cosPhi, sinPhi, ER, JR, JTh, JPs, NormRadius
    real(DP) :: Xyz_D(3), b_D(3), Vp_D(3)
    !----------------------------------------------------------------------------
    ! Compute the ionospheric electric field.
    do j = 1, nPsi
       if (j > 1 .and. j < nPsi ) then
          do i = 2, nTheta-1
             sinTheta = sin(Theta(i,j))
             ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/(dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,j+1)-PHI(i,j-1))/(dPsi(j)*Radius*sinTheta)
          end do
          ETh(1,j) = -(PHI(2,j)-PHI(1,j))/(dTheta(1)*Radius)
          EPs(1,j) = EPs(2,j)
          ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/(dTheta(nTheta)*Radius)
          EPs(nTheta,j) = EPs(nTheta-1,j)
       else if (j == 1) then
          do i = 2, nTheta-1
             sinTheta = sin(Theta(i,j))
             ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/(dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,j+1)-PHI(i,nPsi-1))/(dPsi(j)*Radius*sinTheta)
          end do
          ETh(1,j) = -(PHI(2,j)-PHI(1,j))/(dTheta(1)*Radius)
          EPs(1,j) = EPs(2,j)
          ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/(dTheta(nTheta)*Radius)
          EPs(nTheta,j) = EPs(nTheta-1,j)
       else
          do i = 2, nTheta-1
             sinTheta = sin(Theta(i,j))
             ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/(dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,2)-PHI(i,j-1))/(dPsi(j)*Radius*sinTheta)
          end do
          ETh(1,j) = -(PHI(2,j)-PHI(1,j))/(dTheta(1)*Radius)
          EPs(1,j) = EPs(2,j)
          ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/(dTheta(nTheta)*Radius)
          EPs(nTheta,j) = EPs(nTheta-1,j)
       end if
    end do
  
    ! Compute the ionospheric currents convection velocities.
    do j = 1, nPsi
       do i = 1, nTheta
          cosTheta = cos(Theta(i,j))
          sinTheta = sin(Theta(i,j))
          cosPhi = cos(Psi(i,j))
          sinPhi = sin(Psi(i,j))
  
          if (i == nTheta) then
             ER = 0.00
          else
             ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh(i,j)
          end if
  
          Ex(i,j) = ER*sinTheta*cosPhi + ETh(i,j)*cosTheta*cosPhi - EPs(i,j)*sinPhi
          Ey(i,j) = ER*sinTheta*sinPhi + ETh(i,j)*cosTheta*sinPhi + EPs(i,j)*cosPhi
          Ez(i,j) = ER*cosTheta - ETh(i,j)*sinTheta
  
          JR = 0.00
          JTh =  SigmaThTh(i,j)*ETh(i,j) + SigmaThPs(i,j)*EPs(i,j)
          JPs = -SigmaThPs(i,j)*ETh(i,j) + SigmaPsPs(i,j)*EPs(i,j)
  
          IONO_NORTH_JTh(i,j) = JTh
          IONO_NORTH_JPs(i,j) = JPs
  
          Jx(i,j) = JR*sinTheta*cosPhi + JTh*cosTheta*cosPhi - JPs*sinPhi
          Jy(i,j) = JR*sinTheta*sinPhi + JTh*cosTheta*sinPhi + JPs*cosPhi
          Jz(i,j) = JR*cosTheta - JTh*sinTheta
  
          ! Calculate location in Cartesian coordinates
          call dir_to_xyz(SinTheta,CosTheta,SinPhi,CosPhi,Xyz_D)
          Xyz_D = Xyz_D * (IONO_Radius + IONO_Height) / IONO_Radius
          ! Get magnetic field and normalize it to unity
          call get_planet_field(TimeRamElapsed,Xyz_D,'SMG NORM',b_D)
  
          ! Get potential V = E x B/|B^2|
          b_D = b_D/sum(b_D**2)
          Vp_D = cross_product((/Ex(i,j), Ey(i,j), Ez(i,j)/), b_D)
  
          Ux(i,j) = Vp_D(1)
          Uy(i,j) = Vp_D(2)
          Uz(i,j) = Vp_D(3)
  
       end do
    end do
  
  end subroutine ionosphere_currents

!==================================================================================================
  subroutine ionosphere_jouleheating_ionflux(ETh, EPs, SigmaP, Joule, IonNumFlux)
  
    !\
    ! Joule heating is determined by SigmaP * E^2
    !/
  
    use ModRamTiming,    ONLY: TimeRamElapsed
    use ModSceGrids,     ONLY: IONO_nTheta, IONO_nPsi
    use ModSceVariables, ONLY: IONO_NORTH_Theta, IONO_NORTH_Psi, IONO_Radius, Radius, &
                               IONO_Height

    use nrtype, ONLY: DP 

    use CON_planet_field
    use ModCoordTransform, ONLY: sph_to_xyz

    implicit none

    real(DP), intent(inout) :: SigmaP(:,:), ETh(:,:), EPs(:,:), Joule(:,:), IonNumFlux(:,:)  
  
    integer  :: i, j, iHemisphere
    real(DP) :: bIono, B, ratioOH

    real(DP) :: B_D(3), bIono_D(3), XyzIono_D(3), Xyz_tmp(3)
    real(DP), parameter :: height_fast = 4.0e6

    !\
    ! Joule heating is assumed to be equal to Poynting Flux, 
    ! since according to the Poynting Theorem, divergence of S = J*E,
    ! which is integral(S dA) = integral(J*E dV),
    ! but J in the formula is current density with a unit of A/m^3.
    ! In this routine, J is height-integrated (A/m^2), therefore, 
    ! J * E = S in this code.
    !/
  
    Joule(:,:) =0.0
    IonNumFlux(:,:) = 0.0
    ! Joule heating (or Poynting flux) at ionosphere altitude
    Joule = SigmaP * (ETh**2+EPs**2)
  
    !\
    ! Ion number flux at ionosphere altitude, based on Strangeway et al(2005):
    ! At an altitude of 4000km(where the FAST orbit), f = 2.142e7 * (S)^1.265, 
    ! where unit of f is 1/(cm^2 * s), and unit of S is mW/m^2
    !/
  
    ! Mapping to the ionosphere altitude along dipole magnetic field lines
    do i = 1, IONO_nTheta
       do j = 1, IONO_nPsi
          call sph_to_xyz(Radius, IONO_NORTH_Theta(i,j), IONO_NORTH_Psi(i,j), XyzIono_D)
          call get_planet_field(TimeRamElapsed, XyzIono_D, 'SMG', bIono_D)
          bIono = sqrt(sum(bIono_D**2))
  
          call map_planet_field(TimeRamElapsed, XyzIono_D, 'SMG', (height_fast + IONO_Radius),Xyz_tmp, iHemisphere)
          call get_planet_field(TimeRamElapsed, Xyz_tmp, 'SMG', B_D)
          b = sqrt(sum(B_D**2))
  
          ! Flux at ionosphere altitude by mapping both S and f down to the ionosphere 
          ! from 4000km, caution in the units (flux in /m^2/s, jouleheating in W/m^2)
          ! The total ion number flux = No * Vo
          IonNumFlux(i,j) = 2.142e7 * (Joule(i,j)*1.0e3)**1.265 * (b/bIono)**0.265 * 1.0e4
       end do
    end do
  
  end subroutine ionosphere_jouleheating_ionflux

!==================================================================================================
  subroutine ionosphere_solver(Jr, SigmaThTh, SigmaThPs, SigmaPsPs, &
                               dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                               dSigmaPsPs_dTheta, dSigmaThTh_dPsi, &
                               dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                               Theta, Psi, dTheta, dPsi, Phi_C)
  
    !
    ! This subroutine solves for the ionospheric potential PHI
    ! using the field aligned currents Jr and the conductivity tensor Sigma
    ! and its various derivatives as input.
    !
    ! The idea is that the Jr current is balanced by the divergence of the
    ! horizontal currents. The horizontal currents are Sigma . E, 
    ! the height integrated conductivity Sigma is 2x2 antisymmetric matrix
    ! acting on the Theta and Phi components of the electric field.
    ! The electric field is assumed to be a potential field: E = Grad Phi.
    !
    ! We are solving 
    !
    !    Div( Sigma . Grad Phi ) = - Jr
    ! 
    ! in spherical coordinates. 
    !
    ! This leads to a penta-diagonal linear equation:
    !  
    !  C_A*Phi(i,j) + C_B*Phi(i-1,j) + C_C*Phi(i+1,j) +
    !                 C_D*Phi(i,j-1) + C_E*Phi(i,j+1) = RHS
    !
    ! To avoid division by zero at the poles, the equation is 
    ! multiplied by (sin(Theta)*Radius)**2, thus the 
    ! RHS = Jr * (sin(Theta)*Radius)**2
    !
    ! The linear elliptic PDE for the
    ! electric field potential PHI is defined on the domain 
    ! 0 < Theta < ThetaMax  for the northern hemisphere, and
    ! ThetaMin < Theta < PI for the southern hemisphere), and
    ! 0 < Psi < 2 PI.
    !
    ! The following boundary conditions are applied:
    !
    !      PHI(ThetaMax,Psi) = PHI(ThetaMin,Psi) = 0,
    !
    !      PHI(Theta,0) = PHI(Theta,2*PI).
    ! 
    ! There is no boundary at the poles, but to avoid numerical 
    ! difficulties, the cell value at the pole is replaced with the 
    ! average of the first neighbors:
    !
    !      PHI(0,Psi)  = average( PHI(dTheta, Psi) )
    !
    !      PHI(PI,Psi) = average( PHI(PI-dTheta, Psi) )
    !
    ! where dTheta is the grid resolution at the poles.
    !/

    use ModSceGrids,     ONLY: IONO_nTheta, IONO_nPsi
    use ModSceVariables, ONLY: nThetaUsed, LatBoundary, HighLatBoundary, nThetaSolver, &
                               PhiIono_Weimer, nX, C_A, C_B, C_C, C_D, C_E, north, &
                               Radius, MaxIteration, UseWeimer, DoPrecond, USeInitialGuess, &
                               PhiOld_CB, cpcp_north, Tolerance, iHighBnd

    use nrtype, ONLY: DP, pio2_d, pi_d, cRadtoDeg

    use ModLinearSolver, ONLY: gmres, bicgstab, prehepta, Uhepta, Lhepta

    implicit none
 
    real(DP), intent(inout) :: Jr(:,:), SigmaThTh(:,:), SigmaThPs(:,:), SigmaPsPs(:,:), &
                               dSigmaThTh_dTheta(:,:), dSigmaThPs_dTheta(:,:), &
                               dSigmaPsPs_dTheta(:,:), dSigmaThTh_dPsi(:,:), &
                               dSigmaThPs_dPsi(:,:), dSigmaPsPs_dPsi(:,:), &
                               Theta(:,:), Psi(:,:), Phi_C(:,:), dTheta(:), dPsi(:)
    
    integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1
  
    ! Local variables
    integer  :: i, j, k, iMin, iMax, iI, nIteration, iError, iBlock
    real(DP) :: TermTheta2, TermTheta1, TermPsi2, TermPsi1, sn, cs, sn2, Residual, &
                PhiMax, PhiMin

    real(DP), allocatable :: lat_weimer(:,:), mlt_weimer(:,:), & 
                             dTheta2(:), dPsi2(:), SinTheta_I(:), CosTheta_I(:)

    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------
    allocate(lat_weimer(nTheta,nPsi), mlt_weimer(ntheta,nPsi), & 
             dTheta2(nTheta), dPsi2(nPsi), SinTheta_I(nTheta), CosTheta_I(nTheta))

    if (north) iBlock = 1
    if (.not. north) iBlock = 2

    ! Count the points above the latitude boundary
    nThetaUsed = count(abs(pio2_d-Theta(1:nTheta,1)) > LatBoundary)
    
    ! The pole is a boundary point
    iMin = nTheta - floor(maxval(HighLatBoundary))
    iMax = nThetaUsed+1

    do j=1, nPsi-1
       iHighBnd(j) = nTheta - floor(HighLatBoundary(j)) 
    end do
  
    nThetaSolver = iMax - iMin + 1
  
    ! get the weimer potential for the high polar cap region
    lat_weimer(1:nTheta,:) = 90. - Theta(1:nTheta,:)*cRadToDeg
    mlt_weimer(1:nTheta,:) = Psi(1:nTheta,:)*12.0/pi_d
    where(mlt_weimer <= 12)
       mlt_weimer = mlt_weimer + 12 
    elsewhere
       mlt_weimer = mlt_weimer - 12  
    end where

    ! symmetric for southern/northern hemisphere: lat_weimer is positive.
    call iono_potential_weimer(abs(lat_weimer(1:nTheta,1:nPsi-1)), &
                                   mlt_weimer(1:nTheta,1:nPsi-1), &
                                   PhiIono_weimer(1:nTheta,1:nPsi-1))
    PhiIono_weimer(:,nPsi) = PhiIono_weimer(:,1)
    write(*,*)'max min(weimer):',maxval(PhiIono_weimer), minval(PhiIono_weimer)

    ! Calculate (Delta Psi)^2 (note that dPsi = 2 Delta Psi)
    dPsi2 = (dPsi/2.0)**2
  
    ! Calculate (Delta Theta)^2  (dTheta = 2 Delta Theta for i=1 and nTheta)
    dTheta2        = (dTheta/2.0)**2
    dTheta2(1)     = 4*dTheta2(1)
    dTheta2(nTheta)= 4*dTheta2(nTheta)
  
    SinTheta_I = sin(Theta(:,1))
    CosTheta_I = cos(Theta(:,1))
  
    nX = nPsiUsed*nThetaSolver
 
    allocate(x(nX), y(nX), rhs(nX), b(nX), Bnd_I(nX), d_I(nX), e_I(nX), &
             e1_I(nX), f_I(nX), f1_I(nX))
  
    do j = 1, nPsiUsed
       do i= iMin, iMax
          ! set the matrix to identity for regions with Weimer potential
          if (north .and. i <= iHighBnd(j) .or. .not. north .and. i >= iHighBnd(j))then
             C_A(i,j) = 1.
             C_B(i,j) = 0.
             C_C(i,j) = 0.
             C_D(i,j) = 0.
             C_E(i,j) = 0.
          else
             sn  = SinTheta_I(i)
             cs  = CosTheta_I(i)
             sn2 = sn**2
             
             ! Central difference coefficients for second and first derivatives
             TermTheta2 = SigmaThTh(i,j)*sn2/dTheta2(i)
             TermTheta1 = (SigmaThTh(i,j)*sn*cs + dSigmaThTh_dTheta(i,j)*sn2 &
                         - dSigmaThPs_dPsi(i,j)*sn) / dTheta(i)
             TermPsi2 = SigmaPsPs(i,j) / dPsi2(j)            
             TermPsi1 = (dSigmaThPs_dTheta(i,j)*sn + dSigmaPsPs_dPsi(i,j)) / dPsi(j)
             
             ! Form the complete matrix
             C_A(i,j) = -2.0 * (TermTheta2 + TermPsi2)
             C_B(i,j) =         TermTheta2 - TermTheta1
             C_C(i,j) =         TermTheta2 + TermTheta1
             C_D(i,j) =         TermPsi2   - TermPsi1
             C_E(i,j) =         TermPsi2   + TermPsi1        
          end if
       end do
    enddo
  
    ! Fill in the diagonal vectors
    iI = 0
    do j = 1, nPsiUsed
       do i=iMin,iMax
          iI = iI + 1
          ! The right-hand-side is Jr * Radius^2 sin^2(Theta).
          b(iI)    = Jr(i,j)*(Radius*SinTheta_I(i))**2
          d_I(iI)  = C_A(i,j)
          e_I(iI)  = C_B(i,j)
          f_I(iI)  = C_C(i,j)
          e1_I(iI) = C_D(i,j)
          f1_I(iI) = C_E(i,j)
  
          if(i == iMin)     e_I(iI)  = 0.0
          if(i == iMax)     f_I(iI)  = 0.0
          if(j == 1)        e1_I(iI) = 0.0
          if(j == nPsiUsed) f1_I(iI) = 0.0
       end do
    end do
 
    ! Save original RHS for checking b = 0 !!! test for zero RHS
    Rhs = b
  
    ! move nonlinear part of operator to RHS
    x = 0.0
    UseWeimer = .True.
    DoPrecond = .False.
    call matvec_ionosphere(x, Bnd_I, nX)
    b = b - Bnd_I
    UseWeimer = .False.

    DoPrecond = .False.
    if(DoPrecond)then
       ! A -> LU
       call prehepta(nX,1,nThetaSolver,nX,real(-0.5,kind=8),d_I,e_I,f_I,e1_I,f1_I)
       ! rhs'=U^{-1}.L^{-1}.rhs
       call Lhepta(       nX,1,nThetaSolver,nX,b,d_I,e_I,e1_I)
       call Uhepta(.true.,nX,1,nThetaSolver,nX,b,f_I,f1_I)
    end if
 
    ! Solve A'.x = rhs'
    Residual = Tolerance
    if(UseInitialGuess) then
       iI = 0
       do j = 1, nPsiUsed
          do i=iMin,iMax
             iI = iI + 1
             x(iI) = PhiOld_CB(i,j,iBlock)
          end do
       end do
    else
       x = 0.0
    end if
  
    !select case(NameSolver)
    !case('gmres')
    !   nIteration = MaxIteration
    !   call gmres(matvec_ionosphere, b, x, UseInitialGuess, nX, MaxIteration, Residual, 'abs', nIteration, iError, DoTestMe)
    !case('bicgstab')
       nIteration = 3*MaxIteration
       call bicgstab(matvec_ionosphere, b, x, UseInitialGuess, nX, Residual, 'abs', nIteration, iError, DoTestMe)
    !case default
    !   call CON_stop(NameSub//': unknown NameSolver='//NameSolver)
    !end select
  
    ! Phi_C is the solution within the solved region
    Phi_C(1:iMin,:) = PhiIono_Weimer(1:iMin,:)
    iI = 0
    do j=1, nPsiUsed
       do i = iMin, iMax
          iI = iI + 1
          if (i>=iHighBnd(j)) Phi_C(i,j) = x(iI)     
       end do
    end do
    Phi_C(iMin,:) = PhiIono_Weimer(iMin,:)

    ! Apply periodic boundary condition in Psi direction
    Phi_C(:,nPsi) = Phi_C(:,1)
    PhiMax = maxval(Phi_C)
    PhiMin = minval(Phi_C)
  
    ! Save the solution for next time
    PhiOld_CB(:,:,iBlock) = Phi_C
    do j=1, nPsi-1
       PhiOld_CB(1:iHighBnd(j), j, iBlock) = PhiIono_Weimer(1:iHighBnd(j),j)
    end do
    ! apply periodic boundary condition in Psi direction
    PhiOld_CB(:,nPsi,iBlock) = PhiOld_CB(:,1,iBlock)
  
    ! Apply average condition at north pole
    !Phi_C(1,:) = sum(Phi_C(2,1:nPsiUsed))/nPsiUsed
    cpcp_north = (PhiMax - PhiMin)/1000.0

    deallocate(x, y, b, rhs,Bnd_I, d_I, e_I, f_I, e1_I, f1_I)
    deallocate(lat_weimer, mlt_weimer, dTheta2, dPsi2, SinTheta_I, CosTheta_I)
    return  
  end subroutine ionosphere_solver

!============================================================================
  subroutine matvec_ionosphere(x_I, y_I, n)
 
    use ModRamMain, ONLY: DP
    use ModSceGrids, ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables, ONLY: HighLatBoundary, nThetaSolver, PhiIono_Weimer, &
                               iHighBnd, C_A, C_B, C_C, C_D, C_E, UseWeimer, &
                               nThetaUsed, DoPrecond

    use ModLinearsolver, ONLY: Uhepta, Lhepta

    implicit none
  
    integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1
  
    ! Calculate y = A.x where A is the (pentadiagonal) matrix
  
    integer, intent(in) :: n          ! number of unknowns
    real(DP), intent(in) :: x_I(n)        ! vector of unknowns
    real(DP), intent(out):: y_I(n)        ! y = A.x
  
    integer :: iTheta, iTheta2, iPsi, i, iMin, iMax,j
    real(DP), allocatable :: x_G(:,:) ! 2D array with ghost cells
    !-------------------------------------------------------------------------
    allocate(x_G(0:nTheta+1, 0:nPsi))

    iMin = nTheta - floor(maxval(HighLatBoundary))
    iMax = nThetaUsed+1  
  
    nThetaSolver = iMax - iMin + 1
  
    x_G = 0.0
  
    ! Put 1D vector into 2D solution
    i = 0
    do iPsi = 1, nPsi-1
       do iTheta = iMin, iMax
          i = i+1
          x_G(iTheta, iPsi) = x_I(i)
       enddo
    enddo
  
    !if (UseWeimer) then ! apply the boundary condition at high latitude
    do j=1, nPsi-1
       x_G(iMin-1:iHighBnd(j), j) = PhiIono_Weimer(iMin-1:iHighBnd(j),j)
    end do
    !else
    !   do j=1, nPsi-1
    !      x_G(iMin-1:iHighBnd(j), j) = 0.0
    !   end do
    !end if
    x_G(iMax+1,1:nPsi-1) = 0.0
  
    ! Apply periodic boundary conditions in Psi direction
    x_G(:,nPsi) = x_G(:,1)
    x_G(:,0)    = x_G(:,nPsi-1)
  
    i = 0
    do iPsi = 1, nPsi-1
       do iTheta = iMin, iMax
          i = i+1     
          y_I(i) = C_A(iTheta, iPsi)*x_G(iTheta,   iPsi)   + &
                   C_B(iTheta, iPsi)*x_G(iTheta-1, iPsi)   + &
                   C_C(iTheta, iPsi)*x_G(iTheta+1, iPsi)   + &
                   C_D(iTheta, iPsi)*x_G(iTheta,   iPsi-1) + &
                   C_E(iTheta, iPsi)*x_G(iTheta,   iPsi+1)
       end do
    end do
  
    if (UseWeimer) then ! apply the boundary condition at high latitude
       i = 0
       do iPsi=1, nPsi-1
          do iTheta = iMin, iMax
             i = i + 1
             if (iTheta <= iHighBnd(iPsi)) y_I(i) = 0.0
          end do
       end do
    end if

    ! Preconditioning: y'= U^{-1}.L^{-1}.y
    if (DoPrecond) then
       call Lhepta(       n,1,nThetaSolver,n,y_I,d_I,e_I,e1_I)
       call Uhepta(.true.,n,1,nThetaSolver,n,y_I,    f_I,f1_I)
    end if

    deallocate(x_G)
    return

  end subroutine matvec_ionosphere

!========================================================================
  subroutine iono_potential_weimer(lat, mlt, PhiWeimer)
  
    use ModRamTiming,    ONLY: TimeRamNow 
    use ModRamParams,    ONLY: UseSWMFFIle, NameOmniFile
    use ModRamVariables, ONLY: Kp
    use ModRamConst,     ONLY: RE
    use ModScbVariables, ONLY: tilt
    use ModSceVariables, ONLY: IONO_NORTH_Theta, IONO_NORTH_Psi
    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi

    use w05

    use nrtype, ONLY: DP

    use ModTimeConvert, ONLY: n_day_of_year
    use ModIOUnit, ONLY: UNITTMP_

    implicit none
  
    real(DP), intent(in) :: lat(:,:), mlt(:,:)
    real(DP), intent(inout) :: PhiWeimer(:,:)

    integer :: i,j, maxN, minN
    integer :: doy, ierr, iYear_l, iDoy_l, iHour_l, iMin_l, iLines, isec_l, &
               imsec_l, imonth_l, iday_l, AL_l, SymH_l
    REAL(DP) :: radius, angle, bzimf_l, bndylat, byimf_l, pdyn_l, Nk_l, &
                Vk_l, bTot_l, bximf_l, vx_l, vy_l, vz_l, t_l, maxV, minV
    CHARACTER(LEN = 100) :: header
    !---------------------------------------------------------------------
    OPEN(UNITTMP_, FILE=NameOmniFile, status = 'UNKNOWN', action = 'READ')
    i = 0
    ierr = 0
    DO WHILE (ierr == 0)
       READ(UNITTMP_, *, IOSTAT=ierr) header
       i = i+1
    END DO
    iLines = i-1
    !C PRINT*, 'IP: the file has ', iLines, ' lines.'

    ! Rewind file
    REWIND(UNITTMP_)

    if (UseSWMFFile) then
       print*, 'IP: year, month, day, hour, min, by, bz, v, n'
       Read_SWMF_file_05: DO i = 1, iLines
          read(UNITTMP_,*) iyear_l, imonth_l, iday_l, ihour_l, imin_l, isec_l, imsec_l, &
                           bximf_l, byimf_l, bzimf_l, vx_l, vy_l, vz_l, nk_l, t_l
          if ((TimeRamNow%iYear.eq.iyear_l).and.(TimeRamNow%iMonth.eq.imonth_l).and. &
              (TimeRamNow%iDay.eq.iday_l).and. &
              (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
             vk_l   = SQRT(vx_l**2 + vy_l**2 + vz_l**2)
             write(*,*)iyear_l,imonth_l,iday_l, ihour_l, imin_l, byimf_l, bzimf_l, vk_l, nk_l
             EXIT Read_SWMF_file_05
          END IF
       END DO Read_SWMF_file_05
    else
       PRINT*, 'IP: Year, Day, Hour, Min, Bt, By, Bz, V, N, Pdyn, AL, SYMH'
       Read_OMNI_file_05: DO i = 1, iLines
          READ(UNITTMP_,*) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, &
                           bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
          doy = n_day_of_year(TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay)
          if ((TimeRamNow%iYear.eq.iyear_l).and.(doy.eq.idoy_l).and. &
              (TimeRamNow%iHour.eq.ihour_l).and.(TimeRamNow%iMinute.eq.imin_l)) then
             WRITE(*,*) iYear_l, iDoy_l, iHour_l, iMin_l, bTot_l, byimf_l, &
                        bzimf_l, Vk_l, Nk_l, pdyn_l, AL_l, SymH_l
             EXIT Read_OMNI_file_05
          END IF
       END DO Read_OMNI_file_05
    endif
    CLOSE(UNITTMP_)

    if (bzimf_l.lt.-20) bzimf_l = -20
    if (bzimf_l.gt.20) bzimf_l = 20
    CALL SetModel05(byimf_l, bzimf_l, tilt, Vk_l, Nk_l)

    ! calculate weimer potential from w05 for the high-latitude potential used in the solver
    DO j = 1, Iono_nPsi-1
       DO i = 1, Iono_nTheta
          ! use the weimer 2005 model
          call EpotVal05(lat(i,j), mlt(i,j), 0.0_dp, PhiWeimer(i,j))
          PhiWeimer(i,j) = PhiWeimer(i,j) * 1.0e3 ! in Volts
       END DO
       maxV = maxval(PhiWeimer(:,j))
       maxN = count(PhiWeimer(:,j).gt.0.0)
       minV = minval(PhiWeimer(:,j))
       minN = count(PhiWeimer(:,j).lt.0.0)
       ThetaLoop: DO i = 1, Iono_nTheta
          if (lat(i,j) > 75) cycle ThetaLoop
          if (maxN > minN) then
             if (PhiWeimer(i,j).le.0.1*maxV) PhiWeimer(i,j) = 0.1*maxV
          else
             if (PhiWeimer(i,j).ge.0.1*minV) PhiWeimer(i,j) = 0.1*minV
          endif
       ENDDO ThetaLoop
    END DO

  end subroutine iono_potential_weimer

!==================================================================================================
END MODULE ModSceIono

