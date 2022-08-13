MODULE ModSceIono

  use nrtype, ONLY: DP

  implicit none

  real(DP), allocatable :: x(:), y(:), rhs(:), b(:), Bnd_I(:), d_I(:), e_I(:), &
                           f_I(:), e1_I(:), f1_I(:)

  contains
!==================================================================================================
  subroutine ionosphere_fluxes(iModel)

    !\
    ! Combine the diffusive and discrete precipitating fluxes
    !/

    use ModSceGrids,     ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables, ONLY: IONO_Min_Ave_E, &
                               IONO_Min_EFlux, iono_north_im_eflux, iono_north_im_avee, &
                               iono_north_im_dis_eflux, iono_north_im_dis_avee, &
                               iono_north_im_eflux_diff, &
                               iono_north_ave_e, iono_north_eflux, iono_north_eflux_diff

    use nrtype, ONLY: DP

    implicit none
    
    integer, intent(in) :: iModel
    real(DP), allocatable :: discrete_ef(:,:), discrete_ae(:,:), &
                             diffuse_ef(:,:), diffuse_ae(:,:)
    !---------------------------------------------------------------------------
    ! iModel = 7: IM precipitation flux + Robinson's formula
    ! iModel = 9: IM precipitation flux + GLOW's calculation
     allocate(discrete_ef(Iono_nTheta,Iono_nPsi), discrete_ae(Iono_nTheta,Iono_nPsi), &
             diffuse_ef(Iono_nTheta,Iono_nPsi), diffuse_ae(Iono_nTheta,Iono_nPsi))

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

    where (diffuse_ae < IONO_Min_Ave_E/2) diffuse_ae = IONO_Min_Ave_E/2
    where (discrete_ae < IONO_Min_Ave_E/2) discrete_ae = IONO_Min_Ave_E/2

    ! Let's weight the average energy by the number flux, which is ef/av
    iono_north_ave_e = (diffuse_ef + discrete_ef) &
                      /(diffuse_ef/diffuse_ae + discrete_ef/discrete_ae)

    ! The energy flux should be weighted by the average energy
    iono_north_eflux = (diffuse_ef/diffuse_ae + discrete_ef/discrete_ae) * iono_north_ave_e
    where (iono_north_ave_e < IONO_Min_Ave_E) iono_north_ave_e = IONO_Min_Ave_E

    if (iModel .eq. 9)then  !!! only consider diffuse if coupled to GLOW.
       iono_north_eflux = diffuse_ef ! ergs/cm^2/s/1000. 
       iono_north_ave_e = diffuse_ae ! keV
       iono_north_eflux_diff = iono_north_im_eflux_diff ! differential electron flux (cm^2/s/sr/keV)
    end if
    
    deallocate(discrete_ef, discrete_ae, diffuse_ef, diffuse_ae)
    return
  end subroutine ionosphere_fluxes

!==================================================================================================
  subroutine ionosphere_conductance(Sigma0, SigmaH, SigmaP, SigmaThTh, SigmaThPs, &
                                    SigmaPsPs, dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                                    dSigmaPsPs_dTheta, dSigmaThTh_dPsi, dSigmaThPs_dPsi, &
                                    dSigmaPsPs_dPsi, Eflux, Ave_E, Eflux_diff, Theta, Psi, nTheta, &
                                    nPsi, dTheta, dPsi, x, y, z, iModel)
  
    !\
    ! This subroutine computes the height-integrated field-aligned and
    ! Hall and Pedersen conductances for the ionosphere at each
    ! location of the discretized solution domain.  The gradients of
    ! these quantities are also computed.
    !/

    use ModRamVariables, ONLY: F107, EKEV
    use ModRamMain,        ONLY: PathSceOut
    use ModRamGrids,       ONLY: nE
    use ModRamTiming,      ONLY: TimeRamElapsed, TimeRamNow

    use ModSceGrids,        ONLY: Iono_nTheta, Iono_nPsi
    use ModSceVariables,  ONLY: IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                               StarLightPedConductance, SAVE_NORTH_SigmaH, &
                               SAVE_NORTH_SigmaP, cosThetaTilt, sinThetaTilt, &
                               IONO_Radius, IONO_Height, nzGlow, DoUseFullSpec,&
                               DoSaveGLOWConductivity

    use nrtype, ONLY: DP, cDegtoRad, Pi_d
    
    use ModIOUnit,         ONLY: UnitTmp_
    use CON_axes,          ONLY: transform_matrix
    use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
    use ModRamMpi
    use ModMpi

    implicit none

    integer, intent(in)     :: iModel
    integer, intent(in)     :: nTheta, nPsi
    real(DP), intent(inout) :: Sigma0(:,:), SigmaH(:,:), SigmaP(:,:), SigmaThTh(:,:), &
                               SigmaThPs(:,:), SigmaPsPs(:,:), dSigmaThTh_dTheta(:,:), &
                               dSigmaThPs_dTheta(:,:), dSigmaPsPs_dTheta(:,:), &
                               dSigmaThTh_dPsi(:,:), dSigmaThPs_dPsi(:,:), &
                               dSigmaPsPs_dPsi(:,:), Eflux(:,:), Ave_E(:,:),  Eflux_diff(:,:,:), &
                               Theta(:,:), Psi(:,:), x(:,:), y(:,:), z(:,:), dTheta(:), dPsi(:)
    integer  :: i, j, k
    logical  :: old
    real(DP) :: f107p53, f107p49, cos_limit, meeting_value_p, meeting_value_h, &
                SigmaH_EUV, SigmaP_EUV, SigmaH_SCAT, SigmaP_SCAT, SigmaH_EUV_2, &
                SigmaP_EUV_2, SigmaH_STAR, SigmaP_STAR, sn, cs, sn2, cs2, cs3, &
                cs4, C
    
    real :: ut, ap,  f107r, f107a, f107p, f107y
    real(DP) :: rIono, XyzSmg(3), XyzGeo(3)
    integer  :: iyear, imonth, iday, ihour, iminute, isecond, idoy, ndaymo,ap0,ipoint, idate
    character(len=100) :: NameFile

    real(DP), allocatable :: cos_SZA(:,:)

    real(DP), allocatable :: SigmaH_Glow(:,:), SigmaP_Glow(:,:), &
                             SigmaH_Glow_all(:,:), SigmaP_Glow_all(:,:),&
                             SigmaH_all(:,:), SigmaP_all(:,:),&
                             SigmaH_Particles(:,:), SigmaP_Particles(:,:),&
                             glat(:,:), glong(:,:)
    real(DP), allocatable:: zz(:,:,:), ionrate(:,:,:), eDen(:,:,:),&
                            Pedcond(:,:,:), Hallcond(:,:,:), &
                            zz_all(:,:,:),ionrate_all(:,:,:), eDen_all(:,:,:),&
                            Pedcond_all(:,:,:), Hallcond_all(:,:,:)

    character(len=*), parameter :: NameSub='ionosphere_conductance'
    !-------------------------------------------------------------------------
    allocate(cos_SZA(nTheta,nPsi), SigmaH_Particles(1:IONO_nTheta,1:IONO_NPsi),    &
            SigmaP_Particles(1:IONO_nTheta,1:IONO_NPsi))
 
    cos_SZA = (x*cosTHETATilt-z*sinTHETATilt)/sqrt(x**2 + y**2 + z**2)
  
    ! We are going to need F10.7 ^ 0.53 and F10.7 ^ 0.49 a lot,
    ! So, let's just store them straight away:  
    f107p53 = f107**0.53
    f107p49 = f107**0.49
    cos_limit = cos(70.0*cDegToRad)
    meeting_value_p = f107p49*(0.34*cos_limit+0.93*sqrt(cos_limit))
    meeting_value_h = f107p53*(0.81*cos_limit+0.54*sqrt(cos_limit))

    if (iModel .eq. 9)then
       ! prepare parameters 
       iyear = TimeRamNow%iyear
       imonth= TimeRamNow%iMonth
       iday  = TimeRamNow%iDay
       iHour = TimeRamNow%iHour
       iMinute=TimeRamNow%iMinute
       iSecond=TimeRamNow%iSecond
       ut     = iHour * 3600. + iMinute * 60 + iSecond

       call moda(0, iyear, imonth, iday, idoy, ndaymo)
       idate = (iyear-iyear/100*100)*1000+idoy
       
       ! read in the parameters Ap, F107 (from irifun_2012.f) (ap0 in integer) 
       call apf_only(iyear, imonth, iday, f107r, f107p, f107a, f107y,ap0)
       ap = ap0*1.0

       allocate(zz(1:IONO_nTheta,1:IONO_nPsi,nzGlow),  &
            ionrate(1:IONO_nTheta,1:IONO_nPsi,nzGlow), &
            eDen(1:IONO_nTheta,1:IONO_nPsi,nzGlow),    &
            Pedcond(1:IONO_nTheta,1:IONO_nPsi,nzGlow), &
            Hallcond(1:IONO_nTheta,1:IONO_nPsi,nzGlow),&
            SigmaH_Glow(1:IONO_nTheta,1:IONO_NPsi),    &
            SigmaP_Glow(1:IONO_nTheta,1:IONO_NPsi),    &
            zz_all(1:IONO_nTheta,1:IONO_nPsi,nzGlow),  &
            ionrate_all(1:IONO_nTheta,1:IONO_nPsi,nzGlow), &
            eDen_all(1:IONO_nTheta,1:IONO_nPsi,nzGlow),    &
            Pedcond_all(1:IONO_nTheta,1:IONO_nPsi,nzGlow), &
            Hallcond_all(1:IONO_nTheta,1:IONO_nPsi,nzGlow),&
            SigmaH_Glow_all(1:IONO_nTheta,1:IONO_NPsi),    & 
            SigmaP_Glow_all(1:IONO_nTheta,1:IONO_NPsi),    &
            SigmaH_all(1:IONO_nTheta,1:IONO_NPsi),    & 
            SigmaP_all(1:IONO_nTheta,1:IONO_NPsi))
    end if

    
    iPoint = -1
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

          if (iModel .ne. 9)then ! not use GLOW
   
             !\
             ! Use Robinson's Formula to convert the Ave_E and E_Flux to SigmaP and SigmaH
             !/
             SigmaP_Particles(i,j) = 40.0 * Ave_E(i,j) / (16.0 + Ave_E(i,j)*Ave_E(i,j))  &
                  * sqrt(EFlux(i,j)*1000.0)
  
             SigmaH_Particles(i,j) = 0.45 * (Ave_E(i,j)**0.85) * SigmaP_Particles(i,j)
             
             SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
                  SigmaH_SCAT*SigmaH_SCAT + &
                  SigmaH_STAR*SigmaH_STAR + &
                  SigmaH_Particles(i,j)*SigmaH_Particles(i,j))
             
             SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV + &
                  SigmaP_SCAT*SigmaP_SCAT + &
                  SigmaP_STAR*SigmaP_STAR + &
                  + SigmaP_Particles(i,j)*SigmaP_Particles(i,j))

           else
             !\
             ! use GLOW model for the conductance calculation including EUV etc.
             ! For auroral: either with the single EFlux&Ave_E to assume Maxwellian or with
             ! the FullSpetrum of the Eflux_differential. 
             ! Pass the following: time, location, ap, f107, ef, ev 
             !/ 
             rIono = (IONO_Radius + IONO_Height)/IONO_Radius
             call sph_to_xyz(rIono, Theta(i,j), Psi(i,j), xyzSmg)
             XyzGeo = matmul(transform_matrix(TimeRamElapsed, 'SMG','GEO'), XyzSmg)
             call xyz_to_sph(XyzGeo, rIono,glat(i,j),glong(i,j)) ! (r, theta, phi) 

             glong(i,j) = glong(i,j) * 180./Pi_d
             glat(i,j) = 90 - glat(i,j) * 180./Pi_d !(north/south the same form) 
             
             !\ 
             ! parallize the glow calculation for the 2D points over ionosphere
             !/ 
             if(nProc>1)then
                iPoint = iPoint + 1
                if (mod(iPoint, nProc) /=iProc)CYCLE
             end if
             !\
             ! pass the energy flux and characteristic energy (half of the mean energy)
             ! EFlux*1000 (ergs/cm^2/s); EFlux_Diff(/cm^2/s/sr/keV); Ave_E: keV 
             !/
             if (EFlux(i,j)*1000. > 0.0001 .and. Ave_E(i,j)/2.*1000. >1)then
                call glow_aurora_conductance(idate, real(ut,DP), glat(i,j), glong(i,j), &
                     SigmaP_Glow(i,j), SigmaH_Glow(i,j), &
                     EFlux(i,j)*1000., Ave_E(i,j)/2., EFlux_Diff(i,j,:), EkeV(:), nE, &
                     real(ap,kind=8), real(f107r,DP), real(f107p,DP),&
                     real(f107a,kind=8), DoUseFullSpec, &
                     zz(i,j,1:nzGlow), ionrate(i,j,1:nzGlow), eDen(i,j,1:nzGlow), &
                     Pedcond(i,j,1:nzGlow), Hallcond(i,j,1:nzGlow), nzGlow)
             else
                SigmaP_Glow(i,j) = 0.0
                SigmaH_Glow(i,j) = 0.0
             end if

             ! If the glow model provides the solar flux & auroral flux related conductance.                                                 
             ! SigmaH(i,j) = SigmaH_Glow(i,j)                                                                                             
             ! SigmaP(i,j) = SigmaP_Glow(i,j)                                                                                             
             ! else add the EUV-conductance here as follows.
             
             SigmaH(i,j) = sqrt(SigmaH_EUV*SigmaH_EUV + &
                  SigmaH_SCAT*SigmaH_SCAT + &
                  SigmaH_STAR*SigmaH_STAR + &
                  SigmaH_Glow(i,j)*SigmaH_Glow(i,j))
             
             SigmaP(i,j) = sqrt(SigmaP_EUV*SigmaP_EUV + &
                  SigmaP_SCAT*SigmaP_SCAT + &
                  SigmaP_STAR*SigmaP_STAR + &
                  SigmaP_Glow(i,j)*SigmaP_Glow(i,j))
             
          end if

       enddo
    enddo


    if (imodel .eq. 9 .and. nProc > 1)then
       ! MPI reduce to the head node (iproc=0)                                                                                           
       
       call MPI_reduce(SigmaP, SigmaP_all, nTheta*nPsi, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(SigmaH, SigmaH_all, nTheta*nPsi, MPI_REAL, MPI_SUM, 0, iComm, iError)
       
       call MPI_reduce(SigmaP_Glow, SigmaP_Glow_all, nTheta*nPsi, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(SigmaH_Glow, SigmaH_Glow_all, nTheta*nPsi, MPI_REAL, MPI_SUM, 0, iComm, iError)
       
       call MPI_reduce(zz,      zz_all,      nTheta*nPsi*nzGlow, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(ionrate, ionrate_all, nTheta*nPsi*nzGlow, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(eDen,    eDen_all,    nTheta*nPsi*nzGlow, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(Pedcond, Pedcond_all, nTheta*nPsi*nzGlow, MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(Hallcond,Hallcond_all,nTheta*nPsi*nzGlow, MPI_REAL, MPI_SUM, 0, iComm, iError)

       if(iProc==0)then
          
          SigmaP = SigmaP_all
          SigmaH = SigmaH_all
          
          SigmaP_Glow = SigmaP_Glow_all
          SigmaH_Glow = SigmaH_Glow_all
          
          zz = zz_all
          ionrate = ionrate_all
          eDen = eDen_all
          Pedcond = Pedcond_all
          Hallcond= Hallcond_all
       end if
       !bcast to other processors
       call MPI_bcast(SigmaP, nTheta*nPsi,MPI_REAL,0,iComm,iError)
       call MPI_bcast(SigmaH, nTheta*nPsi,MPI_REAL,0,iComm,iError)

       !! write out the conductance into file                                                                                               
       if (iProc==0 .and. DoSaveGLOWConductivity)then
          if (mod(TimeRamElapsed, 300.0) .eq. 0)then
             write(namefile, '(a,i6.6,a)')PathSceOut//"Conductance_",nint(TimeRamElapsed/300.),".dat"
             open(UnitTmp_, file=trim(namefile),status='unknown')
             write(UnitTmp_, '(a, i3.3)')'nHeight: ', nzGlow
             write(UnitTmp_, '(a, i4.4,1x,i2.2,1x, i2.2,1x,i2.2,1x,i2.2))')'Time: ', iyear, imonth,&
                  iday, ihour, iminute
             write(unitTmp_, '(a)')'Theta Psi Glat Glon EFlux Emean SigamP_all SigmaH_all SigmaP_Glow SigamH_Glow SigmaP_R SigmaH_R'
             write(unitTmp_, '(a)')'zz ionization_rate Ne Pedconductivity Halconductivity'
             
             do j=1, nPsi
                do i=1, nTheta
                   
                   write(UnitTmp_,'(1x, 10f9.4)')&
                        Theta(i,j), Psi(i,j), glat(i,j), glong(i,j), EFlux(i,j)*1000, Ave_E(i,j), &
                        SigmaP(i,j),SigmaH(i,j), SigmaP_Glow(i,j), SigmaH_Glow(i,j)
                   
                   do k=1,nzGlow
                      write(UnitTmp_,'(f6.2, 4e12.3)')zz(i,j,k), ionrate(i,j,k), eDen(i,j,k), &
                           Pedcond(i,j,k), Hallcond(i,j,k)
                   end do
                end do
             end do
             close(UnitTmp_)
          end if
       end if
    end if
    
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
  
    if(iModel .eq. 9)deallocate(zz,ionrate,eDen,PedCond,HallCond,SigmaH_GLOW,SigmaP_GLOW,&
         zz_all,ionrate_all,eDen_all,PedCond_all,HallCond_all,SigmaH_GLOW_all,SigmaP_GLOW_all,&
         SigmaH_all,SigmaP_all)
    
    deallocate(cos_SZA,SigmaH_Particles,SigmaP_Particles)
    return

  end subroutine ionosphere_conductance

!==================================================================================================
  subroutine ionosphere_solver(Jr, SigmaThTh, SigmaThPs, SigmaPsPs, &
                               dSigmaThTh_dTheta, dSigmaThPs_dTheta, &
                               dSigmaPsPs_dTheta, dSigmaThTh_dPsi, &
                               dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                               Theta, Psi, dTheta, dPsi, Phi_C)
  
    ! Modified from SWMF/Ridley_serials IE solver. By Yu & Toth, 2016.
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
                               Radius, MaxIteration, UseWeimer, DoPrecond, UsePreconditioner,USeInitialGuess, &
                               PhiOld_CB, cpcp, Tolerance, iHighBnd, NameSolver

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
    integer  :: i, j, k, iMin, iMax, iI, nIteration, iError
    real(DP) :: TermTheta2, TermTheta1, TermPsi2, TermPsi1, sn, cs, sn2, Residual, &
                PhiMax, PhiMin

    real(DP), allocatable :: lat_weimer(:,:), mlt_weimer(:,:), & 
                             dTheta2(:), dPsi2(:), SinTheta_I(:), CosTheta_I(:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'ionosphere_solver'
    !-------------------------------------------------------------------------
    allocate(lat_weimer(nTheta,nPsi), mlt_weimer(ntheta,nPsi), & 
             dTheta2(nTheta), dPsi2(nPsi), SinTheta_I(nTheta), CosTheta_I(nTheta))

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

    DoPrecond = UsePreconditioner
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
             x(iI) = PhiOld_CB(i,j)
          end do
       end do
    else
       x = 0.0
    end if
  
    select case(NameSolver)
    case('gmres')
       nIteration = MaxIteration
       call gmres(matvec_ionosphere, b, x, UseInitialGuess, nX, MaxIteration, Residual, 'abs', nIteration, iError, DoTestMe)
    case('bicgstab')
       nIteration = 3*MaxIteration
       call bicgstab(matvec_ionosphere, b, x, UseInitialGuess, nX, Residual, 'abs', nIteration, iError, DoTestMe)
    case default
       call CON_stop(NameSub//': unknown NameSolver='//NameSolver)
    end select
  
    ! Phi_C is the solution within the solved region
    Phi_C(1:iMin,:) = 0.0 !PhiIono_Weimer(1:iMin,:)
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
    PhiOld_CB(:,:) = Phi_C
    do j=1, nPsi-1
       PhiOld_CB(1:iHighBnd(j), j) = PhiIono_Weimer(1:iHighBnd(j),j)
    end do
    ! apply periodic boundary condition in Psi direction
    PhiOld_CB(:,nPsi) = PhiOld_CB(:,1)
  
    ! Apply average condition at north pole
    !Phi_C(1,:) = sum(Phi_C(2,1:nPsiUsed))/nPsiUsed
    cpcp = (PhiMax - PhiMin)/1000.0

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
  
    if (UseWeimer) then ! apply the boundary condition at high latitude
       do j=1, nPsi-1
          x_G(iMin-1:iHighBnd(j), j) = PhiIono_Weimer(iMin-1:iHighBnd(j),j)
       end do
    else
       do j=1, nPsi-1
          x_G(iMin-1:iHighBnd(j), j) = 0.0
       end do
    end if
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
subroutine glow_aurora_conductance(idate, ut, glat, glong, &
     SigmaP, SigmaH, ef_in, ec_in, ef_diff,ekeV_diff, nE,  &
     ap,f107,f107p, f107a, flux_spec, z, ionrate, eDen,      &
     Pedcond, Hallcond, nz)

  use ModGlowBasic, ONLY: glowbasic_ram
  use nrtype,       ONLY: pi_d

  implicit none

  real(DP),    intent(in) :: glat, glong, ef_in, ec_in, &
                             f107,f107p,f107a, ut, ap
  integer,     intent(in) :: idate, nE, nz
  logical,     intent(in) :: flux_spec
  real(DP),    intent(in) :: ef_diff(nE), ekev_diff(nE)
  real(DP),    intent(out):: SigmaP, SigmaH
  real(DP), dimension(nz), intent(out) :: z, ionrate, eDen, Pedcond, Hallcond
  real(DP) :: ef, ec
  real(DP) :: logef_diff(nE), logec_diff(nE), ef_diff_tmp(nE)
  ! ---------------------------------------------------------------------                                                         
  ! calculate the conductance from glow model                                                                                     
  ef    = ef_in       ! ergs/cm^2/s                                                                                               
  ec    = ec_in*1.0e3 ! convert to eV                                                                                             

  if (flux_spec .eqv. .true.)then
     ef_diff_tmp = ef_diff

     if (abs(glat - 55.026)< 0.1 .and. abs(glong - 344.779) < 0.1)then
        write(*,*)'glat:',glat, 'glong:',glong
        write(*,*)'ef_diff:',ef_diff
     end if
                                                              ! Glow asks for /cm^2/s/eV for the flux                             
     logef_diff = log10(ef_diff_tmp/1000.*pi_d)             ! convert /cm^2/s/sr/keV to /cm^2/s/eV (integrate over pitch-angle)    
     logec_diff = log10(ekeV_diff*1.0e3)                   ! convert keV to eV                                                    
  end if

  if (flux_spec .eqv. .false.)then
     ! pass only flux level and characteristic energy                                                                             
     call glowbasic_ram(idate,ut, glat, glong, ap, f107, f107p, f107a, &
          ef, ec, SigmaP, SigmaH, nE, z, ionrate, eDen, Pedcond, Hallcond, nz)
  else
     ! pass the full flux spectra                                                                                                 
     call glowbasic_ram(idate, ut, glat, glong, ap, f107, f107p, f107a, &
          ef, ec, SigmaP, SigmaH, nE, z, ionrate, eDen, Pedcond, Hallcond, nz, &
          logec_diff, logef_diff)
  end if

end subroutine glow_aurora_conductance

!==================================================================================================
END MODULE ModSceIono

