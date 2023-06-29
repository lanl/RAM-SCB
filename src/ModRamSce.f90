!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamSce

  implicit none

  contains
!==================================================================================================
  subroutine calculate_precip_flux_jr(IS, nTheta, nPhi, rIono, energy_flux, ave_e, &
                                      num_flux, dis_energy_flux, dis_ave_e, &
                                      Jr, high_latboundary)

    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: nR, nT, nE, nPa
    use ModRamVariables, ONLY: WMU, MU, UPA, EKEV, XNE, WE, LZ, FLUX, PPERE, PPARE, &
                               DL1
    use ModScbGrids,     ONLY: nthe, npsi, nzeta, nXRaw, nYRaw, nXRawExt
    use ModScbVariables, ONLY: x, y, z, paraj, nThetaEquator, r0Start

    use ModRamGSL,       ONLY: GSL_Interpolation_2D
    use ModScbFunctions, ONLY: Extap

    use nrtype, ONLY: cElectronMass, cElectronCharge, pi_d, pio2_d, cRadtoDeg

    implicit none

    integer, intent(in)  :: IS, nTheta, nPhi
    real(DP), intent(in) :: rIono

    real(DP), intent(inout) :: energy_flux(:,:), ave_e(:,:), num_flux(:,:), Jr(:,:), &
                               dis_energy_flux(:,:), dis_ave_e(:,:), high_latboundary(:)

    integer  :: i, j, k, kk, l, nTheta_north, GSLerr, nn, iMin, d
    real(DP) :: rr1, thangle, thangleOnIono, minl, dTheta, dPhi,Rm, eV, efactor, press, &
                dydummy, radius, angle, beta=0.1, Nele, jpar, coordPoint(2)

    integer,  allocatable :: idx(:)
    real(DP), allocatable :: ave_flux(:,:,:), f(:,:,:,:), &
                             ave_fluxExt(:,:,:), num_fluxeq(:,:), NeExt(:,:), &
                             PPerEExt(:,:), PParEExt(:,:), rRawExt(:), aRawExt(:), &
                             colatGrid(:,:), lonGrid(:,:), F0(:,:), energy_flux_iono(:,:), &
                             num_flux_iono(:,:), ave_e_iono(:,:), NeEq(:,:), PParEEq(:,:), &
                             PPerEEq(:,:), dis_num_flux_iono(:,:), dis_ave_e_iono(:,:), &
                             dis_energy_flux_iono(:,:), ave_fluxEq(:,:,:), xTemp(:,:), &
                             yTemp(:,:)

    real(DP), allocatable :: ave_fluxtmp(:,:), ave_fluxEQtmp(:,:), energy_fluxtmp(:,:), &
                             ave_etmp(:,:), num_fluxtmp(:,:), jr_tmp(:,:), &
                             dis_energy_fluxtmp(:,:), dis_ave_etmp(:,:), &
                             dis_num_fluxtmp(:,:), high_latboundary_tmp(:), &
                             colat(:), lon(:), temp(:), rGrid(:,:), aGrid(:,:)
    !---------------------------------------------------------------------------
    !\
    ! by Yiqun Yu @ 2015
    ! Calculate energy flux and average energy from the precipitation number 
    ! flux at the ionosphere altitude.
    !/
    ! -- 1:  map the precipitation flux at the equator down to the ionosphere 
    !     altitude. Interpolate into the ionospheric grids.
    ! -- 2:  calculate the energy flux and ave. energy.
    !            <E> = \int{f(E)EdE}/\int{f(E)dE}
    !            jE  = \int{f(E)EdE}
    !            Here, f(E) is the energy differential flux, already integrated 
    !                   over the pitch angle.
    ! This is used to pass to IE for the computation of height-integrated Hall 
    ! and Pedeson Conductances, based on the Robinson 1987 formula. 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !\ Y. Yu, 2013, to couple with IE (Region2 Jr)
    ! interpolate the paraj onto the ionosphere grids
    ! like in the swmfiono_II <-- PhiIono(1:npsi,2:zeta) <-- PhiIono(nIeTheta,
    ! nIePhi)
    !/
    ! paraJr is positive when along the field line into the north pole (this is 
    ! opposite to that in GM)
    !----------------------------------------------------------------------------
    nTheta_north = int(nTheta/2)+1

    ALLOCATE(idX(nPhi))
    ALLOCATE(ave_flux(nR,nT,nE), f(nR,nT,NE,NPA), &
             ave_fluxext(nXRawExt,nYRaw,nE), num_fluxeq(nR,nT), NeExt(nXRawExt,nYRaw), &
             PperEExt(nXRawExt,nYRaw), PparEExt(nXRawExt,nYRaw), rRawExt(nXRawExt), aRawExt(nYRaw), &
             colatGrid(npsi,nzeta), lonGrid(npsi,nzeta), F0(npsi,nzeta), &
             energy_flux_iono(npsi,nzeta), num_flux_iono(npsi,nzeta), ave_e_iono(npsi,nzeta), &
             NeEq(npsi,nzeta), PParEEq(npsi,nzeta), PPerEEq(npsi,nzeta), &
             dis_num_flux_iono(npsi,nzeta), dis_ave_e_iono(npsi,nzeta), &
             dis_energy_flux_iono(npsi,nzeta), ave_fluxEq(npsi,nzeta,nE), &
             colat(nTheta_north), lon(nPhi), xTemp(npsi,nzeta), yTemp(npsi,nzeta))
     allocate(Energy_fluxtmp(nTheta,nPhi), Ave_etmp(nTheta,nPhi), Num_fluxtmp(nTheta,nPhi), &
              dis_Energy_fluxtmp(nTheta,nPhi), dis_Ave_etmp(nTheta,nPhi), Jr_tmp(nTheta,nPhi), &
              high_latboundary_tmp(1:nPhi/2+1), temp(nXRawExt), rGrid(npsi,nzeta), &
              aGrid(npsi,nzeta))

    energy_flux = 0.0
    ave_e = 0.0
    num_flux = 0.0
    Jr = 0.0
    high_latboundary = 0.5*pi_d

    ! FLUX is saved from ram_all for each species before F2 is evolved.
    F(2:NR,1:nT-1,2:NE,2:NPA) = FLUX(IS,2:NR,1:NT-1,2:NE,2:NPA)
    F(2:NR,nT,2:NE,2:NPA) = F(2:NR,1,2:NE,2:NPA)
    F(2:NR,:,2:NE,1)      = F(2:NR,:,2:NE,2)

    ! calculate the averaged flux at the equator for the precipitation flux (Jordanova 1997).
    ! averaged over the pitch angle (in the loss cone), so to remove the angle dependence.
    ! then for the low altitude flux, it has no pitch angle dependence (i.e., isotropic now).
    do i=2, nR
       do j=1, nT
          num_fluxeq(i,j) = 0.0
          do k=2,nE
             ave_flux(i,j,k) = 0.0
             do l=upa(i), npa
                ave_flux(i,j,k) = ave_flux(i,j,k)+f(i,j,k,l)*WMU(l)
             end do
             ave_flux(i,j,k) = ave_flux(i,j,k)/(mu(NPA)-mu(UPA(i)))
             num_fluxeq(i,j) = num_fluxeq(i,j) + pi_d*ave_flux(i,j,k)*WE(IS,k)
          end do
       end do
    end do

    ! Interpolate from RAM -> SCB
    ! RAM grid + extension
    rRawExt(1:nXRaw+1) = LZ(1:nxRaw+1)
    do j = nXRaw+2, nXRawExt
       rRawExt(j) = rRawExt(j-1) + DL1
    enddo
    DO k = 1, nYRaw ! starting from midnight (MLT=0)
       aRawExt(k) = 24.* REAL(k-1)/REAL(nYRaw-1)
    END DO
    aRawExt = aRawExt * 360./24 * pi_d / 180._dp ! Convert to Radians

    ave_fluxExt(1:nXRaw+1,:,:) = ave_flux(1:nXRaw+1,:,:)
    NeExt(1:nXRaw+1,:) = XNE(1:nXRaw+1,:)
    PPerEExt(1:nXRaw+1,:) = PPerE(1:nXRaw+1,:)
    PParEExt(1:nXRaw+1,:) = PParE(1:nXRaw+1,:)
    ! Extrapolation (for now just set equal)
    do i = nXRaw+2, nXRawExt
       ! Extrapolate average flux
       ave_fluxExt(i,:,:) = ave_fluxExt(i-1,:,:)
       ! Extrapolate electron density
       NeExt(i,:) = NeExt(i-1,:)
       ! Extrapolate perpendicular pressure
       PPerEExt(i,:) = PPerEExt(i-1,:)
       ! Extrapolate parallel pressure
       PParEExt(i,:) = PParEExt(i-1,:)
    end do
    where(ave_fluxExt .le. 1.0e-31) ave_fluxExt = 1.0e-31

    ! SCB Grid
    DO k = 2, nzeta
       DO j = 1, npsi
          radius = SQRT((x(nThetaEquator,j,k))**2 + y(nThetaEquator,j,k)**2)
          angle = ASIN(y(nThetaEquator,j,k) / radius) + pi_d
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .GE.0)) &
               angle = 2*pi_d - ASIN(y(nThetaEquator,j,k) / radius)
          IF ((x(nThetaEquator,j,k) .LE. 0) .AND. (y(nThetaEquator,j,k) .LE.0)) &
               angle = - ASIN(y(nThetaEquator,j,k) / radius)
          rGrid(j,k) = radius
          aGrid(j,k) = angle ! 0 degree at midnight, index starting from noon.
       END DO
    END DO

    CALL GSL_Interpolation_2D(rRawExt, aRawExt, NeExt(1:nXRawExt,1:nYRaw), &
                              rGrid(1:npsi,2:nzeta), aGrid(1:npsi,2:nzeta), NeEQ(1:npsi,2:nzeta), &
                              GSLerr)
    where(NeEQ.le.0.0) NeEQ = 0.0

    CALL GSL_Interpolation_2D(rRawExt, aRawExt, PPerEExt(1:nXRawExt,1:nYRaw), &
                              rGrid(1:npsi,2:nzeta), aGrid(1:npsi,2:nzeta), PPerEEQ(1:npsi,2:nzeta), &
                              GSLerr)
    where(PperEEQ.le.0.0) PperEEQ = 0.0

    CALL GSL_Interpolation_2D(rRawExt(1:nXRawExt), aRawExt(1:nYRaw), PParEExt(1:nXRawExt,1:nYRaw), &
                              rGrid(1:npsi,2:nzeta), aGrid(1:npsi,2:nzeta), PParEEQ(1:npsi,2:nzeta), &
                              GSLerr)
    where(PparEEQ.le.0.0) PparEEQ = 0.0

    DO k = 1, nE
       CALL GSL_Interpolation_2D(rRawExt(1:nXRawExt), aRawExt(1:nYRaw), ave_fluxExt(1:nXRawExt,1:nYRaw,k), &
                                 rGrid(1:npsi,2:nzeta), aGrid(1:npsi,2:nzeta), ave_fluxEQ(1:npsi,2:nzeta,k), &
                                 GSLerr)
       where(ave_fluxEQ(:,:,k).le.0.0) ave_fluxEQ(:,:,k) = 0.0
    END DO

    dis_num_flux_iono = 0.0
    dis_energy_flux_iono = 0.0
    dis_ave_e_iono = 0.0
    ! calculate the discreate number flux and energy flux at the equator in the scb grids
    do i=1,npsi
       do j=2,nzeta
          ! (electron) ppar in normalized unit, so *pnormal to nPa
          ! NeEQ in 1/cm^3
          ! bj in muA/m^2
          ! at the low altitude boundary
          press = (2*PperEEQ(i,j)+PparEEQ(i,j))/3.0*0.16*1.0e-9 !keV/cm^3-->nPa--> Pa
          Nele = NeEQ(i,j)*1.0e6 !1/m^3 at equator...

          ! bj(1,i,j) at the pole; bj(ntheequator,i,j) at the equator (!=0?)
          ! !units of paraj is muA/m^2
          jpar = abs(paraj(i,j))*1.0e-6 ! A/m^2 (paraj at the boundary)

          ! following Raeder 2001 (use parameter at ionosphere altitude
          if (paraj(i,j) .lt. 0 .and. press .gt. 0)then ! upward FACs region in the northern
             eV = cElectronCharge**2*Nele**1.5/sqrt(2*pi_d*cElectronMass*press)*jpar
             dis_energy_flux_iono(i,j) = 4*eV*jpar/1.6e-12         ! J/m^2/s--> keV/cm^2/s
             dis_ave_e_iono(i,j) = cElectronCharge*eV/1.6e-16      ! J-->keV
          end if
       end do
    end do

    !\
    ! Calculate the energy flux and averaged energy at 200km (along the B):
    !  -- precip_flux(at 200km) at all directions is equal to ave_flux (at equator).
    !  -- (the ave_flux is the same for all the pitch angles; 
    !  -- precip_flux(90deg) = ave_flux(at edge of the loss cone)
    !  -- isotropic now across the 200km atm. plane (same at r0Start, along B line)
    !  -- convert to the flux perpendicular to thep plane: flux*cos(alpha). 
    !  -- integrate over solid anle (--> *pi) and energy (--> energy flux)
    !/
    do i=1, npsi
       do j=2, nzeta
          energy_flux_iono(i,j) = 0.0
          num_flux_iono(i,j) = 0.0
          ave_e_iono(i,j) = 0.0
          do k=1, nE
             energy_flux_iono(i,j) = energy_flux_iono(i,j) + pi_d*ave_fluxEQ(i,j,k)*EKEV(IS,k)*WE(IS,k)
             num_flux_iono(i,j) = num_flux_iono(i,j) + pi_d*ave_fluxEQ(i,j,k)*WE(IS,k)
          end do
          if (num_flux_iono(i,j) .eq. 0.0) then
             num_flux_iono(i,j) = 1.0e-31
             ave_e_iono(i,j) = 1.0e-31
          else
             ave_e_iono(i,j) = energy_flux_iono(i,j)/num_flux_iono(i,j)
          end if
       end do
    end do

    ! the longitude is the same as in angleGrid (index start from the noon, 0 degree at midnight)
    lonGrid(:,2:nzeta) = aGrid(:,2:nzeta)

    ! Now this is mapped along the B field lines to the ionosphere altitude (~200km as the loss cone is calculated there)
    ! assume around the Earth surface (r0Start): the first grid point in 3D equli. code?
    ! the northest point, r0Start, not at the ionosphere altitude. So need to map to IonoAltitude 
    DO k = 2, nzeta
       DO j = 1, npsi
          rr1 = sqrt(x(nthe,j,k)**2+y(nthe,j,k)**2+z(nthe,j,k)**2)
          thangle = asin(z(nthe,j,k)/rr1) ! latitude at the northest point
          thangleOnIono = acos(sqrt( (cos(thangle))**2 * rIono/rr1))!r0Start )) ! latitude at IonoAltitude
          colatGrid(j,k) = 0.5*pi_d - thangle!OnIono
       END DO
    END DO

    ! now interpolate the scb spatial grid into the ionospheric grids.
    ! (nR, nT) --> (colatgrid, longrid)--> (colat, lon)
    dPhi = 2.0*pi_d/real(nPhi-1)
    dTheta = pi_d/real(nTheta-1)

    ! only find the north hemsiphere grid for interpolation
    colat = 0.0
    lon   = 0.0
    do i=2, nPhi ! Longitude goes from 0 to 360, index start from midnight, lon=0 is at midnight
       lon(i) = lon(i-1) + dPhi
    end do
    do i=2, nTheta_north ! Colat goes from 0 to 90.
       colat(i) = colat(i-1) + dTheta ! all the hemisphere
    end do

    xTemp = rIono*sin(colatGrid(:,:)) * cos(lonGrid(:,:))
    yTemp = rIono*sin(colatGrid(:,:)) * sin(lonGrid(:,:))
    do j=1, nPhi      
       ! find the closest longitude
       idx(j) = 2
       minl = abs(lon(j) - lonGrid(1,2))
       do k=2,nzeta
          if ( abs(lon(j) - lonGrid(1,k)) .lt. minl)then
             minl = abs(lon(j) - lonGrid(1,k))
             idx(j) = k
          end if
       end do
       high_latboundary(j) = 0.5*pi_d
       iMin = 0
       do i=1, nTheta_north
          if ((colat(i).lt.minval(colatGrid(:,idx(j)))).or. &
              (colat(i).gt.maxval(colatGrid(:,idx(j))))) then
             ! inside the polar inner or outer boundary of scb don't do the interpolation
             energy_flux(i,j)     = 0.0_dp
             ave_e(i,j)           = 0.0_dp
             num_flux(i,j)        = 0.0_dp
             dis_energy_flux(i,j) = 0.0_dp
             dis_ave_e(i,j)       = 0.0_dp
             Jr(i,j)              = 0.0_dp
             if (colat(i) .lt. minval(colatGrid(:,idx(j))))then
                if (high_latboundary(j) .gt. (0.5*pi_d - colat(i))) then
                   high_latboundary(j) = 0.5*pi_d-colat(i)
                end if
                iMin = i
             end if
          else
             coordPoint(1) = rIono*sin(colat(i)) * cos(lon(j))
             coordPoint(2) = rIono*sin(colat(i)) * sin(lon(j))
             call GSL_Interpolation_2D(xTemp, yTemp, energy_flux_iono, &
                                       coordPoint(1), coordPoint(2), energy_flux(i,j), &
                                       GSLerr)
             call GSL_Interpolation_2D(xTemp, yTemp, ave_e_iono, &
                                       coordPoint(1), coordPoint(2), ave_e(i,j), &
                                       GSLerr)
             call GSL_Interpolation_2D(xTemp, yTemp, num_flux_iono, &
                                       coordPoint(1), coordPoint(2), num_flux(i,j), &
                                       GSLerr)
             call GSL_Interpolation_2D(xTemp, yTemp, dis_energy_flux_iono, &
                                       coordPoint(1), coordPoint(2), dis_energy_flux(i,j), &
                                       GSLerr)
             call GSL_Interpolation_2D(xTemp, yTemp, dis_ave_e_iono, &
                                       coordPoint(1), coordPoint(2), dis_ave_e(i,j), &
                                       GSLerr)
             call GSL_Interpolation_2D(xTemp, yTemp, paraj, &
                                       coordPoint(1), coordPoint(2), Jr(i,j), &
                                       GSLerr)
             !do kk = 1, nE
             !   call GSL_Interpolation_2D(xTemp, yTemp, ave_fluxEQ(:,:,kk), &
             !                             coordPoint(1), coordPoint(2), diff_flux_iono(i,j,kk), &
             !                             GSLerr)
             !enddo
          end if
       end do
       Jr(iMin+1,j) = 0.0_dp
       energy_flux(iMin+1,j) = 0.0_dp
       ave_e(iMin+1,j) = 0.0_dp
       num_flux(iMin+1,j) = 0.0_dp
       dis_energy_flux(iMin+1,j) = 0.0_dp
       dis_ave_e(iMin+1,j) = 0.0_dp
    end do
    ! Jr is positive into the north pole, opposite to that in GM; to be consistent: -1*Jr
    Jr(1:nTheta_north,:) = -Jr(1:nTheta_north,:)

    ! map to the south hemisphere
    do i=1, nTheta_north
       energy_flux(nTheta_north+i-1,:)      = energy_flux(nTheta_north-i+1,:)
       ave_e(nTheta_north+i-1,:)            = ave_e(nTheta_north-i+1,:)
       num_flux(nTheta_north+i-1,:)         = num_flux(nTheta_north-i+1,:)
       dis_energy_flux(nTheta_north+i-1,:)  = dis_energy_flux(nTheta_north-i+1,:)
       dis_ave_e(nTheta_north+i-1,:)        = dis_ave_e(nTheta_north-i+1,:)
       !diff_flux_iono(nTheta_north+i-1,:,:) = diff_flux_iono(nTheta_north-i+1,:,:)
       Jr(nTheta_north+i-1,:)               = Jr(nTheta_north-i+1,:)
    end do

    ! in order to couple to the ionosphere module that starts the index at
    ! noon with phi=0, shift the array about 180 degree
    Energy_fluxtmp(:, 1:nPhi/2+1) = Energy_flux(:, 1:nPhi/2+1)
    Energy_flux(:,1:nPhi/2+1)     = Energy_flux(:,nPhi/2+1:nPhi)
    Energy_flux(:,nPhi/2+1:nPhi)  = Energy_fluxtmp(:,1:nPhi/2+1)
    Ave_etmp(:, 1:nPhi/2+1)       = Ave_e(:, 1:nPhi/2+1)
    Ave_e(:,1:nPhi/2+1)           = Ave_e(:,nPhi/2+1:nPhi)
    Ave_e(:,nPhi/2+1:nPhi)        = Ave_etmp(:,1:nPhi/2+1)
    Num_fluxtmp(:, 1:nPhi/2+1)    = Num_flux(:, 1:nPhi/2+1)
    Num_flux(:,1:nPhi/2+1)        = Num_flux(:,nPhi/2+1:nPhi)
    Num_flux(:,nPhi/2+1:nPhi)     = Num_fluxtmp(:,1:nPhi/2+1)

    dis_Energy_fluxtmp(:, 1:nPhi/2+1) = dis_Energy_flux(:, 1:nPhi/2+1)
    dis_Energy_flux(:,1:nPhi/2+1)     = dis_Energy_flux(:,nPhi/2+1:nPhi)
    dis_Energy_flux(:,nPhi/2+1:nPhi)  = dis_Energy_fluxtmp(:,1:nPhi/2+1)
    dis_Ave_etmp(:, 1:nPhi/2+1)       = dis_Ave_e(:, 1:nPhi/2+1)
    dis_Ave_e(:,1:nPhi/2+1)           = dis_Ave_e(:,nPhi/2+1:nPhi)
    dis_Ave_e(:,nPhi/2+1:nPhi)        = dis_Ave_etmp(:,1:nPhi/2+1)

    Jr_tmp(:,1:nPhi/2+1) = Jr(:,1:nPhi/2+1)
    Jr(:,1:nPhi/2+1)     = Jr(:,nPhi/2+1:nPhi)
    Jr(:,nPhi/2+1:nPhi)  = Jr_tmp(:,1:nPhi/2+1)

    high_latboundary_tmp(1:nPhi/2+1) = high_latboundary(1:nPhi/2+1)
    high_latboundary(1:nPhi/2+1)     = high_latboundary(nPhi/2+1:nPhi)
    high_latboundary(nPhi/2+1:nPhi)  = high_latboundary_tmp(1:nPhi/2+1)

    DEALLOCATE(ave_flux, f, ave_fluxext, num_fluxeq, NeExt, &
               PperEExt, PparEExt, rRawExt, aRawExt, colatGrid, lonGrid, F0, &
               energy_flux_iono, num_flux_iono, ave_e_iono, NeEq, PParEEq, PPerEEq, &
               dis_num_flux_iono, dis_ave_e_iono, dis_energy_flux_iono, ave_fluxEq, &
               colat, lon, xTemp, yTemp)
    deallocate(Energy_fluxtmp, Ave_etmp, Num_fluxtmp, dis_energy_fluxtmp, dis_ave_etmp, &
               Jr_tmp, high_latboundary_tmp, temp, rGrid, aGrid)
    DEALLOCATE(idX)
    return

  end subroutine calculate_precip_flux_jr

!==================================================================================================
END MODULE ModRamSce
