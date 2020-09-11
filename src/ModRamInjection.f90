module ModRamInjection

  use ModRamMain,      ONLY: DP

  implicit none

  contains

!==========================================================================
  subroutine load_injection_file(S, fname, fluxOut)

     use ModRamMain,      ONLY: PathRamIn
     use ModRamConst,     ONLY: Re, Mp, Q
     use ModRamGrids,     ONLY: nR, nT, nE, nPa, nS, dR, dPhi
     use ModRamVariables, ONLY: Lz, Phi, EkeV, Mu, dE, dMu, PA, RFactor, &
                                FFactor, FNHS, species
     use ModScbVariables, ONLY: bnormal
 
     use ModIoUnit,       ONLY: UNITTMP_

     use nrtype, ONLY: pi_d, twopi_d

     implicit none

     character(len=*), intent(in) :: fname
     integer, intent(in) :: S
     real(DP), dimension(:,:,:,:), intent(inout) :: fluxOut

     integer :: i, j, k, l, n, iR, iT, iE, iPa, iError
     character(len=200) filename, HEADER
     real(DP) :: x_i, y_i, z_i, d_i, e_i, a_i, factor, den
     real(DP) :: theta, p0, radius, r0, b_r, b_t, b_p, bmag, bmin
     real(DP) :: vpar, vper, alpha, alpha_0, energy, den1
     real(DP) :: XGEO, YGEO, ZGEO, XGSW, YGSW, ZGSW
     real(DP), allocatable :: num_density(:,:,:,:)

     allocate(num_density(nR,nT,nE,nPa))

     num_density = 0.0
     den1 = -1.0

     filename = trim(PathRamIn)//fname
     OPEN(UNITTMP_,FILE=trim(filename),STATUS='OLD')
     READ(UNITTMP_,'(A)') HEADER
     read_file: do
        ! Current file format assumes x, y, and z are in Re
        ! energy is in eV, and pitch angle is in degrees
        read(UNITTMP_,*,IOSTAT=IError) xGEO, yGEO, zGEO, e_i, a_i, d_i
        if (iError.ne.0) exit read_file

        ! Convert from GEO coordinates to SM coordinates
        call GEOGSW_08(XGEO,YGEO,ZGEO,xGSW,yGSW,zGSW,1)
        CALL SMGSW_08(x_i,y_i,z_i,xGSW,yGSW,zGSW,-1)

        r0 = sqrt(x_i**2 + y_i**2 + z_i**2) ! Re
        p0 = atan2(y_i,x_i)                 ! rad
        if (p0 < 0._dp) p0 = p0 + twopi_d
        energy  = e_i/1000.0                ! keV
        alpha_0 = a_i*pi_d/180.0            ! rad
        den     = d_i

        ! Find the correct R, T, E, Pa bin
        iR = 1
        do i = 2,nR
          if ((r0 > Lz(i-1)).and.(r0 < Lz(i))) then
             iR = i
             exit
          endif
        enddo

        iT = 1
        do i = 2,nT
          if ((p0 > Phi(i-1)).and.(p0 < Phi(i))) then
             iT = i
             exit
          endif
        enddo

        iE = 1
        do i = 2,nE
          if ((energy > EkeV(i-1)).and.(energy < EkeV(i))) then
             iE = i
             exit
          endif
        enddo

        iPa = 1
        do i = 2,nPa
          if ((cos(alpha_0) > Mu(i-1)).and.(cos(alpha_0) < Mu(i))) then
             iPa = i
             exit
          endif
        enddo

        ! Add the number density to the grid assuming a maxwellian distribution in energy 
        factor = 4.0E6*EkeV(iE)
        den = 1._dp/sqrt(species(S)%s_mass)*exp(-1._dp*EkeV(iE)/energy)*factor*den*energy**(-1.5)
        !den = FFactor(S,iR,iE,iPa)*FNHS(iR,iT,iPa)*den/dE(iE)/dMu(iPa)
        num_density(iR,iT,iE,iPa) = num_density(iR,iT,iE,iPa) + den*FFactor(S,iR,iE,iPa)*FNHS(iR,iT,iPa)
     enddo read_file
     CLOSE(UNITTMP_)

     fluxOut(:,:,:,:) = num_density

     deallocate(num_density)
     return

  end subroutine load_injection_file

end module ModRamInjection
