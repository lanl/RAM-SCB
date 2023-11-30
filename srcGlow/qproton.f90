! Subroutine QPROTON calculates ionization rates from energetic proton precipitation
! using the parameterization of Fang et al. (2010, 2013).

! Driver interface to QPROT converts units to mks and loops over altitude and energy

! Initial version of wrapper routine, Stan Solomon, 2/2018
! Removed dissociative ionization and re-ordered species ionization rates, Stan Solomon, 2/2018
! Converted input to array of proton flux as a function of energy, Stan Solomon, 2/2018

! Input:
!    nmaj  : number of major species (typically 3, O, O2, and N2)
!    jmax  : number of altitude levels
!    npbins: number of energy bins
!    pflux : array of number flux of incident protons in each energy bin [cm^-2 sec^-1 eV^-1]
!    pener : energy grid, bin centers (eV)
!    pdel  : energy bin widths (eV)
!    zz    : altitude array [cm]
!    ztn   : neutral temperature array [K]
!    zo    : atomic oxygen density array [cm^-3]
!    zo2   : molecular oxygen density array [cm^-3]
!    zn2   : molecular nitrogen density array [cm^-3]
!    zns   : atomic nitrogen density array [cm^-3]

! Output: 
!    pia   : ionization rate array for each major species (O, O2, N2) at each altitude [cm^-3/sec]

! Internal variables used by QPROT:
!    fmono : scalar energy flux of monoenergetic incident protons [eV/cm^2/sec]
!    emono : scalar energy of monoenergetic incident protons [eV]
!    h     : scale height [m]
!    rhomasstot  : total mass density [kg/m^3]
!    qtot  : ionization rate [1/cm^3/sec]


  subroutine qproton (nmaj, jmax, npbins, pflux, pener, pdel, zz, ztn, zo, zo2, zn2, zns, pia)

    implicit none

    integer, intent(in) :: nmaj, jmax, npbins
    real, intent(in)    :: pflux(npbins), pener(npbins), pdel(npbins)
    real, intent(in)    :: zz(jmax), ztn(jmax), zo(jmax), zo2(jmax), zn2(jmax), zns(jmax)
    real, intent(out)   :: pia(nmaj,jmax)

    integer :: j,n
    real :: fmono, emono
    real :: alt, rhoo, rhoo2, rhon2, rhons, rhototal, rhomasstot, meanmass, g, h, qtot, delta
    real, parameter :: avon = 6.0221413d23  ! Avogadro's #, atoms or molecules/mol
    real, parameter :: bolc = 1.3806488d-23 ! Boltzmann's constant, m^2 kg s^-2 K^-1
    real, parameter :: g0 = 9.80665         ! standard gravity, m^2/s
    real, parameter :: re = 6371009.        ! Mean Earth radius, m

    pia(:,:) = 0.

    ! Loop over altitudes:

    do j=1,jmax

      ! Convert units to mks:

      alt=zz(j)*1e-2
      rhoo=zo(j)*1e6
      rhoo2=zo2(j)*1e6
      rhon2=zn2(j)*1e6
      rhons=zns(j)*1e6
      rhototal = rhoo + rhoo2 + rhon2 + rhons                         !  #/m^3
      rhomasstot=(rhoo*.016+rhoo2*.032+rhon2*.028+rhons*.014) / avon  !  kg/m^3

      ! Calculate scale height:

      meanmass = rhomasstot/rhototal
      g = g0*(re**2 / (re+alt)**2)
      h = (bolc*ztn(j))/(meanmass*g)                                  !  m

      ! Loop over energies:

        do n=1,npbins

          fmono = pflux(n) * pener(n) * pdel(n)
          emono = pener(n)

          ! Call qprot to calculate total ionization rate qtot at each altitude:

          call qprot(fmono,emono,h,rhomasstot,qtot)

          ! Increment species-specific proton ioniztion rate using qtot:
          ! Note that we are neglecting direct ionization of N.

          delta = 0.92*rhoN2 + 1.50*rhoO2 + 0.56*rhoO
          pia(1,j) = pia(1,j) + qtot * (0.56*rhoO)  / delta
          pia(2,j) = pia(2,j) + qtot * (1.50*rhoO2) / delta
          pia(3,j) = pia(3,j) + qtot * (0.92*rhoN2) / delta

        enddo

    enddo

  end subroutine qproton

! -----------------------------------------------------

! Subroutine QPROT

! Calculates total ionization rate due to proton precipitation for a single altitude and energy.
! Range of validity is 100 eV to 1 MeV; ionization rate is set to zero if outside of range.

! Algorithm provided by Xiaohua Fang (see Fang et al., 2010, 2013).
! All intellectual property is that of Dr. Xiaohua Fang (CU/LASP).

! Untested version obtained from Ryan McGranaghan, 6/2016
! Initial adaptation for testing, Stan Solomon, 6/2016
! Refactored to free-form Fortran 90, Stan Solomon, 2/2018
! Moved inside wrapper routine, Stan Solomon, 2/2018

! Note mixed units: scale heigh and density input is mks
!                   proton flux input is cgs
!                   output ionization rate is cgs
 
! Input:
!    fmono : energy flux of incident protons [eV/cm^2/sec]
!    emono : monoenergy of incident protons  [eV]
!    h     : scale height [m]
!    rhomasstot  : total mass density [kg/m^3]

! Output: 
!    qtot  : ionization rate [1/cm^3/sec]
  
! References:
!    Fang, X., C. Randall, D. Lummerzheim, W. Wang, G. Lu, S. Solomon, and R. Frahm (2010),
!    Parameterization of monoenergetic electron impact ionization, Geophys. Res. Lett., 37, 
!    L22106, doi:10.1029/2010GLO45406.
 
!    Fang, X., D. Lummerzheim, and C. Jackman (2013), Proton impact ionization and a fast 
!    calculation method, J. Geophys. Res. Space Physics, 118, 5369-5378, doi:10.1002/jgra.50484

 
  SUBROUTINE QPROT (fmono, emono, h, rhomasstot, qtot)

    implicit none

    real, intent(in) :: fmono, emono 
    real, intent(in) :: h, rhomasstot
    real, intent(out) :: qtot
    real :: p(4,12), c(12), si, so, lne
 
    DATA  p / 2.55050e0,   2.69476e-1, -2.58425e-1,   4.43190e-2, &
             6.39287e-1,  -1.85817e-1, -3.15636e-2,   1.01370e-2, &
              1.63996e0,   2.43580e-1,  4.29873e-2,   3.77803e-2, &
            -2.13479e-1,   1.42464e-1,  1.55840e-2,   1.97407e-3, &
            -1.65764e-1,   3.39654e-1, -9.87971e-3,   4.02411e-3, &
            -3.59358e-2,   2.50330e-2, -3.29365e-2,   5.08057e-3, &
            -6.26528e-1,    1.46865e0,  2.51853e-1,  -4.57132e-2, &
              1.01384e0,   5.94301e-2, -3.27839e-2,   3.42688e-3, &
            -1.29454e-6,  -1.43623e-1,  2.82583e-1,   8.29809e-2, &
            -1.18622e-1,   1.79191e-1,  6.49171e-2,  -3.99715e-3, &
              2.94890e0,  -5.75821e-1,  2.48563e-2,   8.31078e-2, &
            -1.89515e-1,   3.53452e-2,  7.77964e-2,  -4.06034e-3/

    if ((emono > 99.) .and. (emono < 1.01e6)) then
      lne = log(emono*1.e-3)
      c=exp(p(1,:)+p(2,:)*lne+p(3,:)*(lne**2)+p(4,:)*(lne**3))
      si=log(7.5e3/emono*(rhomasstot*h*1000.)**9.e-1)
      so=c(1)*exp(c(2)*si-c(3)*exp(c(4)*si))   &
        +c(5)*exp(c(6)*si-c(7)*exp(c(8)*si))   &
        +c(9)*exp(c(10)*si-c(11)*exp(c(12)*si))
      qtot=so*fmono/h/3.5e3
    else
      qtot=0.
    endif

  END SUBROUTINE QPROT
