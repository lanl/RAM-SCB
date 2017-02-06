SUBROUTINE input

  ! Initializes various quantities

  USE nrtype
  USE Module1

  IMPLICIT NONE

  INTEGER :: k, j, i, nc
  REAL(DP) :: xpsitot, psis, xpl, phi, aa
  REAL(DP), PARAMETER :: pow = 1.0_dp
  REAL(DP), PARAMETER :: TINY = 1.E-15_dp

  bzero = 1.0_dp

  NMAXimum = MAX(nthe, npsi, nzeta) 

  ! Pressure
  pressurequot = 1._dp  ! 1 is quiet time; used for some quiet-time computations
  xzero = 6.6_dp

  ! Additional parameters
  nZetaMidnight = (nzeta+3) / 2  
  nThetaEquator = nthe/2 + 1

  re               = 1.0_dp
  xzero3           = xzero**3

  nthem = nthe - 1
  nthep = nthe + 1
  npsim = npsi - 1
  npsip = npsi + 1
  nzetm = nzeta - 1
  nzetp = nzeta + 1
  dphi = twopi_d / REAL(nzeta - 1, DP)

  ! Define (theta, rho, zeta) computational coordinates
  DO i = 1, nthe
     thetaVal(i) = pi_d * REAL(i - 1, DP) / REAL(nthe - 1, DP)
  END DO

  DO j = 1, npsi
     rhoVal(j) = REAL(j - 1, DP) / REAL(npsi - 1, DP)
  END DO

  DO k                 = 1, nzeta
     phi               = REAL(k-2, DP) * dphi
     zetaVal(k)        =  phi   
  END DO

  chiVal = (thetaVal + constTheta * SIN(2._dp*thetaVal)) 

  IF (iUseSavedData /= 1) THEN
 !C    IF (rank == 0) OPEN(10, file = 'zeta_values', action = 'write', status = 'replace')
     DO k  = 1, nzetp
        phi               = REAL(k-2, DP) * dphi
        alphaVal(k)       = phi + constz*SIN(phi) ! Concentration at midnight
        fzet(k)           = 1._dp + constz*COS(phi)
        fzetp(k)          = - constz * SIN(phi)

        !Concentration at dusk-midnight (3pi/4), good for storm computations            
        !alphaVal(k)         = phi + constz*(SIN(phi+0.18_dp*pi_d) - SIN(0.18_dp*pi_d))
        !fzet(k)           = 1._dp + constz*COS(phi+0.18_dp*pi_d)
        !fzetp(k)          = - constz*SIN(phi+0.18_dp*pi_d)

!C        IF (rank == 0) WRITE(10, '(I3, 1X, 3G14.7)') k, alphaVal(k), fzet(k), fzetp(k)
     END DO
!C     IF (rank == 0) CLOSE(10)
  END IF

  alphaValInitial = alphaVal

  !cc...dipole field line in terms of polar coordinate given by r=x*cos(theta)**2
  !c..  B = grad(psi) X grad(alpha), psi = -x0**3 * cos(theta)**2 / r

  !c... x0 is the distance from earth in the equatorial plane
  !cc..  B = (x0/r)**3 * [ cos(theta) the-dir - 2.* sin(theta) r-dir]
  !cc..  normalize B(r=xzero)=B0=1.0 at equator with x0=xzero
  !cc..  define psival on j grids

  !cc.. bnormal is Earth's dipole magnetic field at xzero in equator (in nT)
  !cc.. enormal is the normalization unit for the electric field (to give E-field in mV/m)
  !cc.. for xzero = 6.6 R_E, bnormal = 107.83 nT
  !cc.. pnormal = bnormal**2 in nPa; for xzero = 6.6 R_E, pnormal = 9.255 nPa
  bnormal        = 0.31_dp / xzero3 * 1.E5_dp
  enormal        = bnormal * 6.4
  pnormal        = 7.958E-4_dp * bnormal*bnormal  ! 0.01/(4pi) in front of bnormal**2

  !cc.. p0 is in unit of pnormal
  !cc  For Earth's surface dipole field B_D=0.31e-4 T, R_E=6.4e6 m, permeability = 4.*pi*1.e-7 H/m
  !cc  The unit conversion constant pjconst = 0.0134
  !cc  The current is in unit of (microA/m**2) by multiplying with pjconst
  pjconst      = 1.e6_dp * 0.31E-4_dp / (xzero3 * 4._dp * pi_d * 1.E-7_dp * 6.4E6_dp)

  psiin        = -xzero3/xpsiin
  psiout       = -xzero3/xpsiout
  psitot       = psiout-psiin

  aa           = -xzero3

  xpsitot = xpsiout - xpsiin

  !c  Need to make sure that psival is a monotonically increasing function of j
  !c  define psival grids that correspond to dipole psivals for j=1 and j=npsi
  !c  define psival grids that correspond to equal equatorial distance grids in the midnight sector

  !c...  the (i, j, k) coordinate system is related to (theta,psis,zeta) flux coordinate
  !c...   by theta = (i-1) * pi / (nthe-1), 0 < theta < pi
  !c...   by zeta = phi + delta(i,j,k), 0 < zeta < 2*pi
  !c...   by psis=(j-1) / (npsi-1), 0 < psis < 1.0

  !c...      psiin < psival < psiout

  !c...  f(j) = d(psival)/d(psis)
  !c..   fp = d(f)/d(psis)

  !C  IF (iUseSavedData /= 1) THEN  ! If data not from saved previous configuration
  !C IF (rank == 0) OPEN(10, file = 'psi_values', action = 'write', status = 'replace')

  DO j = 1, npsi
     psis = REAL(j-1, DP) / REAL(npsi-1, DP)
     xpl = xpsiin + xpsitot * psis**pow
     !C        xpl = xpsiin + xpsitot * 0.5_dp * (psis + sin(0.5_dp*pi_d*psis))
     ! choice of psival(j) grids ! It has to be the same one defined in the tracing program !!!
     psival(j) = aa / xpl 

     !C   if (j > 1) then 
     f(j) = (-aa / xpl**2) * xpsitot * pow * psis**(pow-1.)
     !C          fp(j) = (2. * aa / xpl**3) * (pow*xpsitot*psis**(pow-1))**2 - aa/xpl**2 * pow*(pow-1.)*xpsitot*psis**(pow-2.) 
     fp(j) = 0._dp ! If pow = 1
     !C  end if

     !        f(j) = (-aa/xpl**2) * xpsitot * 0.5_dp * (1._dp + 0.5_dp*pi_d * cos(0.5_dp*pi_d*psis))
     !        fp(j)  = 0._dp ! Not needed for this calculation

  !C   IF (rank == 0) THEN
  !C      WRITE(10, '(I3, 1X, 3G14.7)') j, psival(j), f(j), fp(j); CALL FLUSH(10)
        !C WRITE(*, '(I3, 1X, 3G14.7)') j, psival(j), f(j), fp(j); CALL FLUSH(6)
  !C   END IF
  END DO
  !C     f(1) = f(2)
  !C     fp(1) = fp(2)

 !C IF (rank == 0) CLOSE(10)
  !C  END IF

  ! Various quantities for the finite difference schemes
  dr =               1._dp/REAL(npsi - 1, dp)      ! d(rho)
  dt =               pi_d/REAL(nthe - 1, dp)        ! d(theta)
  dpPrime = twopi_d / REAL(nzeta -1, dp)   ! d(zeta)  (because dp has been defined as 8, double prec.)
  rdr = 1._dp / dr
  rdt = 1._dp / dt
  rdp = 1._dp / dpPrime
  rdrsq = rdr**2
  rdtsq = rdt**2
  rdpsq = rdp**2
  rdr2 = 0.5_dp * rdr
  rdt2 = 0.5_dp * rdt
  rdp2 = 0.5_dp * rdp
  rdr4 = 0.25_dp * rdr
  rdt4 = 0.25_dp * rdt
  rdp4 = 0.25_dp * rdp
  rdrdp4 = 0.25_dp * rdr * rdp
  rdpdt4 = 0.25_dp * rdp * rdt
  rdtdr4 = 0.25_dp * rdt * rdr


  RETURN 

END SUBROUTINE input
