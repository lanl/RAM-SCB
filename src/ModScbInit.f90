MODULE ModScbInit
  ! Contains subroutines for initializattion of SCB component
  
  use ModScbVariables, ONLY: f, fp, rhoVal, chiVal, thetaVal, zetaVal, fzet, &
                             fzetp, xzero, xzero3, dphi, bnormal, pnormal, &
                             enormal, pjconst, psiin, psiout, psitot, bzero, &
                             constZ, constTheta, pressurequot, re1, xpsiin, &
                             xpsiout, r0Start, byimfglobal, bzimfglobal, &
                             pdynglobal, blendGlobal, blendGlobalInitial, &
                             p_dyn, by_imf, bz_imf, dst_global, wTsyg, tilt, &
                             nZetaMidnight, nThetaEquator, nMaximum
  
  implicit none
  
  contains
!==============================================================================
  SUBROUTINE scb_init
  
    use ModScbParams,    ONLY: blendAlphaInit, blendPsiInit
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap
    use ModScbVariables, ONLY: blendAlpha, blendPsi, alphaVal, alphaValInitial, &
                               psiVal, xpsiout, xpsiin
  
    USE ModScbIO,     ONLY: computational_domain
    use nrtype, ONLY: DP, pi_d, twopi_d

    IMPLICIT NONE
  
    INTEGER :: k, j, i, nc
  
    REAL(DP) :: xpsitot, psis, xpl, phi, aa, dphi
  
    REAL(DP), PARAMETER :: pow = 1.0_dp, TINY = 1.E-15_dp

    call computational_domain
  
    blendAlpha = blendAlphaInit
    blendPsi = blendPsiInit
  
    bzero = 1.0_dp
  
    NMAXimum = MAX(nthe, npsi, nzeta)
  
    ! Pressure
    pressurequot = 1._dp  ! 1 is quiet time; used for some quiet-time computations
    xzero = 6.6_dp
  
    ! Additional parameters
    nZetaMidnight = (nzeta+3) / 2
    nThetaEquator = nthe/2 + 1
  
    re1 = 1.0_dp
    xzero3 = xzero**3
  
    dphi  = twopi_d/REAL(nzeta-1, DP)
  
    ! Define (theta, rho, zeta) computational coordinates
    DO i = 1, nthe
       thetaVal(i) = pi_d * REAL(i-1, DP)/REAL(nthe-1, DP)
    END DO
  
    DO j = 1, npsi
       rhoVal(j) = REAL(j-1, DP)/REAL(npsi-1, DP)
    END DO
  
    DO k = 1, nzeta
       phi = REAL(k-2, DP) * dphi
       zetaVal(k) =  phi
    END DO
  
    chiVal = (thetaVal + constTheta * SIN(2._dp*thetaVal))
  
    DO k = 1, nzetap
       phi         = REAL(k-2, DP) * dphi
       alphaVal(k) = phi + constz*SIN(phi) ! Concentration at midnight
       fzet(k)     = 1._dp + constz*COS(phi)
       fzetp(k)    = - constz * SIN(phi)
  
       !Concentration at dusk-midnight (3pi/4), good for storm computations            
       !alphaVal(k)         = phi + constz*(SIN(phi+0.18_dp*pi_d) -
       !SIN(0.18_dp*pi_d))
       !fzet(k)           = 1._dp + constz*COS(phi+0.18_dp*pi_d)
       !fzetp(k)          = - constz*SIN(phi+0.18_dp*pi_d)
    END DO
    alphaValInitial = alphaVal
  
    !cc...dipole field line in terms of polar coordinate given by
    !r=x*cos(theta)**2
    !c..  B = grad(psi) X grad(alpha), psi = -x0**3 * cos(theta)**2 / r
  
    !c... x0 is the distance from earth in the equatorial plane
    !cc..  B = (x0/r)**3 * [ cos(theta) the-dir - 2.* sin(theta) r-dir]
    !cc..  normalize B(r=xzero)=B0=1.0 at equator with x0=xzero
    !cc..  define psival on j grids
  
    !cc.. bnormal is Earth's dipole magnetic field at xzero in equator (in nT)
    !cc.. enormal is the normalization unit for the electric field (to give
    !E-field in mV/m)
    !cc.. for xzero = 6.6 R_E, bnormal = 107.83 nT
    !cc.. pnormal = bnormal**2 in nPa; for xzero = 6.6 R_E, pnormal = 9.255 nPa
    bnormal = 0.31_dp / xzero3 * 1.E5_dp
    enormal = bnormal * 6.4
    pnormal = 7.958E-4_dp * bnormal*bnormal  ! 0.01/(4pi) in front of bnormal**2
  
    !cc.. p0 is in unit of pnormal
    !cc  For Earth's surface dipole field B_D=0.31e-4 T, R_E=6.4e6 m, permeability
    != 4.*pi*1.e-7 H/m
    !cc  The unit conversion constant pjconst = 0.0134
    !cc  The current is in unit of (microA/m**2) by multiplying with pjconst
    pjconst = 1.e6_dp * 0.31E-4_dp / (xzero3 * 4._dp * pi_d * 1.E-7_dp * 6.4E6_dp)
  
    psiin  = -xzero3/xpsiin
    psiout = -xzero3/xpsiout
    psitot = psiout-psiin
  
    aa = -xzero3
  
    xpsitot = xpsiout - xpsiin
  
    !c  Need to make sure that psival is a monotonically increasing function of j
    !c  define psival grids that correspond to dipole psivals for j=1 and j=npsi
    !c  define psival grids that correspond to equal equatorial distance grids in
    !the midnight sector
  
    !c...  the (i, j, k) coordinate system is related to (theta,psis,zeta) flux
    !coordinate
    !c...   by theta = (i-1) * pi / (nthe-1), 0 < theta < pi
    !c...   by zeta = phi + delta(i,j,k), 0 < zeta < 2*pi
    !c...   by psis=(j-1) / (npsi-1), 0 < psis < 1.0
  
    !c...      psiin < psival < psiout
  
    !c...  f(j) = d(psival)/d(psis)
    !c..   fp = d(f)/d(psis)
  
    DO j = 1, npsi
       psis = REAL(j-1, DP) / REAL(npsi-1, DP)
       xpl = xpsiin + xpsitot * psis**pow
       psival(j) = aa / xpl
       f(j) = (-aa / xpl**2) * xpsitot * pow * psis**(pow-1.)
       fp(j) = 0._dp ! If pow = 1
    END DO
  
    RETURN
  
  END SUBROUTINE
  
  END MODULE ModScbInit
