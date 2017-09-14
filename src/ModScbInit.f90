MODULE ModScbInit
  ! Contains subroutines for initializattion of SCB component
  
  use ModScbVariables!, ONLY: f, fp, rhoVal, chiVal, thetaVal, zetaVal, fzet, &
                     !        fzetp, xzero, xzero3, dphi, bnormal, pnormal, &
                     !        enormal, pjconst, psiin, psiout, psitot, bzero, &
                     !        constZ, constTheta, pressurequot, re1, xpsiin, &
                     !        xpsiout, r0Start, byimfglobal, bzimfglobal, &
                     !        pdynglobal, blendGlobal, blendGlobalInitial, &
                     !        p_dyn, by_imf, bz_imf, dst_global, wTsyg, tilt, &
                     !        nZetaMidnight, nThetaEquator, nMaximum
  
  implicit none
  
  contains
!==============================================================================
  subroutine scb_allocate

    use ModRamGrids, ONLY: nR, nT, nPa
    use ModScbGrids

    implicit none

!!!!! Allocate Arrays
! Main SCB Variables
    ALLOCATE(jacobian(nthe,npsi,nzeta), bX(nthe,npsi,nzeta), bY(nthe,npsi,nzeta), &
             bZ(nthe,npsi,nzeta), bXIntern(nthe,npsi,nzeta), bYIntern(nthe,npsi,nzeta), &
             bZIntern(nthe,npsi,nzeta), x(nthe,npsi,nzeta+1), y(nthe,npsi,nzeta+1), &
             z(nthe,npsi,nzeta+1), xx(nthe,npsi,nzeta+1), yy(nthe,npsi,nzeta+1), &
             bf(nthe,npsi,nzeta+1), bsq(nthe,npsi,nzeta+1), fluxVolume(npsi,nzeta))
! ModScbEuler Variables
    ALLOCATE(psi(nthe,npsi,nzeta+1), psiSav1(nthe,npsi,nzeta+1), psiSav2(nthe,npsi,nzeta+1), &
             alfa(nthe,npsi,nzeta+1), alfaSav1(nthe,npsi,nzeta+1), alfaSav2(nthe,npsi,nzeta+1), &
             alphaVal(nzeta+1), alphaValInitial(nzeta+1), psiVal(npsi))
! ModScbInit Variables
    ALLOCATE(f(npsi), fp(npsi), rhoVal(npsi), chiVal(nthe), thetaVal(nthe), zetaVal(nzeta), &
             fzet(nzeta+1), fzetp(nzeta+1))
! ModScbRun Variables
    ALLOCATE(radEqMidNew(npsi), gradRhoSq(nthe,npsi,nzeta), gradZetaSq(nthe,npsi,nzeta), &
             gradRhoGradZeta(nthe,npsi,nzeta), gradRhoGradTheta(nthe,npsi,nzeta), &
             gradThetaGradZeta(nthe,npsi,nzeta), gradPsiGradAlpha(nthe,npsi,nzeta), &
             dPperdAlpha(nthe,npsi,nzeta), dBBdAlpha(nthe,npsi,nzeta), dBBdPsi(nthe,npsi,nzeta), &
             dPPerdPsi(nthe,npsi,nzeta), dPPardAlpha(nthe,npsi,nzeta), dPPardPsi(nthe,npsi,nzeta), &
             dPPerdTheta(nthe,npsi,nzeta), dPPerdRho(nthe,npsi,nzeta), dPPerdZeta(nthe,npsi,nzeta), &
             dBsqdAlpha(nthe,npsi,nzeta), dBsqdPsi(nthe,npsi,nzeta), dBsqdTheta(nthe,npsi,nzeta), &
             gradThetaSq(nthe,npsi,nzeta), bfInitial(nthe,npsi,nzeta+1), pressure3D(nthe,npsi,nzeta+1), &
             bj(nthe,npsi,nzeta+1), phij(nthe,npsi,nzeta+1), ppar(nthe,npsi,nzeta+1), &
             pper(nthe,npsi,nzeta+1), tau(nthe,npsi,nzeta+1), sigma(nthe,npsi,nzeta+1), &
             dPdAlpha(nthe,npsi,nzeta+1), dPdPsi(nthe,npsi,nzeta+1), dSqPdAlphaSq(nthe,npsi,nzeta+1), &
             dSqPdPsiSq(nthe,npsi,nzeta+1))
! Extra Variables
    ALLOCATE(EXInd(nthe,npsi,nzeta+1), EYInd(nthe,npsi,nzeta+1), EZInd(nthe,npsi,nzeta+1), &
             EXConv(nthe,npsi,nzeta+1), EYConv(nthe,npsi,nzeta+1), EZConv(nthe,npsi,nzeta+1), &
             PhiIono(npsi,nzeta+1), dPhiIonodAlpha(npsi,nzeta+1), dPhiIonodBeta(npsi,nzeta+1), &
             radGrid(npsi,nzeta), angleGrid(npsi,nzeta), ratioEq(npsi,nzeta), dela(nzeta+1), &
             factor(nzeta+1), curvaturePsi(nthe,npsi,nzeta), curvatureGradPsi(nthe,npsi,nzeta), &
             gradAlphaSq(nthe,npsi,nzeta), radialCurvature(nthe,npsi,nzeta), dAlphadT(nthe,npsi,nzeta), &
             dBetadT(nthe,npsi,nzeta), bpar(nthe,npsi,nzeta+1), bper(nthe,npsi,nzeta+1), &
             jParDirect(nthe,npsi,nzeta+1))
! Mixed Variables
    ALLOCATE(hdens_Cart(NR,NT,NPA), h_Cart(NR,NT,NPA), h_Cart_interp(NR,NT,NPA), I_Cart(NR,NT,NPA), &
             I_Cart_interp(NR,NT,NPA), bZEq_Cart(NR,NT), flux_vol_Cart(NR,NT), radRaw(0:NR), &
             azimRaw(NT))
!!!!!

  nthem = nthe-1; nthep = nthe+1
  npsim = npsi-1; npsip = npsi+1
  nzetam = nzeta-1; nzetap = nzeta+1
  nq = nthe + 2
  ny = npsi + 2
  na = nzeta + 2
  nm = 2*nq
  nyp2 = ny + 2

  dr      = 1._dp/REAL(npsi - 1, dp)
  dt      = pi_d/REAL(nthe - 1, dp)
  dpPrime = 2*pi_d / REAL(nzeta -1, dp)
  rdr     = 1._dp / dr
  rdt     = 1._dp / dt
  rdp     = 1._dp / dpPrime
  rdrsq   = rdr**2
  rdtsq   = rdt**2
  rdpsq   = rdp**2
  rdr2    = 0.5_dp * rdr
  rdt2    = 0.5_dp * rdt
  rdp2    = 0.5_dp * rdp
  rdr4    = 0.25_dp * rdr
  rdt4    = 0.25_dp * rdt
  rdp4    = 0.25_dp * rdp
  rdrdp4  = 0.25_dp * rdr * rdp
  rdpdt4  = 0.25_dp * rdp * rdt
  rdtdr4  = 0.25_dp * rdt * rdr

  nAzimRAM = NT
  nXRaw    = NR-1
  nXRawExt = NR+1
  nYRaw    = NT

  end subroutine scb_allocate

!==============================================================================
  subroutine scb_deallocate

    implicit none

!!!!! Deallocate Arrays
! Main SCB Variables
    DEALLOCATE(jacobian, bX, bY, bZ, bXIntern, bYIntern, bZIntern, x, y, &
               z, xx, yy, bf, bsq, fluxVolume)
! ModScbEuler Variables
    DEALLOCATE(psi, psiSav1, psiSav2, alfa, alfaSav1, alfaSav2, alphaVal, &
               alphaValInitial, psiVal)
! ModScbInit Variables
    DEALLOCATE(f, fp, rhoVal, chiVal, thetaVal, zetaVal, fzet, fzetp)
! ModScbRun Variables
    DEALLOCATE(radEqMidNew, gradRhoSq, gradZetaSq, gradRhoGradZeta, gradRhoGradTheta, &
               gradThetaGradZeta, gradPsiGradAlpha, dPperdAlpha, dBBdAlpha, dBBdPsi, &
               dPPerdPsi, dPPardAlpha, dPPardPsi, dPPerdTheta, dPPerdRho, dPPerdZeta, &
               dBsqdAlpha, dBsqdPsi, dBsqdTheta, gradThetaSq, bfInitial, pressure3D, &
               bj, phij, ppar, pper, tau, sigma, dPdAlpha, dPdPsi, dSqPdAlphaSq, &
               dSqPdPsiSq)
! Extra Variables
    DEALLOCATE(EXInd, EYInd, EZInd, EXConv, EYConv, EZConv, PhiIono, dPhiIonodAlpha, &
               dPhiIonodBeta, radGrid, angleGrid, ratioEq, dela, factor, curvaturePsi, &
               curvatureGradPsi, gradAlphaSq, radialCurvature, dAlphadT, dBetadT, bpar, &
               bper, jParDirect)
! Mixed Variables
    DEALLOCATE(hdens_Cart, h_Cart, h_Cart_interp, I_Cart, I_Cart_interp, bZEq_Cart, &
               flux_vol_Cart, radRaw, azimRaw)
!!!!!

  end subroutine scb_deallocate

!==============================================================================
  SUBROUTINE scb_init
  
    use ModScbParams,    ONLY: blendAlphaInit, blendPsiInit
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables, ONLY: blendAlpha, blendPsi, alphaVal, alphaValInitial, &
                               psiVal, xpsiout, xpsiin

    use ModScbEuler, ONLY: psiges, alfges
    use ModScbIO,    ONLY: computational_domain

    use nrtype, ONLY: DP, pi_d, twopi_d

    IMPLICIT NONE
  
    INTEGER :: k, j, i, nc
  
    REAL(DP) :: xpsitot, psis, xpl, phi, aa, dphi
  
    REAL(DP), PARAMETER :: pow = 1.0_dp, TINY = 1.E-15_dp

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
  
    DO k = 1, nzeta+1
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
 
    call computational_domain

    !c  Need to make sure that psival is a monotonically increasing function of
    !j
    !c  define psival grids that correspond to dipole psivals for j=1 and j=npsi
    !c  define psival grids that correspond to equal equatorial distance grids
    !in
    !the midnight sector
    psiin   = -xzero3/xpsiin
    psiout  = -xzero3/xpsiout
    psitot  = psiout-psiin
    xpsitot = xpsiout - xpsiin
    DO j = 1, npsi
       psis = REAL(j-1, DP) / REAL(npsi-1, DP)
       xpl = xpsiin + xpsitot * psis**pow
       psival(j) = -xzero3 / xpl
       f(j) = (xzero3 / xpl**2) * xpsitot * pow * psis**(pow-1.)
       fp(j) = 0._dp ! If pow = 1
    END DO

    call psiges
    call alfges

    RETURN
  
  END SUBROUTINE
  
  END MODULE ModScbInit
