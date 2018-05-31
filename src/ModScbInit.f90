!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbInit
  ! Contains subroutines for initializattion of SCB component
  
  implicit none
  
  contains
!==============================================================================
  subroutine scb_allocate

    use ModRamGrids, ONLY: nR, nT, nPa
    use ModScbGrids
    use ModScbVariables
    implicit none

!--- SCE Components
    ALLOCATE(paraj(npsi,nzeta+1))
    paraj = 0.0
!---

!!!!! Allocate Arrays
! Main SCB Variables
    ALLOCATE(jacobian(nthe,npsi,nzeta), bX(nthe,npsi,nzeta), bY(nthe,npsi,nzeta), &
             bZ(nthe,npsi,nzeta), bXIntern(nthe,npsi,nzeta), bYIntern(nthe,npsi,nzeta), &
             bZIntern(nthe,npsi,nzeta), x(nthe,npsi,nzeta+1), y(nthe,npsi,nzeta+1), &
             z(nthe,npsi,nzeta+1), xx(nthe,npsi,nzeta+1), yy(nthe,npsi,nzeta+1), &
             bf(nthe,npsi,nzeta+1), bsq(nthe,npsi,nzeta+1), fluxVolume(npsi,nzeta))
    jacobian = 0.0; bX = 0.0; bY = 0.0; bZ = 0.0; bXIntern = 0.0; bYIntern = 0.0
    bZIntern = 0.0; x = 0.0; y = 0.0; z = 0.0; xx = 0.0; yy = 0.0; bf = 0.0; bsq = 0.0
    fluxVolume = 0.0
!
    ALLOCATE(vecx(nthe,npsi,nzeta), vecr(nthe,npsi,nzeta), vecd(nthe,npsi,nzeta), &
             vec1(nthe,npsi,nzeta), vec2(nthe,npsi,nzeta), vec3(nthe,npsi,nzeta), &
             vec4(nthe,npsi,nzeta), vec6(nthe,npsi,nzeta), vec7(nthe,npsi,nzeta), &
             vec8(nthe,npsi,nzeta), vec9(nthe,npsi,nzeta))
    vecx = 0.0; vecr = 0.0; vecd = 0.0; vec1 = 0.0; vec2 = 0.0; vec3 = 0.0
    vec4 = 0.0; vec6 = 0.0; vec7 = 0.0; vec8 = 0.0; vec9 = 0.0;
! ModScbEuler Variables
    ALLOCATE(psi(nthe,npsi,nzeta+1), psiSav1(nthe,npsi,nzeta+1), psiSav2(nthe,npsi,nzeta+1), &
             alfa(nthe,npsi,nzeta+1), alfaSav1(nthe,npsi,nzeta+1), alfaSav2(nthe,npsi,nzeta+1), &
             alphaVal(nzeta+1), alphaValInitial(nzeta+1), psiVal(npsi), psiPrev(nthe,npsi,nzeta+1), &
             alfaPrev(nthe,npsi,nzeta+1))
    psi = 0.0; psiSav1 = 0.0; psiSav2 = 0.0; alfa = 0.0; alfaSav1 = 0.0;
    alfaSav2 = 0.0; alphaVal = 0.0; alphaValInitial = 0.0; psiVal = 0.0;
    psiPrev = 0.0; alfaPrev = 0.0
! ModScbInit Variables
    ALLOCATE(f(npsi), fp(npsi), rhoVal(npsi), chiVal(nthe), thetaVal(nthe), zetaVal(nzeta), &
             fzet(nzeta+1), fzetp(nzeta+1), chi(nthe,npsi,nzeta+1))
    f = 0.0; fp = 0.0; rhoVal = 0.0; chiVal = 0.0; thetaVal = 0.0; zetaVal = 0.0
    fzet = 0.0; fzetp = 0.0; chi = 0.0;
! Convergence Variables
    ALLOCATE(Jx(nthe,npsi,nzeta), Jy(nthe,npsi,nzeta), Jz(nthe,npsi,nzeta), &
             GradPx(nthe,npsi,nzeta), GradPy(nthe,npsi,nzeta), GradPz(nthe,npsi,nzeta), &
             JxBx(nthe,npsi,nzeta), JxBy(nthe,npsi,nzeta), JxBz(nthe,npsi,nzeta), &
             JCrossB(nthe,npsi,nzeta), GradP(nthe,npsi,nzeta), jGradRho(nthe,npsi,nzeta), &
             jGradZeta(nthe,npsi,nzeta), jGradTheta(nthe,npsi,nzeta))
    Jx = 0.0; Jy = 0.0; Jz = 0.0; GradPx = 0.0; GradPy = 0.0; GradPz = 0.0
    JxBx = 0.0; JxBy = 0.0; JxBz = 0.0; JCrossB = 0.0; GradP = 0.0
    jGradRho = 0.0; jGradZeta = 0.0; jGradTheta = 0.0
! Derivative Variables
    ALLOCATE(gradRhoSq(nthe,npsi,nzeta), gradZetaSq(nthe,npsi,nzeta), &
             gradRhoGradZeta(nthe,npsi,nzeta), gradRhoGradTheta(nthe,npsi,nzeta), &
             gradThetaGradZeta(nthe,npsi,nzeta), gradThetaSq(nthe,npsi,nzeta), &
             gradRhoX(nthe,npsi,nzeta), gradZetaX(nthe,npsi,nzeta), gradThetaX(nthe,npsi,nzeta), &
             gradRhoY(nthe,npsi,nzeta), gradZetaY(nthe,npsi,nzeta), gradThetaY(nthe,npsi,nzeta), &
             gradRhoZ(nthe,npsi,nzeta), gradZetaZ(nthe,npsi,nzeta), gradThetaZ(nthe,npsi,nzeta), &
             derivXRho(nthe,npsi,nzeta), derivXZeta(nthe,npsi,nzeta), derivXTheta(nthe,npsi,nzeta), &
             derivYRho(nthe,npsi,nzeta), derivYZeta(nthe,npsi,nzeta), derivYTheta(nthe,npsi,nzeta), &
             derivZRho(nthe,npsi,nzeta), derivZZeta(nthe,npsi,nzeta), derivZTheta(nthe,npsi,nzeta))
    gradRhoSq = 0.0; gradZetaSq = 0.0; gradRhoGradZeta = 0.0; gradRhoGradTheta = 0.0
    gradThetaGradZeta = 0.0; gradThetaSq = 0.0; gradRhoX = 0.0; gradZetaX = 0.0; gradThetaX = 0.0
    gradRhoY = 0.0; gradZetaY = 0.0; gradThetaY = 0.0; gradRhoZ = 0.0; gradZetaZ = 0.0; gradThetaZ = 0.0
    derivXRho = 0.0; derivXZeta = 0.0; derivXTheta = 0.0; derivYRho = 0.0; derivYZeta = 0.0
    derivYTheta = 0.0; derivZRho = 0.0; derivZZeta = 0.0; derivZTheta = 0.0
!
! ModScbRun Variables
    ALLOCATE(radEqMidNew(npsi), gradPsiGradAlpha(nthe,npsi,nzeta), &
             dPperdAlpha(nthe,npsi,nzeta), dBBdAlpha(nthe,npsi,nzeta), dBBdPsi(nthe,npsi,nzeta), &
             dPPerdPsi(nthe,npsi,nzeta), dPPardAlpha(nthe,npsi,nzeta), dPPardPsi(nthe,npsi,nzeta), &
             dPPerdTheta(nthe,npsi,nzeta), dPPerdRho(nthe,npsi,nzeta), dPPerdZeta(nthe,npsi,nzeta), &
             dBsqdAlpha(nthe,npsi,nzeta), dBsqdPsi(nthe,npsi,nzeta), dBsqdTheta(nthe,npsi,nzeta), &
             bfInitial(nthe,npsi,nzeta+1), pressure3D(nthe,npsi,nzeta+1), &
             bj(nthe,npsi,nzeta+1), phij(nthe,npsi,nzeta+1), ppar(nthe,npsi,nzeta+1), &
             pper(nthe,npsi,nzeta+1), tau(nthe,npsi,nzeta+1), sigma(nthe,npsi,nzeta+1), &
             dPdAlpha(nthe,npsi,nzeta+1), dPdPsi(nthe,npsi,nzeta+1), dSqPdAlphaSq(nthe,npsi,nzeta+1), &
             dSqPdPsiSq(nthe,npsi,nzeta+1),dBsqdRho(nthe,npsi,nzeta),dBsqdZeta(nthe,npsi,nzeta))
    radEqMidNew = 0.0; gradPsiGradAlpha = 0.0
    dPperdAlpha = 0.0; dBBdAlpha = 0.0; dBBdPsi = 0.0; dPPerdPsi = 0.0
    dPPardAlpha = 0.0; dPPardPsi = 0.0; dPPerdTheta = 0.0; dPPerdRho = 0.0
    dPPerdZeta = 0.0; dBsqdAlpha = 0.0; dBsqdPsi = 0.0; dBsqdTheta = 0.0
    gradThetaSq = 0.0; bfInitial = 0.0; pressure3D = 0.0; bj = 0.0; phij = 0.0
    ppar = 0.0; pper = 0.0; tau = 0.0; sigma = 0.0; dPdAlpha = 0.0; dPdPsi = 0.0
    dSqPdAlphaSq = 0.0; dSqPdPsiSq = 0.0; dBsqdRho = 0.0; dBsqdZeta = 0.0
! Extra Variables
    ALLOCATE(EXInd(nthe,npsi,nzeta+1), EYInd(nthe,npsi,nzeta+1), EZInd(nthe,npsi,nzeta+1), &
             EXConv(nthe,npsi,nzeta+1), EYConv(nthe,npsi,nzeta+1), EZConv(nthe,npsi,nzeta+1), &
             PhiIono(npsi,nzeta+1), dPhiIonodAlpha(npsi,nzeta+1), dPhiIonodBeta(npsi,nzeta+1), &
             radGrid(npsi,nzeta), angleGrid(npsi,nzeta), ratioEq(npsi,nzeta), dela(nzeta+1), &
             factor(nzeta+1), curvaturePsi(nthe,npsi,nzeta), curvatureGradPsi(nthe,npsi,nzeta), &
             gradAlphaSq(nthe,npsi,nzeta), radialCurvature(nthe,npsi,nzeta), dAlphadT(nthe,npsi,nzeta), &
             dBetadT(nthe,npsi,nzeta), bpar(nthe,npsi,nzeta+1), bper(nthe,npsi,nzeta+1), &
             jParDirect(nthe,npsi,nzeta+1))
    EXInd = 0.0; EYInd = 0.0; EZInd = 0.0; EXConv = 0.0; EYConv = 0.0; EZConv = 0.0
    PhiIono = 0.0; dPhiIonodAlpha = 0.0; dPhiIonodBeta = 0.0; radGrid = 0.0
    angleGrid = 0.0; ratioEq = 0.0; dela = 0.0; factor = 0.0; curvaturePsi = 0.0
    curvatureGradPsi = 0.0; gradAlphaSq = 0.0; radialCurvature = 0.0; dAlphadT = 0.0
    dBetadT = 0.0; bpar = 0.0; bper = 0.0; jParDirect = 0.0
! Mixed Variables
    ALLOCATE(hdens_Cart(NR,NT,NPA), h_Cart(NR,NT,NPA), h_Cart_interp(NR,NT,NPA), I_Cart(NR,NT,NPA), &
             I_Cart_interp(NR,NT,NPA), bZEq_Cart(NR,NT), flux_vol_Cart(NR,NT), radRaw(0:NR), &
             azimRaw(NT))
    hdens_Cart = 0.0; h_Cart = 0.0; h_Cart_interp = 0.0; I_Cart = 0.0
    I_Cart_interp = 0.0; bZEq_Cart = 0.0; flux_vol_Cart = 0.0; radRaw = 0.0
    azimRaw = 0.0
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
  nXRawExt = NR+floor(1.5/(5./nR))
  nYRaw    = NT

  end subroutine scb_allocate

!==============================================================================
  subroutine scb_deallocate

    use ModScbVariables

    implicit none

!--- SCE Components
    DEALLOCATE(paraj)
!---

!!!!! Deallocate Arrays
! Main SCB Variables
    DEALLOCATE(jacobian, bX, bY, bZ, bXIntern, bYIntern, bZIntern, x, y, &
               z, xx, yy, bf, bsq, fluxVolume)
    DEALLOCATE(vecx, vecr, vecd, vec1, vec2, vec3, vec4, vec6, vec7, vec8, vec9)
! ModScbEuler Variables
    DEALLOCATE(psi, psiSav1, psiSav2, alfa, alfaSav1, alfaSav2, alphaVal, &
               alphaValInitial, psiVal, psiPrev, alfaPrev)
! ModScbInit Variables
    DEALLOCATE(f, fp, rhoVal, chiVal, thetaVal, zetaVal, fzet, fzetp, chi)
! Convergence Variables
    DEALLOCATE(Jx, Jy, Jz, GradPx, GradPy, GradPz, JxBx, JxBy, JxBz, JCrossB, GradP, &
               jGradRho, jGradZeta, jGradTheta)
! Derivative Variables
    DEALLOCATE(gradRhoSq, gradZetaSq, gradRhoGradZeta, gradRhoGradTheta, &
               gradThetaGradZeta, gradThetaSq, gradRhoX, gradZetaX, gradThetaX, &
               gradRhoY, gradZetaY, gradThetaY, gradRhoZ, gradZetaZ, gradThetaZ, &
               derivXRho, derivXZeta, derivXTheta, derivYRho, derivYZeta, derivYTheta, &
               derivZRho, derivZZeta, derivZTheta)
! ModScbRun Variables
    DEALLOCATE(radEqMidNew, gradPsiGradAlpha, dPperdAlpha, dBBdAlpha, dBBdPsi, &
               dPPerdPsi, dPPardAlpha, dPPardPsi, dPPerdTheta, dPPerdRho, dPPerdZeta, &
               dBsqdAlpha, dBsqdPsi, dBsqdTheta, bfInitial, pressure3D, &
               bj, phij, ppar, pper, tau, sigma, dPdAlpha, dPdPsi, dSqPdAlphaSq, &
               dSqPdPsiSq, dBsqdRho, dBsqdZeta)
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
 
    use ModRamParams,    ONLY: IsRestart
    use ModRamTiming,    ONLY: Dts, DtsMin
    use ModScbParams,    ONLY: blendAlphaInit, blendPsiInit
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables, ONLY: blendAlpha, blendPsi, alphaVal, alphaValInitial, &
                               psiVal, xpsiout, xpsiin, r0Start, nZetaMidnight, &
                               nThetaEquator, bzero, pressurequot, xzero, xzero3, &
                               re1, bnormal, pnormal, enormal, pjconst, zetaVal, &
                               thetaVal, rhoVal, nMaximum

    use ModRamFunctions, ONLY: ram_sum_pressure
    use ModRamScb,       ONLY: computehI
    use ModScbRun,       ONLY: scb_run
    use ModScbEuler,     ONLY: psiges, alfges
    use ModScbIO,        ONLY: computational_domain

    use nrtype, ONLY: DP, pi_d, twopi_d

    implicit none
  
    INTEGER :: k, j, i, nc
  
    REAL(DP) :: xpsitot, psis, xpl, phi, aa, dphi
  
    r0Start = 1.0

    ! Additional parameters
    nZetaMidnight = (nzeta+1)/2 + 1
    nThetaEquator = (nthe+1)/2

    blendAlpha = blendAlphaInit
    blendPsi = blendPsiInit

!!!!! Set Normalization Constants
    bzero = 1.0_dp

    NMAXimum = MAX(nthe, npsi, nzeta)

    ! Pressure
    pressurequot = 1._dp  ! 1 is quiet time; used for some quiet-time computations
    xzero = 6.6_dp

    re1 = 1.0_dp
    xzero3 = xzero**3
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
    pnormal = bnormal*bnormal/(4._dp * pi_d * 1.E-7_dp)*1.E-9_dp  ! Pb = B^2/(2*permeability)

    !cc.. p0 is in unit of pnormal
    !cc  For Earth's surface dipole field B_D=0.31e-4 T, R_E=6.4e6 m,
    !permeability=4.*pi*1.e-7 H/m
    !cc  The unit conversion constant pjconst = 0.0134
    !cc  The current is in unit of (microA/m**2) by multiplying with pjconst
    pjconst = 1.e6_dp * 0.31E-4_dp / (xzero3 * 4._dp * pi_d * 1.E-7_dp * 6.4E6_dp)
!!!!!

!!!!! Set computational coordinates grid
    ! Define (theta, rho, zeta) computational coordinates
    DO i = 1, nthe
       thetaVal(i) = pi_d * REAL(i-1, DP)/REAL(nthe-1, DP)
    END DO

    DO j = 1, npsi
       rhoVal(j) = REAL(j-1, DP)/REAL(npsi-1, DP)
    END DO

    dphi  = twopi_d/REAL(nzeta-1, DP)
    DO k = 1, nzeta
       phi = REAL(k-2, DP) * dphi
       zetaVal(k) =  phi
    END DO
!!!!!

    RETURN
  
  END SUBROUTINE

END MODULE ModScbInit
