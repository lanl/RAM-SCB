Module ModScbVariables

  use nrtype, ONLY: DP
  use ModScbGrids

  implicit none
  save

! Figure out
REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: &
    EXInd,             &
    EYInd,             &
    EZInd,             &
    EXConv,            &
    EYConv,            &
    EZConv
real(DP), dimension(npsi,nzeta+1)      :: &
    PhiIono,           &
    dPhiIonodAlpha,    &
    dPhiIonodBeta
real(DP), dimension(NR,NT,NPA)         :: &
    h_Cart,            &
    h_Cart_interp,     &
    I_Cart,            &
    I_Cart_interp,     &
    hdens_Cart
real(DP), dimension(NR,NT)             :: &
    bZEq_Cart,         &
    flux_vol_Cart
real(DP), dimension(npsi,nzeta)        :: &
    radGrid,           &
    angleGrid,         &
    ratioEq
real(DP), dimension(nzeta+1)           :: &
    dela,              &
    factor
REAL(DP) :: radRaw(0:nR), azimRaw(nYRaw)
REAL(DP) :: sumb, sumdb
REAL(DP), DIMENSION(nthe,npsi,nzeta)   :: &
    curvaturePsi,      &
    curvatureGradPsi,  &
    gradAlphaSq,       &
    radialCurvature,   &
    dAlphadT,          &
    dBetadT
REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: &
    bpar,              &
    bper,              &
    jParDirect
REAL(DP) :: &
    normDiff,          &
    normDiffPrev,      &
    normJxB,           &
    normGradP
!

! MAIN SCB VARIABLES
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jacobian, bX, bY, bZ, bxintern, &
                                          byintern, bzintern
  REAL(DP), DIMENSION(nthe,npsi,nzetap) :: x, y, z, xx, yy, bf, bsq
!

! ModScbEuler Variables
  REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: psi, psiSav1, psiSav2, alfa, &
                                            alfaSav1, alfaSav2
  REAL(DP), DIMENSION(nzeta+1) :: alphaVal, alphaValInitial
  REAL(DP), DIMENSION(npsi) :: psiVal
  REAL(DP) :: blendAlpha, blendPsi, diffmx, rjac, decreaseConvAlpha, &
              decreaseConvPsi, errorAlpha, errorAlphaPrev, errorPsi, &
              errorPsiPrev
  integer :: iAlphaMove, iPsiMove
!

! ModScbInit Variables
  REAL(DP), DIMENSION(npsi) :: f, fp, rhoVal
  REAL(DP), DIMENSION(nthe) :: chiVal, thetaVal
  REAL(DP), DIMENSION(nzeta) :: zetaVal
  REAL(DP), DIMENSION(nzetap) :: fzet, fzetp
  REAL(DP) :: xzero, xzero3, dphi, bnormal, pnormal, enormal, pjconst, psiin, &
              psiout, psitot, bzero, constZ, constTheta, pressurequot, re1, &
              xpsiin, xpsiout, r0Start, byimfglobal, bzimfglobal, pdynglobal, &
              blendGlobal, blendGlobalInitial
  REAL(DP) :: p_dyn                = 0._dp, &
              by_imf               = 0._dp, &
              bz_imf               = 0._dp, &
              dst_global           = 0._dp, &
              wTsyg(1:6)           = 0._dp, &
              tilt                 = 0._dp
  integer :: nZetaMidnight, nThetaEquator, nMaximum
!

! ModScbRun Variables
  REAL(DP), DIMENSION(npsi) :: radEqMidNew
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: gradRhoSq, gradZetaSq, gradRhoGradZeta, &
            gradRhoGradTheta, gradThetaGradZeta, gradPsiGradAlpha, dPperdAlpha, &
            dBBdAlpha, dBBdPsi, dPPerdPsi, dPPardAlpha, dPPardPsi, dPPerdTheta, &
            dPPerdRho, dPPerdZeta, dBsqdAlpha, dBsqdPsi, dBsqdTheta, gradThetaSq
  REAL(DP), DIMENSION(nthe,npsi,nzetap) :: bfInitial, pressure3D, bj, phij, ppar, &
            pper, tau, sigma, dPdAlpha, dPdPsi, dSqPdAlphaSq, dSqPdPsiSq
  REAL(DP) :: DstDPS, DstDPSInsideGeo, DstBiot, DstBiotInsideGeo

  INTEGER :: kmax, nisave, nitry, iteration, iConvGlobal, lconv
!

End Module ModScbVariables
