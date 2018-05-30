!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================
Module ModScbVariables
  use nrtype, ONLY: DP


  implicit none

  integer :: left, right, nFail = 0
  logical :: SORFail = .false., hICalc = .true.

! Figure out
  REAL(DP), ALLOCATABLE :: EXInd(:,:,:), EYInd(:,:,:), EZInd(:,:,:), EXConv(:,:,:), &
                           EYConv(:,:,:), EZConv(:,:,:), PhiIono(:,:), dPhiIonodAlpha(:,:), &
                           dPhiIonodBeta(:,:), radGrid(:,:), angleGrid(:,:), ratioEq(:,:), &
                           dela(:), factor(:), curvaturePsi(:,:,:), curvatureGradPsi(:,:,:), &
                           gradAlphaSq(:,:,:), radialCurvature(:,:,:), dAlphadT(:,:,:), &
                           dBetadT(:,:,:), bpar(:,:,:), bper(:,:,:), jParDirect(:,:,:)
!
  real(DP), ALLOCATABLE :: h_Cart(:,:,:), h_Cart_interp(:,:,:), I_Cart(:,:,:), &
                           I_Cart_interp(:,:,:), hdens_Cart(:,:,:), bZEq_Cart(:,:), &
                           flux_vol_Cart(:,:), radRaw(:), azimRaw(:), fluxVolume(:,:)
!
  REAL(DP) :: sumb, sumdb, normDiff, normDiffPrev, normJxB, normGradP
!

!
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: vecd, vec1, vec2, vec3, vec4, vec6, &
                                               vec7, vec8, vec9, vecx, vecr
!

!--- SCE Components
  real(DP), allocatable :: paraj(:,:)
!---

! MAIN SCB VARIABLES
  REAL(DP), ALLOCATABLE :: jacobian(:,:,:), bX(:,:,:), bY(:,:,:), bZ(:,:,:), bxintern(:,:,:), &
                           byintern(:,:,:), bzintern(:,:,:), x(:,:,:), y(:,:,:), z(:,:,:), &
                           xx(:,:,:), yy(:,:,:), bf(:,:,:), bsq(:,:,:)
!

! ModScbEuler Variables
  REAL(DP), ALLOCATABLE :: psi(:,:,:), psiSav1(:,:,:), psiSav2(:,:,:), alfa(:,:,:), &
                           alfaSav1(:,:,:), alfaSav2(:,:,:), alphaVal(:), alphaValInitial(:), &
                           psiVal(:), psiPrev(:,:,:), alfaPrev(:,:,:)
  REAL(DP) :: blendAlpha, blendPsi, diffmx, rjac, decreaseConvAlpha, &
              decreaseConvPsi, errorAlpha, errorAlphaPrev, errorPsi, &
              errorPsiPrev
  integer :: iAlphaMove, iPsiMove
!

! ModScbInit Variables
  REAL(DP), ALLOCATABLE :: f(:), fp(:), rhoVal(:), chiVal(:), thetaVal(:), zetaVal(:), &
                           fzet(:), fzetp(:), chi(:,:,:)
  REAL(DP) :: xzero, xzero3, dphi, bnormal, pnormal, enormal, pjconst, psiin, &
              psiout, psitot, bzero, constZ, constTheta, pressurequot, re1, &
              xpsiin, xpsiout, r0Start, byimfglobal, bzimfglobal, pdynglobal, &
              blendGlobal, blendGlobalInitial
  REAL(DP) :: p_dyn      = 0._dp, &
              by_imf     = 0._dp, &
              bz_imf     = 0._dp, &
              dst_global = 0._dp, &
              wTsyg(1:6) = 0._dp, &
              tilt       = 0._dp
  integer :: nZetaMidnight, nThetaEquator, nMaximum
!

! Convergence Variables
  REAL(DP), ALLOCATABLE :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:), GradPx(:,:,:), GradPy(:,:,:), &
                           GradPz(:,:,:), JxBx(:,:,:), JxBy(:,:,:), JxBz(:,:,:), &
                           JCrossB(:,:,:), GradP(:,:,:), jGradRho(:,:,:), jGradZeta(:,:,:), &
                           jGradTheta(:,:,:)
!
! Derivative Variables
  REAL(DP), ALLOCATABLE :: gradRhoSq(:,:,:), gradZetaSq(:,:,:), gradRhoGradZeta(:,:,:), &
                           gradRhoGradTheta(:,:,:), gradThetaGradZeta(:,:,:), gradThetaSq(:,:,:), &
                           gradRhoX(:,:,:), gradZetaX(:,:,:), gradThetaX(:,:,:), &
                           gradRhoY(:,:,:), gradZetaY(:,:,:), gradThetaY(:,:,:), &
                           gradRhoZ(:,:,:), gradZetaZ(:,:,:), gradThetaZ(:,:,:), &
                           derivXRho(:,:,:), derivXZeta(:,:,:), derivXTheta(:,:,:), &
                           derivYRho(:,:,:), derivYZeta(:,:,:), derivYTheta(:,:,:), &
                           derivZRho(:,:,:), derivZZeta(:,:,:), derivZTheta(:,:,:)
!

! ModScbRun Variables
  REAL(DP), ALLOCATABLE :: radEqMidNew(:), gradPsiGradAlpha(:,:,:), &
                           dPperdAlpha(:,:,:), dBBdAlpha(:,:,:), dBBdPsi(:,:,:), &
                           dPPerdPsi(:,:,:), dPPardAlpha(:,:,:), dPPardPsi(:,:,:), &
                           dPPerdTheta(:,:,:), dPPerdRho(:,:,:), dPPerdZeta(:,:,:), &
                           dBsqdAlpha(:,:,:), dBsqdPsi(:,:,:), dBsqdTheta(:,:,:), &
                           bfInitial(:,:,:), pressure3D(:,:,:), bj(:,:,:), phij(:,:,:), &
                           ppar(:,:,:), pper(:,:,:), tau(:,:,:), sigma(:,:,:), &
                           dPdAlpha(:,:,:), dPdPsi(:,:,:), dSqPdAlphaSq(:,:,:), &
                           dSqPdPsiSq(:,:,:), dBsqdRho(:,:,:), dBsqdZeta(:,:,:)
  REAL(DP) :: DstDPS, DstDPSInsideGeo, DstBiot, DstBiotInsideGeo
  INTEGER :: kmax, nisave, nitry, iteration, iConvGlobal, lconv
!

End Module ModScbVariables
