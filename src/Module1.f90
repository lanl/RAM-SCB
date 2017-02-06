
MODULE Module1

  USE nrtype, ONLY : pi_d, SP, DP ! Definition of single/double precision type
  USE ModRamMain, ONLY : NR, NT

  IMPLICIT NONE

  !  nthe is the number of theta grid points 
  !  npsi is the number of psi surfaces 
  !  nzeta is the number of zeta surfaces 

  ! Grid size - must match the computational boundary input file (*.cdf)
  ! INTEGER, PARAMETER :: nthe = 42, npsi=25, nzeta=31, &
  INTEGER, PARAMETER :: nthe=51, npsi=35, nzeta=61, &
  ! INTEGER, PARAMETER :: nthe=61, npsi = 41, nzeta = 71, &
       nq = nthe + 2, ny = npsi + 2, na = nzeta + 2, nm = 2 * nq, &
       nyp2 = ny + 2 

  ! Number of LTs in RAM for coupling in RAM-SCB
  INTEGER, PARAMETER :: nAzimRAM = NT 
  ! Number of radial positions in RAM for coupling in RAM-SCB
  INTEGER, PARAMETER :: nXRaw = NR-1 
  INTEGER, PARAMETER :: nYRaw = nAzimRAM

  REAL(DP), PARAMETER :: mu0 = 4._dp*pi_d*1.E-7_dp, BEarth = 0.31_dp*1.E-4_dp, REarth = 6.4_dp*1.E6_dp 

  REAL(SP) :: rHour
  REAL(DP) :: pressurequot

  INTEGER :: method=2, tsygcorrect, isSorDetailNeeded, isEnergDetailNeeded, isFBDetailNeeded, &
       nmaximum, iCorrectedPressure, iPressureChoice, iAMR, iOutput, iOuterMethod, lconv,inest,nitry, &
       isotropy, iCountPressureCall = 0, iUseSavedData, iteration, isw1, itout, numit, &
       iderr, itooff, jorgn, icurm, npsim,npsip, nthem,nthep,nthe2,nThetaEquator,  &
       nZetaMidnight,nzetm,nzetp,nimax,isym,nisave,isw2, iReduceAnisotropy, kMax

  REAL(DP) :: xEqMidNew(npsi), &
       radEqMidNew(npsi), g(npsi), f(npsi), fp(npsi), fzet(nzeta+1), fzetp(nzeta+1), psival(npsi), &
       alphaval(nzeta+1), alphaValInitial(nzeta+1), chival(nthe), thetaVal(nthe), &
       rhoVal(npsi), zetaVal(nzeta), &
       dela(nzeta+1), factor(nzeta+1)

  REAL(DP), DIMENSION(npsi,nzeta+1), SAVE :: PhiIono ! SAVE attribute to retain it if SW data is bad (for Weimer model)

  REAL(DP), DIMENSION(npsi,nzeta+1) :: dPhiIonodAlpha, dPhiIonodbeta

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jacobian, gradRhoSq, gradZetaSq, &
       gradRhoGradZeta, gradRhoGradTheta, gradThetaGradZeta, &
       curvaturePsi, gradPsiGradAlpha, curvatureGradPsi, &
       gradAlphaSq, radialCurvature

  REAL(DP), DIMENSION(nthe,npsi,nzeta+1) :: ppar, pper, &
       sigma, tau, pressure3D, EXInd, EYInd, EZInd, EXConv, EYConv, EZConv, &
       psi,psiSav1,psiSav2,alfa,alfaSav1,alfaSav2, x,y,z, xx, yy, &
       bsq, bf, bfInitial, bpar,bper, &
       phij, dpdAlpha, dpdPsi, dSqPdAlphaSq, dSqPdPsiSq, &
       bj, jParDirect

  REAL(dp), DIMENSION(nthe,npsi,nzeta) :: dPperdAlpha, dBBdAlpha, &
       dBBdPsi, dPperdPsi, dPpardAlpha, dPpardPsi,dPperdTheta, dPperdRho, dPperdZeta, &
       dBsqdAlpha, dBsqdPsi, bX, bY, bZ, bxintern, byintern, bzintern, dBsqdTheta, dAlphadT, dBetadT

  REAL(DP), DIMENSION(152,202,48) :: alpha_cyl, beta_cyl

  REAL(DP) :: dphi, rjac, mjac, aisw2, sumdb, sumb, diffmx, xmx, &
       dr,dt,dpPrime,rdr,rdt,rdp,rdrsq,rdtsq,rdpsq,  &
       rdr2,rdt2,rdp2,rdr4,rdt4,rdp4,rdrdp4,rdpdt4,rdtdr4 ,omega,delta,eps,  &
       xmin,xmaj,epsb, xzero3,re,aguess, bzero, xzero,xplmin,xplmax, &
       xpsiout,xpsiin,bcfactor,bzimf, psiout,psiin,psitot,psimax,psimin,psimm, &
       pnormal, bnormal, enormal, pjconst, constz, constTheta

  REAL(DP) :: r0Start, blendInitial, blendGlobal, blendGlobalInitial, decreaseConv, &
       decreaseConvAlphaMin=5e-1, decreaseConvAlphaMax=5e-1, decreaseConvPsiMin=5e-1, &
       decreaseConvPsiMax=5e-1, blendAlphaInit=0.20, blendPsiInit=0.20, &
       decreaseConvAlpha, decreaseConvPsi, blendAlpha, blendPsi, &
       thresh=1.1, damp=0.9, relax=1.2, blendMin=0.01, blendMax=0.5 
  
  REAL(DP) :: start_time, end_time
  REAL(DP) :: DstDPS, DstDPSInsideGeo, DstBiot, DstBiotInsideGeo

  REAL(DP) :: ratioEq(npsi,nzeta)

  CHARACTER*2 :: Day, DayP1, HourInteger, HourIntegerP1, MinChar
  CHARACTER*5 :: HourDecimal, HourDecShort
  CHARACTER*4 :: suffixFilePress
  CHARACTER*500 :: fileNamePressure
  CHARACTER*200 :: prefixIn, prefixOut, prefixRAMOut, prefixPres
  INTEGER :: iDay, iHour, iHourAbs, iMin, iHourBegin, iHourChoice, iCallRoutine, iWantAlphaExtrapolation, iAzimOffset, &
       iLossCone, k9, k21, iInducedE, iConvE

  INTEGER :: numProc, rank, iTilt
  INTEGER :: iConvGlobal, iSm, iSm2, iSWMF, iAlphaMove, iPsiMove, iST3, iInterpMethod

  REAL(DP) :: p_dyn = 0._dp, by_imf = 0._dp, bz_imf = 0._dp, dst_global = 0._dp, wTsyg(1:6) = 0._dp, tilt = 0._dp ! Initialize

  REAL(DP) :: normDiff, normDiffPrev, normJxB, normGradP

  CHARACTER(LEN=2) :: Hour
  CHARACTER(LEN=6) :: event
  REAL(DP) :: pdynGlobal, byimfGlobal, bzimfGlobal
  INTEGER :: iCompDomain, iDumpRAMFlux, iConvergenceMethod, iElectric = 0

  INTEGER :: nrelax
  REAL(DP) :: radGrid(npsi,nzeta), angleGrid(npsi,nzeta)
  REAL(DP) :: errorAlpha, errorAlphaPrev, errorPsi, errorPsiPrev

  ! Persistent file units:
  INTEGER :: iUnitDst, iUnitLog

END MODULE Module1
