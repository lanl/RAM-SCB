!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamRestart

  implicit none
  
  contains
  !==========================================================================
  subroutine write_restart
    use ModRamParams, ONLY: TimedRestarts
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRestartOut, niter
    use ModRamFunctions, ONLY: RamFileName
    use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamStart, TimeRamNow, DtsNext, &
                               TOld
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP
    use ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap
    use ModScbParams,    ONLY: constZ, constTheta
    use ModScbVariables, ONLY: x, y, z, alphaVal, psiVal, chiVal, xpsiout, &
                               xpsiin, kmax, thetaVal, f, fzet, &
                               zetaVal, rhoVal
    !!!! Module Subroutines/Functions
    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts
    !!!! Share Modules
    use ModIOUnit, ONLY: UNITTMP_
    !!!! NetCdf Modules
    use netcdf

    implicit none
    
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iPParTVar, &
               iPPerTVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, &
               iTOldVar, iCValVar, icTVar, icZVar, ipsiInVar, ipsiOutVar, &
               ikMaxVar, idAdRVar, idBdPVar, iThetaVar, iRhoVar, iZetaVar

    integer :: nRDim, nTDim, nEDim, nPaDim, nSDim, nThetaDim, nPsiDim, &
               nZetaDim, nRPDim, nZetaPDim
    integer, parameter :: iDeflate = 2, yDeflate = 1

    character(len=200) :: NameFile

    character(len=*), parameter :: NameSub='write_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    if(DoTest)write(*,'(a,f11.2)') 'RAM-SCB: Writing restarts at t=',&
         TimeRamElapsed

    ! Write ascii portion of restart.
    if (TimedRestarts) then
       NameFile=RamFileName(PathRestartOut//'/restart_info','txt',TimeRamNow)
    else
       NameFile=PathRestartOut//'/restart_info.txt'
    endif
    open(unit=UnitTMP_, file=trim(NameFile), status='replace')
    write(UnitTMP_, *) 'TIMING:'
    write(UnitTMP_, '(a, i4.4, 2i2.2, 1x, 3i2.2)')'Start (YYYYMMDD HHMMSS)= ', &
         TimeRamStart%iYear, TimeRamStart%iMonth, TimeRamStart%iDay, &
         TimeRamStart%iHour, TimeRamStart%iMinute,TimeRamStart%iSecond
    write(UnitTMP_, '(a,f15.4)')'Elapsed (seconds)      = ', TimeRamElapsed
    write(UnitTMP_, '(a, i15)') 'Completed Iterations   = ', nIter
    write(UnitTMP_, *)'GRID:'
    write(UnitTMP_, '(a, 4i3)') 'nR, nL, nE, nPA        = ', nR, nT, nE, nPA
    close(unitTMP_)

    ! OPEN FILE
    if (TimedRestarts) then
       NameFile = RamFileName(PathRestartOut//'/restart','nc',TimeRamNow)
    else
       NameFile = PathRestartOut//'/restart.nc'
    endif
    iStatus = nf90_create(trim(NameFile), nf90_HDF5, iFileID)
    call ncdf_check(iStatus, NameSub)
    call write_ncdf_globatts(iFileID)

    ! CREATE DIMENSIONS
    iStatus = nf90_def_dim(iFileID, 'nR',     nR,     nRDim)
    iStatus = nf90_def_dim(iFileID, 'nRP',    nR+1,   nRPDim)
    iStatus = nf90_def_dim(iFileID, 'nT',     nT,     nTDim)
    iStatus = nf90_def_dim(iFileID, 'nE',     nE,     nEDim)
    iStatus = nf90_def_dim(iFileID, 'nPa',    nPa,    nPaDim)
    iStatus = nf90_def_dim(iFileID, 'nS',     4,      nSDim)
    iStatus = nf90_def_dim(iFileID, 'nTheta', nthe,   nThetaDim)
    iStatus = nf90_def_dim(iFileID, 'nPsi',   npsi,   nPsiDim)
    iStatus = nf90_def_dim(iFileID, 'nZeta',  nzeta,  nZetaDim)
    iStatus = nf90_def_dim(iFileID, 'nZetaP', nzetap, nZetaPDim)

    ! START DEFINE MODE
    !! FLUXES
    !!! Electron Flux
    iStatus = nf90_def_var(iFileID, 'FluxE', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxEVar)
    iStatus = nf90_def_var_deflate(iFileID, iFluxEVar, 0, yDeflate, iDeflate)

    iStatus = nf90_put_att(iFileID, iFluxEVar, 'title', &
                           'This is an example title')

    !!! Proton Flux
    iStatus = nf90_def_var(iFileID, 'FluxH', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxHVar)
    iStatus = nf90_def_var_deflate(iFileID, iFluxHVar, 0, yDeflate, iDeflate)

    !!! Oxygen Flux
    iStatus = nf90_def_var(iFileID, 'FluxO', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxOVar)
    iStatus = nf90_def_var_deflate(iFileID, iFluxOVar, 0, yDeflate, iDeflate)

    !!! Helium Flux
    iStatus = nf90_def_var(iFileID, 'FluxHe', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxHeVar)
    iStatus = nf90_def_var_deflate(iFileID, iFluxHeVar, 0, yDeflate, iDeflate)

    !!! Boundary Flux
    !iStatus = nf90_def_var(iFileID, 'FGEOS', nf90_double, &
    !                       (/nSdim,nTDim,nEDim,nPadim/), iGEOVar)
    !iStatus = nf90_def_var_deflate(iFileID, iGEOVar, 0, 1, iDeflate)

    !!! 3D Flux (needed for satellite files)
    !iStatus = nf90_def_var(iFileID, 'Flux3D', nf90_double, &
    !                       (/nSdim,nPsiDim,nZetaDim,nEDim,nPADim/), iFlux3DVar)
    !iStatus = nf90_def_var_deflate(iFileID, iFlux3DVar, 0, 1, iDeflate)

    !! PRESSURES
    !!! Total Parallel Pressures
    iStatus = nf90_def_var(iFileID, 'PParT', nf90_double, &
                           (/nSdim,nRDim,nTDim/), iPParTVar)
    iStatus = nf90_def_var_deflate(iFileID, iPParTVar, 0, yDeflate, iDeflate)

    !!! Total Perpendicular Pressures
    iStatus = nf90_def_var(iFileID, 'PPerT', nf90_double, &
                           (/nSdim,nRDim,nTDim/), iPPerTVar)
    iStatus = nf90_def_var_deflate(iFileID, iPPerTVar, 0, yDeflate, iDeflate)

    !! MAGNETIC FIELD
    !!! Bx
    !iStatus = nf90_def_var(iFileID, 'Bx', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaDim/), iBxVar)
    !iStatus = nf90_def_var_deflate(iFileID, iBxVar, 0, 1, iDeflate)

    !!! By
    !iStatus = nf90_def_var(iFileID, 'By', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaDim/), iByVar)
    !iStatus = nf90_def_var_deflate(iFileID, iByVar, 0, 1, iDeflate)

    !!! Bz
    !iStatus = nf90_def_var(iFileID, 'Bz', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaDim/), iBzVar)
    !iStatus = nf90_def_var_deflate(iFileID, iBzVar, 0, 1, iDeflate)

    !!! Total B
    !iStatus = nf90_def_var(iFileID, 'B', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaPDim/), iBfVar)
    !iStatus = nf90_def_var_deflate(iFileID, iBfVar, 0, 1, iDeflate)

    !! ELECTRIC FIELD/POTENTIAL
    !!! Previous Electric Potential
    iStatus = nf90_def_var(iFileID, 'VTOL', nf90_double, &
                           (/nRPDim,nTDim/), iVTOLVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTOLVar, 0, yDeflate, iDeflate)

    !!! Next Electric Potential
    iStatus = nf90_def_var(iFileID, 'VTN', nf90_double, &
                           (/nRPDim,nTDim/), iVTNVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTNVar, 0, yDeflate, iDeflate)

    !!! Current Electric Potential
    iStatus = nf90_def_var(iFileID, 'VT', nf90_double, &
                           (/nRPDim,nTDim/), iVTVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTVar, 0, yDeflate, iDeflate)


    !! hI OUTPUTS
    !!! H
    iStatus = nf90_def_var(iFileID, 'FNHS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iHVar)
    iStatus = nf90_def_var_deflate(iFileID, iHVar, 0, yDeflate, iDeflate)

    !!! Hbn
    iStatus = nf90_def_var(iFileID, 'BOUNHS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iBHVar)
    iStatus = nf90_def_var_deflate(iFileID, iBHVar, 0, yDeflate, iDeflate)

    !!! I
    iStatus = nf90_def_var(iFileID, 'FNIS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iIVar)
    iStatus = nf90_def_var_deflate(iFileID, iIVar, 0, yDeflate, iDeflate)

    !!! Ibn
    iStatus = nf90_def_Var(iFileID, 'BOUNIS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iBIVar)
    iStatus = nf90_def_var_deflate(iFileID, iBIVar, 0, yDeflate, iDeflate)

    !!! BNES
    iStatus = nf90_def_var(iFileID, 'BNES', nf90_double, &
                           (/nrPDim,nTDim/), iBNESVar)
    iStatus = nf90_def_var_deflate(iFileID, iBNESVar, 0, yDeflate, iDeflate)

    !!! Neutral Density
    iStatus = nf90_def_var(iFileID, 'HDNS', nf90_double, &
                           (/nRPDim,nTDim,nPaDim/), iHDNSVar)
    iStatus = nf90_def_var_deflate(iFileID, iHDNSVar, 0, yDeflate, iDeflate)

    !!!
    iStatus = nf90_def_var(iFileID, 'EIR', nf90_double, &
                           (/nrPDim,nTDim/), iEIRVar)
    iStatus = nf90_def_var_deflate(iFileID, iEIRVar, 0, yDeflate, iDeflate)

    !!!
    iStatus = nf90_def_var(iFileID, 'EIP', nf90_double, &
                           (/nrPDim,nTDim/), iEIPVar)
    iStatus = nf90_def_var_deflate(iFileID, iEIPVar, 0, yDeflate, iDeflate)

    !!! db/dt
    iStatus = nf90_def_var(iFileID, 'dBdt', nf90_double, &
                           (/nrPDim,nTDim/), iDBDTVar)
    iStatus = nf90_def_var_deflate(iFileID, iDBDTVar, 0, yDeflate, iDeflate)

    !!! dI/dt
    iStatus = nf90_def_var(iFileID, 'dIdt', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iDIDTVar)
    iStatus = nf90_def_var_deflate(iFileID, iDIDTVar, 0, yDeflate, iDeflate)

    !!! dIbn/dt
    iStatus = nf90_def_var(iFileID, 'dIbndt', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iDIBNVar)
    iStatus = nf90_def_var_deflate(iFileID, iDIBNVar, 0, yDeflate, iDeflate)

    !! GRID OUTPUTS
    !!! X
    iStatus = nf90_def_var(iFileID, 'x', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iXVar)
    iStatus = nf90_def_var_deflate(iFileID, iXVar, 0, yDeflate, iDeflate)

    !!! Y
    iStatus = nf90_def_var(iFileID, 'y', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iYVar)
    iStatus = nf90_def_var_deflate(iFileID, iYVar, 0, yDeflate, iDeflate)

    !!! Z
    iStatus = nf90_def_var(iFileID, 'z', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iZVar)
    iStatus = nf90_def_var_deflate(iFileID, iZVar, 0, yDeflate, iDeflate)

    !!! Theta (Along the field line)
    iStatus = nf90_def_var(iFileID, 'theta', nf90_double, &
                           (/nThetaDim/), iThetaVar)
    iStatus = nf90_def_var_deflate(iFileID, iThetaVar, 0, yDeflate, iDeflate)

    !!! Rho (Radial)
    iStatus = nf90_def_var(iFileID, 'rho', nf90_double, &
                           (/nPsiDim/), iRhoVar)
    iStatus = nf90_def_var_deflate(iFileID, iRhoVar, 0, yDeflate, iDeflate)

    !!! Zeta (Azimuthal)
    iStatus = nf90_def_var(iFileID, 'zeta', nf90_double, &
                           (/nZetaDim/), iZetaVar)
    iStatus = nf90_def_var_deflate(iFileID, iZetaVar, 0, yDeflate, iDeflate)

    !! ALPHA/BETA
    !!! AlphaVal (psiVal)
    iStatus = nf90_def_var(iFileID, 'alphaVal', nf90_double, &
                           (/nPsiDim/), iAValVar)
    iStatus = nf90_def_var_deflate(iFileID, iAValVar, 0, yDeflate, iDeflate)

    !!! BetaVal (alphaVal)
    iStatus = nf90_def_var(iFileID, 'betaVal', nf90_double, &
                           (/nZetaPDim/), iBValVar)
    iStatus = nf90_def_var_deflate(iFileID, iBValVar, 0, yDeflate, iDeflate)

    !!! ChiVal
    iStatus = nf90_def_var(iFileID, 'chiVal', nf90_double, &
                           (/nThetaDim/), iCValVar)
    iStatus = nf90_def_var_deflate(iFileID, iCValVar, 0, yDeflate, iDeflate)

    !!! Alpha (psi)
    !iStatus = nf90_def_var(iFileID, 'alpha', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaPDim/), iAlphaVar)
    !iStatus = nf90_def_var_deflate(iFileID, iAlphaVar, 0, 1, iDeflate)

    !!! Beta (alpha)
    !iStatus = nf90_def_var(iFileID, 'beta', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaPDim/), iBetaVar)
    !iStatus = nf90_def_var_deflate(iFileID, iBetaVar, 0, 1, iDeflate)

    !!! Chi
    !iStatus = nf90_def_var(iFileID, 'chi', nf90_double, &
    !                       (/nThetaDim,nPsiDim,nZetaPDim/), iChiVar)
    !iStatus = nf90_def_var_deflate(iFileID, iChiVar, 0, 1, iDeflate)

    !!! dAlpha/dRho (f)
    iStatus = nf90_def_var(iFileID, 'dAlphadRho', nf90_double, &
                           (/nPsiDim/), idAdRVar)
    iStatus = nf90_def_var_deflate(iFileID, idAdRVar, 0, yDeflate, iDeflate)

    !!! dBeta/dPhi (fzet)
    iStatus = nf90_def_var(iFileID, 'dBetadPhi', nf90_double, &
                           (/nZetaPDim/), idBdPVar)
    iStatus = nf90_def_var_deflate(iFileID, idBdPVar, 0, yDeflate, iDeflate)

    !! MISC
    !!! Index of pitch angles (needed for satellite files)
    !iStatus = nf90_def_var(iFileID, 'IndexPA', nf90_int, &
    !                       (/nThetaDim,nPsiDim,nZetaDim,nPaDim/), iPaIndexVar)
    !iStatus = nf90_def_var_deflate(iFileID, iPaIndexVar, 0, 1, iDeflate)

    !!! Previous calculated time step
    iStatus = nf90_def_var(iFileID, 'DtNext', nf90_double, iDtVar)

    !!! Time when SCB was last called
    iStatus = nf90_def_var(iFileID, 'TOld', nf90_double, iTOldVar)

    !!!
    iStatus = nf90_def_var(iFileID, 'constTheta', nf90_double, icTVar)
    iStatus = nf90_def_var(iFileID, 'constZ',     nf90_double, icZVar)
    iStatus = nf90_def_var(iFileID, 'xpsiin',     nf90_double, ipsiInVar)
    iStatus = nf90_def_var(iFileID, 'xpsiout',    nf90_double, ipsiOutVar)
    iStatus = nf90_def_var(iFileID, 'kmax',       nf90_int,    ikMaxVar)

    ! END DEFINE MODE
    iStatus = nf90_enddef(iFileID)
    call ncdf_check(iStatus, NameSub)


    ! START WRITE MODE
    !! FLUXES
    iStatus = nf90_put_var(iFileID, iFluxEVar,  F2(1,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHVar,  F2(2,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHeVar, F2(3,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxOVar,  F2(4,:,:,:,:))
    !iStatus = nf90_put_var(iFileID, iGEOVar,    FGEOS(:,:,:,:))
    !iStatus = nf90_put_var(iFileID, iFlux3DVar, FLUX3DEQ(:,:,:,:,:))

    !! PRESSURES
    iStatus = nf90_put_var(iFileID, iPParTVar, PParT(:,:,:))
    iStatus = nf90_put_var(iFileID, iPPerTVar, PPerT(:,:,:))

    !! MAGNETIC FIELD
    !iStatus = nf90_put_var(iFileID, iBxVar, bX(:,:,:))
    !iStatus = nf90_put_var(iFileID, iByVar, bY(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBzVar, bZ(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBfVar, bf(:,:,:))

    !! ELECTRIC FIELD
    iStatus = nf90_put_var(iFileID, iVTOLVar, VTOL(:,:))
    iStatus = nf90_put_var(iFileID, iVTNVar,  VTN(:,:))
    iStatus = nf90_put_var(iFileID, iVTVar,   VT(:,:))

    !! hI OUTPUTS
    iStatus = nf90_put_var(iFileID, iHVar,    FNHS(:,:,:))
    iStatus = nf90_put_var(iFileID, iBHVar,   BOUNHS(:,:,:))
    iStatus = nf90_put_var(iFileID, iIVar,    FNIS(:,:,:))
    iStatus = nf90_put_var(iFileID, iBIVar,   BOUNIS(:,:,:))
    iStatus = nf90_put_var(iFileID, iBNESVar, BNES(:,:))
    iStatus = nf90_put_var(iFileID, iHDNSVar, HDNS(:,:,:))
    iStatus = nf90_put_var(iFileID, iEIRVar,  EIR(:,:))
    iStatus = nf90_put_var(iFileID, iEIPVar,  EIP(:,:))
    iStatus = nf90_put_var(iFileID, iDBDTVar, dBdt(:,:))
    iStatus = nf90_put_var(iFileID, iDIDTVar, dIdt(:,:,:))
    iStatus = nf90_put_var(iFileID, iDIBNVar, dIbndt(:,:,:))

    !! GRID OUTPUTS
    iStatus = nf90_put_var(iFileID, iXVar, x(:,:,:))
    iStatus = nf90_put_var(iFileID, iYVar, y(:,:,:))
    iStatus = nf90_put_var(iFileID, iZVar, z(:,:,:))
    iStatus = nf90_put_var(iFileID, iThetaVar, thetaVal(:))
    iStatus = nf90_put_var(iFileID, iRhoVar,   rhoVal(:))
    iStatus = nf90_put_var(iFileID, iZetaVar,  zetaVal(:))

    !! ALPHA/BETA
    iStatus = nf90_put_var(iFileID, iAValVar,  psival(:))
    iStatus = nf90_put_var(iFileID, iBValVar,  alphaval(:))
    iStatus = nf90_put_var(iFileID, iCValVar,  chiVal(:))
    !iStatus = nf90_put_var(iFileID, iAlphaVar, psi(:,:,:))
    !iStatus = nf90_put_var(iFileID, iBetaVar,  alfa(:,:,:))
    !iStatus = nf90_put_var(iFileID, iChiVar,   chi(:,:,:))
    iStatus = nf90_put_var(iFileID, idAdRVar,  f(:))
    iStatus = nf90_put_var(iFileID, idBdPVar,  fzet(:))

    !! MISC
    !iStatus = nf90_put_var(iFileID, iPaIndexVar, indexPA(:,:,:,:))
    iStatus = nf90_put_var(iFileID, iDtVar, DtsNext)
    iStatus = nf90_put_var(iFileID, iTOldVar, TOld)
    iStatus = nf90_put_var(iFileID, icTVar, constTheta)
    iStatus = nf90_put_var(iFileID, icZVar, constZ)
    iStatus = nf90_put_var(iFileID, ipsiInVar, xpsiin)
    iStatus = nf90_put_var(iFileID, ipsiOutVar, xpsiout)
    iStatus = nf90_put_var(iFileID, ikmaxVar, kmax)

    ! END WRITE MODE
    call ncdf_check(iStatus, NameSub)
    
    ! CLOSE FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_restart

  !==========================================================================
  subroutine read_restart
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, PathRestartIn
    use ModRamFunctions, ONLY: RamFileName
    use ModRamTiming,    ONLY: DtsNext, TOld
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP, PParH, PPerH, PParO, &
                               PPerO, PParHe, PPerHe, PParE, PPerE
    use ModScbGrids,     ONLY: npsi
    use ModScbParams,    ONLY: constTheta, constZ
    use ModScbVariables, ONLY: x, y, z, alphaVal, psiVal, chiVal, xpsiout, &
                               xpsiin, left, right, kmax, thetaVal, f, fp, &
                               fzet, fzetp, zetaVal, rhoVal
    !!!! Module Subroutines/Functions
    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts
    !!!! NetCdf Modules
    use netcdf

    implicit none
    
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iPParTVar, &
               iPPerTVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, &
               iTOldVar, iCValVar, icTVar, icZVar, ipsiInVar, ipsiOutVar, &
               ikMaxVar, idAdRVar, idBdPVar, iThetaVar, iRhoVar, iZetaVar

    character(len=*), parameter :: NameSub='read_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
   
    left = 1
    right = npsi

    write(*,*) 'Loading restart files.'
    ! LOAD RESTART FILE
    iStatus = nf90_open(trim(PathRestartIn//'/restart.nc'),nf90_nowrite,iFileID)
    call ncdf_check(iStatus, NameSub)

    ! GET VARIABLE IDS
    !! FLUXES
    iStatus = nf90_inq_varid(iFileID, 'FluxE',  iFluxEVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxH',  iFluxHVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxHe', iFluxHeVar)
    iStatus = nf90_inq_varid(iFileID, 'FluxO',  iFluxOVar)
    !iStatus = nf90_inq_varid(iFileID, 'FGEOS',  iGEOVar)
    !iStatus = nf90_inq_varid(iFileID, 'Flux3D', iFlux3DVar)

    !! PRESSURES
    iStatus = nf90_inq_varid(iFileID, 'PParT', iPParTVar)
    iStatus = nf90_inq_varid(iFileID, 'PPerT', iPPerTVar)

    !! MAGNETIC FIELD
    !iStatus = nf90_inq_varid(iFileID, 'Bx', iBxVar)
    !iStatus = nf90_inq_varid(iFileID, 'By', iByVar)
    !iStatus = nf90_inq_varid(iFileID, 'Bz', iBzVar)
    !iStatus = nf90_inq_varid(iFileID, 'B',  iBfVar)

    !! ELECTRIC FIELD
    iStatus = nf90_inq_varid(iFileID, 'VTOL', iVTOLVar)
    iStatus = nf90_inq_varid(iFileID, 'VTN',  iVTNVar)
    iStatus = nf90_inq_varid(iFileID, 'VT',   iVTVar)

    !! hI OUTPUTS
    iStatus = nf90_inq_varid(iFileID, 'FNHS',   iHVar)
    iStatus = nf90_inq_varid(iFileID, 'BOUNHS', iBHVar)
    iStatus = nf90_inq_varid(iFileID, 'FNIS',   iIvar)
    iStatus = nf90_inq_varid(iFileID, 'BOUNIS', iBIVar)
    iStatus = nf90_inq_varid(iFileID, 'BNES',   iBNESVar)
    iStatus = nf90_inq_varid(iFileiD, 'HDNS',   iHDNSVar)
    iStatus = nf90_inq_varid(iFileID, 'EIR',    iEIRVar)
    iStatus = nf90_inq_varid(iFileID, 'EIP',    iEIPVar)
    iStatus = nf90_inq_varid(iFileID, 'dBdt',   iDBDTVar)
    iStatus = nf90_inq_varid(iFileID, 'dIdt',   iDIDTVar)
    iStatus = nf90_inq_varid(iFileID, 'dIbndt', iDIBNVar)

    !! GRID OUTPUTS
    iStatus = nf90_inq_varid(iFileID, 'x', iXVar)
    iStatus = nf90_inq_varid(iFileID, 'y', iYVar)
    iStatus = nf90_inq_varid(iFileID, 'z', iZVar)
    iStatus = nf90_inq_varid(iFileID, 'theta', iThetaVar)
    iStatus = nf90_inq_varid(iFileID, 'rho',   iRhoVar)
    iStatus = nf90_inq_varid(iFileID, 'zeta',  iZetaVar)

    !! ALPHA/BETA
    iStatus = nf90_inq_varid(iFileID, 'alphaVal', iAValVar)
    iStatus = nf90_inq_varid(iFileID, 'betaVal',  iBValVar)
    iStatus = nf90_inq_varid(iFileID, 'chiVal',   iCValVar)
    !iStatus = nf90_inq_varid(iFileID, 'alpha',    iAlphaVar)
    !iStatus = nf90_inq_varid(iFileID, 'beta',     iBetaVar)
    !iStatus = nf90_inq_varid(iFileID, 'chi',      iChiVar)
    iStatus = nf90_inq_varid(iFileID, 'dAlphadRho', idAdRVar)
    iStatus = nf90_inq_varid(iFileID, 'dBetadPhi',  idBdPVar)

    !! MISC
    !iStatus = nf90_inq_varid(iFileID, 'IndexPA', iPaIndexVar)
    iStatus = nf90_inq_varid(iFileID, 'DtNext', iDtVar)
    iStatus = nf90_inq_varid(iFileID, 'TOld', iTOldVar)
    iStatus = nf90_inq_varid(iFileID, 'constTheta', icTVar)
    iStatus = nf90_inq_varid(iFileID, 'constZ',     icZVar)
    iStatus = nf90_inq_varid(iFileID, 'xpsiin',     ipsiInVar)
    iStatus = nf90_inq_varid(iFileID, 'xpsiout',    ipsiOutVar)
    iStatus = nf90_inq_varid(iFileID, 'kmax',       ikmaxVar)

    ! READ DATA
    !! FLUXES
    iStatus = nf90_get_var(iFileID, iFluxEVar,  F2(1,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHVar,  F2(2,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHeVar, F2(3,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxOVar,  F2(4,:,:,:,:))
    !iStatus = nf90_get_var(iFileID, iGEOVar,    FGEOS(:,:,:,:))
    !iStatus = nf90_get_var(iFileID, iFlux3DVar, FLUX3DEQ(:,:,:,:,:))

    !! PRESSURES
    iStatus = nf90_get_var(iFileID, iPParTVar, PParT(:,:,:))
    iStatus = nf90_get_var(iFileID, iPPerTVar, PPerT(:,:,:))
    PPerO = PPerT(4,:,:)
    PParO = PParT(4,:,:)
    PPerHe = PPerT(3,:,:)
    PParHe = PParT(3,:,:)
    PPerE = PPerT(1,:,:)
    PParE = PParT(1,:,:)
    PPerH = PPerT(2,:,:)
    PParH = PParT(2,:,:)

    !! MAGNETIC FIELD
    !iStatus = nf90_get_var(iFileID, iBxVar, bX(:,:,:))
    !iStatus = nf90_get_var(iFileID, iByVar, bY(:,:,:))
    !iStatus = nf90_get_var(iFileID, iBzVar, bZ(:,:,:))
    !iStatus = nf90_get_Var(iFileID, iBfVar, bf(:,:,:))

    !! ELECTRIC FIELD
    iStatus = nf90_get_var(iFileID, iVTOLVar, VTOL(:,:))
    iStatus = nf90_get_var(iFileID, iVTNVar,  VTN(:,:))
    iStatus = nf90_get_var(iFileID, iVTVar,   VT(:,:))

    !! hI OUTPUTS
    iStatus = nf90_get_var(iFileID, iHVar,    FNHS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBHVar,   BOUNHS(:,:,:))
    iStatus = nf90_get_var(iFileID, iIVar,    FNIS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBIVar,   BOUNIS(:,:,:))
    iStatus = nf90_get_var(iFileID, iBNESVar, BNES(:,:))
    iStatus = nf90_get_var(iFileID, iHDNSVar, HDNS(:,:,:))
    iStatus = nf90_get_var(iFileID, iEIRVar,  EIR(:,:))
    iStatus = nf90_get_Var(iFileID, iEIPVar,  EIP(:,:))
    iStatus = nf90_get_var(iFileID, iDBDTVar, dBdt(:,:))
    iStatus = nf90_get_var(iFileID, iDIDTVar, dIdt(:,:,:))
    iStatus = nf90_get_var(iFileID, iDIBNVar, dIbndt(:,:,:))

    !! GRID OUTPUTS
    iStatus = nf90_get_var(iFileID, iXVar, x(:,:,:))
    iStatus = nf90_get_var(iFileID, iYVar, y(:,:,:))
    iStatus = nf90_get_var(iFileID, iZVar, z(:,:,:))
    iStatus = nf90_get_var(iFileID, iThetaVar, thetaVal(:))
    iStatus = nf90_get_var(iFileID, iRhoVar, rhoVal(:))
    iStatus = nf90_get_var(iFileID, iZetaVar, zetaVal(:))

    !! ALPHA/BETA
    iStatus = nf90_get_var(iFileID, iAValVar,  psival(:))
    iStatus = nf90_get_var(iFileID, iBValVar,  alphaval(:))
    iStatus = nf90_get_var(iFileID, iCValVar,  chiVal(:))
    !iStatus = nf90_get_var(iFileID, iAlphaVar, psi(:,:,:))
    !iStatus = nf90_get_var(iFileID, iBetaVar,  alfa(:,:,:))
    !iStatus = nf90_get_var(iFileID, iChiVar,   chi(:,:,:))
    iStatus = nf90_get_var(iFileID, idAdRVar,  f(:))
    iStatus = nf90_get_var(iFileID, idBdPVar,  fzet(:))
    fp = 0._dp
    fzetp = 0._dp

    !! MISC
    !iStatus = nf90_get_var(iFileID, iPaIndexVar, indexPA(:,:,:,:))
    iStatus = nf90_get_var(iFileID, iDtVar, DTsNext)
    iStatus = nf90_get_Var(iFileID, iTOldVar, TOld)
    iStatus = nf90_get_Var(iFileID, icTVar, constTheta)
    iStatus = nf90_get_Var(iFileID, icZVar, constZ)
    iStatus = nf90_get_Var(iFileID, ipsiInVar, xpsiin)
    iStatus = nf90_get_Var(iFileID, ipsiOutVar, xpsiout)
    iStatus = nf90_get_Var(iFileID, ikMaxVar, kMax)

    ! CLOSE RESTART FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine read_restart

!==============================================================================
END MODULE ModRamRestart
