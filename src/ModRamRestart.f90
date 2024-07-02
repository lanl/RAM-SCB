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
    use ModRamGrids,     ONLY: NR, NT, NE, NPA, nS
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP, EKEV, PA, NECR, TAU, species
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbParams,    ONLY: constZ, constTheta
    use ModScbVariables, ONLY: x, y, z, alphaVal, psiVal, chiVal, xpsiout, &
                               xpsiin, kmax, thetaVal, f, fzet, &
                               zetaVal, rhoVal, bX, bY, bZ
    !!!! Module Subroutines/Functions
    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts
    !!!! Share Modules
    use ModIOUnit, ONLY: UNITTMP_
    !!!! NetCdf Modules
    use netcdf

    implicit none
    
    integer :: iFluxVar, iPParTVar, iPPerTVar, iHVar, iBHVar, iEGridVar, iPaGridVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iDtVar, iXVar, iYVar, iZVar, iBxVar, iByVar, iBzVar, &
               iVTNVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, &
               iTOldVar, iCValVar, icTVar, icZVar, ipsiInVar, ipsiOutVar, &
               ikMaxVar, idAdRVar, idBdPVar, iThetaVar, iRhoVar, iZetaVar, &
               iNECRVar, itauVar

    integer :: nRDim, nTDim, nEDim, nPaDim, nSDim, nThetaDim, nPsiDim, &
               nZetaDim, iS
    integer, parameter :: iDeflate = 2, yDeflate = 1

    character(len=999) :: NameFile

    character(len=*), parameter :: NameSub='write_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    if(DoTest)write(*,'(a,f11.2)') 'RAM-SCB: Writing restarts at t=',&
         TimeRamElapsed

    ! Write ascii portion of restart.
    if (TimedRestarts) then
       NameFile=RamFileName(trim(PathRestartOut)//'/restart_info','txt',TimeRamNow)
    else
       NameFile=trim(PathRestartOut)//'/restart_info.txt'
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
       NameFile = RamFileName(trim(PathRestartOut)//'/restart','nc',TimeRamNow)
    else
       NameFile = trim(PathRestartOut)//'/restart.nc'
    endif
    iStatus = nf90_create(trim(NameFile), nf90_HDF5, iFileID)
    call ncdf_check(iStatus, NameSub)
    call write_ncdf_globatts(iFileID)

    ! CREATE DIMENSIONS
    iStatus = nf90_def_dim(iFileID, 'nR',     nR,     nRDim)
    iStatus = nf90_def_dim(iFileID, 'nT',     nT,     nTDim)
    iStatus = nf90_def_dim(iFileID, 'nE',     nE,     nEDim)
    iStatus = nf90_def_dim(iFileID, 'nPa',    nPa,    nPaDim)
    iStatus = nf90_def_dim(iFileID, 'nS',     4,      nSDim)
    iStatus = nf90_def_dim(iFileID, 'nTheta', nthe,   nThetaDim)
    iStatus = nf90_def_dim(iFileID, 'nPsi',   npsi,   nPsiDim)
    iStatus = nf90_def_dim(iFileID, 'nZeta',  nzeta,  nZetaDim)

    ! START DEFINE MODE
    !! TEMP STUFF
    do iS = 1, nS
        iStatus = nf90_def_var(iFileID, 'EnergyGrid_'//species(iS)%s_name, nf90_double, &
                               (/nEDim/), iEgridVar)
    enddo
    iStatus = nf90_def_var(iFileID, 'PitchAngleGrid', nf90_double, &
                           (/nPaDim/), iPaGridVar)
    !! FLUXES
    do iS = 1, nS
        iStatus = nf90_def_var(iFileID, 'Flux_'//species(iS)%s_name, nf90_double, &
                               (/nRDim, nTDim, nEDim, nPaDim/), iFluxVar)
        iStatus = nf90_def_var_deflate(iFileID, iFluxVar, 0, yDeflate, iDeflate)
    enddo

    !! PRESSURES
    !!! Total Parallel Pressures
    do iS = 1, nS
        !!! Total Parallel Pressures
        iStatus = nf90_def_var(iFileID, 'PParT_'//species(iS)%s_name, nf90_double, &
                               (/nRDim,nTDim/), iPParTVar)
        iStatus = nf90_def_var_deflate(iFileID, iPParTVar, 0, yDeflate, iDeflate)

        !!! Total Perpendicular Pressures
        iStatus = nf90_def_var(iFileID, 'PPerT_'//species(iS)%s_name, nf90_double, &
                               (/nRDim,nTDim/), iPPerTVar)
        iStatus = nf90_def_var_deflate(iFileID, iPPerTVar, 0, yDeflate, iDeflate)
    enddo

    !! MAGNETIC FIELD
    !! Bx
    iStatus = nf90_def_var(iFileID, 'Bx', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iBxVar)
    iStatus = nf90_def_var_deflate(iFileID, iBxVar, 0, 1, iDeflate)

    !! By
    iStatus = nf90_def_var(iFileID, 'By', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iByVar)
    iStatus = nf90_def_var_deflate(iFileID, iByVar, 0, 1, iDeflate)

    !! Bz
    iStatus = nf90_def_var(iFileID, 'Bz', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iBzVar)
    iStatus = nf90_def_var_deflate(iFileID, iBzVar, 0, 1, iDeflate)

    !! ELECTRIC FIELD/POTENTIAL
    !!! Previous Electric Potential
    iStatus = nf90_def_var(iFileID, 'VTOL', nf90_double, &
                           (/nRDim,nTDim/), iVTOLVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTOLVar, 0, yDeflate, iDeflate)

    !!! Next Electric Potential
    iStatus = nf90_def_var(iFileID, 'VTN', nf90_double, &
                           (/nRDim,nTDim/), iVTNVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTNVar, 0, yDeflate, iDeflate)

    !!! Current Electric Potential
    iStatus = nf90_def_var(iFileID, 'VT', nf90_double, &
                           (/nRDim,nTDim/), iVTVar)
    iStatus = nf90_def_var_deflate(iFileID, iVTVar, 0, yDeflate, iDeflate)


    !! hI OUTPUTS
    !!! H
    iStatus = nf90_def_var(iFileID, 'FNHS', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iHVar)
    iStatus = nf90_def_var_deflate(iFileID, iHVar, 0, yDeflate, iDeflate)

    !!! Hbn
    iStatus = nf90_def_var(iFileID, 'BOUNHS', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iBHVar)
    iStatus = nf90_def_var_deflate(iFileID, iBHVar, 0, yDeflate, iDeflate)

    !!! I
    iStatus = nf90_def_var(iFileID, 'FNIS', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iIVar)
    iStatus = nf90_def_var_deflate(iFileID, iIVar, 0, yDeflate, iDeflate)

    !!! Ibn
    iStatus = nf90_def_Var(iFileID, 'BOUNIS', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iBIVar)
    iStatus = nf90_def_var_deflate(iFileID, iBIVar, 0, yDeflate, iDeflate)

    !!! BNES
    iStatus = nf90_def_var(iFileID, 'BNES', nf90_double, &
                           (/nRDim,nTDim/), iBNESVar)
    iStatus = nf90_def_var_deflate(iFileID, iBNESVar, 0, yDeflate, iDeflate)

    !!! Neutral Density
    iStatus = nf90_def_var(iFileID, 'HDNS', nf90_double, &
                           (/nRDim,nTDim,nPaDim/), iHDNSVar)
    iStatus = nf90_def_var_deflate(iFileID, iHDNSVar, 0, yDeflate, iDeflate)

    !!!
    iStatus = nf90_def_var(iFileID, 'EIR', nf90_double, &
                           (/nRDim,nTDim/), iEIRVar)
    iStatus = nf90_def_var_deflate(iFileID, iEIRVar, 0, yDeflate, iDeflate)

    !!!
    iStatus = nf90_def_var(iFileID, 'EIP', nf90_double, &
                           (/nRDim,nTDim/), iEIPVar)
    iStatus = nf90_def_var_deflate(iFileID, iEIPVar, 0, yDeflate, iDeflate)

    !!! db/dt
    iStatus = nf90_def_var(iFileID, 'dBdt', nf90_double, &
                           (/nRDim,nTDim/), iDBDTVar)
    iStatus = nf90_def_var_deflate(iFileID, iDBDTVar, 0, yDeflate, iDeflate)

    !!! dI/dt
    iStatus = nf90_def_var(iFileID, 'dIdt', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iDIDTVar)
    iStatus = nf90_def_var_deflate(iFileID, iDIDTVar, 0, yDeflate, iDeflate)

    !!! dIbn/dt
    iStatus = nf90_def_var(iFileID, 'dIbndt', nf90_double, &
                           (/nRDim,nTDim,NPaDim/), iDIBNVar)
    iStatus = nf90_def_var_deflate(iFileID, iDIBNVar, 0, yDeflate, iDeflate)

    !! GRID OUTPUTS
    !!! X
    iStatus = nf90_def_var(iFileID, 'x', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iXVar)
    iStatus = nf90_def_var_deflate(iFileID, iXVar, 0, yDeflate, iDeflate)

    !!! Y
    iStatus = nf90_def_var(iFileID, 'y', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iYVar)
    iStatus = nf90_def_var_deflate(iFileID, iYVar, 0, yDeflate, iDeflate)

    !!! Z
    iStatus = nf90_def_var(iFileID, 'z', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iZVar)
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
                           (/nZetaDim/), iBValVar)
    iStatus = nf90_def_var_deflate(iFileID, iBValVar, 0, yDeflate, iDeflate)

    !!! ChiVal
    iStatus = nf90_def_var(iFileID, 'chiVal', nf90_double, &
                           (/nThetaDim/), iCValVar)
    iStatus = nf90_def_var_deflate(iFileID, iCValVar, 0, yDeflate, iDeflate)

    !!! dAlpha/dRho (f)
    iStatus = nf90_def_var(iFileID, 'dAlphadRho', nf90_double, &
                           (/nPsiDim/), idAdRVar)
    iStatus = nf90_def_var_deflate(iFileID, idAdRVar, 0, yDeflate, iDeflate)

    !!! dBeta/dPhi (fzet)
    iStatus = nf90_def_var(iFileID, 'dBetadPhi', nf90_double, &
                           (/nZetaDim/), idBdPVar)
    iStatus = nf90_def_var_deflate(iFileID, idBdPVar, 0, yDeflate, iDeflate)

    !! MISC
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

    !! Plasmasphere Variables
    iStatus = nf90_def_var(iFileID, 'necr', nf90_double,  &
                           (/nRDim,nTDim/), iNECRVar)
    iStatus = nf90_def_var(iFileID, 'tau', nf90_double,  &
                           (/nRDim,nTDim/), itauVar)
 
    ! END DEFINE MODE
    iStatus = nf90_enddef(iFileID)
    call ncdf_check(iStatus, NameSub)


    ! START WRITE MODE
    !! TEMP STUFF
    do iS = 1, nS
        iStatus = nf90_inq_varid(iFileID, 'EnergyGrid_'//species(iS)%s_name, iEGridVar)
        iStatus = nf90_put_var(iFileID, iEGridVar, EKEV(iS,:))
    enddo
    iStatus = nf90_put_var(iFileID, iPaGridVar, PA(:))

    !! FLUXES
    do iS = 1, nS
        iStatus = nf90_inq_varid(iFileID, 'Flux_'//species(iS)%s_name, iFluxVar)
        iStatus = nf90_put_var(iFileID, iFluxVar,  F2(iS,:,:,:,:))
    enddo

    !! PRESSURES
    do iS = 1, nS
        iStatus = nf90_inq_varid(iFileID, 'PParT_'//species(iS)%s_name, iPParTVar)
        iStatus = nf90_put_var(iFileID, iPParTVar, PParT(iS,:,:))
        iStatus = nf90_inq_varid(iFileID, 'PPerT_'//species(iS)%s_name, iPPerTVar)
        iStatus = nf90_put_var(iFileID, iPPerTVar, PPerT(iS,:,:))
    enddo

    !! MAGNETIC FIELD
    iStatus = nf90_put_var(iFileID, iBxVar, bX(:,:,:))
    iStatus = nf90_put_var(iFileID, iByVar, bY(:,:,:))
    iStatus = nf90_put_var(iFileID, iBzVar, bZ(:,:,:))

    !! ELECTRIC FIELD
    iStatus = nf90_put_var(iFileID, iVTOLVar, VTOL(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iVTNVar,  VTN(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iVTVar,   VT(2:nR+1,:))

    !! hI OUTPUTS
    iStatus = nf90_put_var(iFileID, iHVar,    FNHS(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iBHVar,   BOUNHS(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iIVar,    FNIS(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iBIVar,   BOUNIS(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iBNESVar, BNES(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iHDNSVar, HDNS(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iEIRVar,  EIR(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iEIPVar,  EIP(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iDBDTVar, dBdt(2:nR+1,:))
    iStatus = nf90_put_var(iFileID, iDIDTVar, dIdt(2:nR+1,:,:))
    iStatus = nf90_put_var(iFileID, iDIBNVar, dIbndt(2:nR+1,:,:))

    !! GRID OUTPUTS
    iStatus = nf90_put_var(iFileID, iXVar, x(:,:,1:nzeta))
    iStatus = nf90_put_var(iFileID, iYVar, y(:,:,1:nzeta))
    iStatus = nf90_put_var(iFileID, iZVar, z(:,:,1:nzeta))
    iStatus = nf90_put_var(iFileID, iThetaVar, thetaVal(:))
    iStatus = nf90_put_var(iFileID, iRhoVar,   rhoVal(:))
    iStatus = nf90_put_var(iFileID, iZetaVar,  zetaVal(:))

    !! ALPHA/BETA
    iStatus = nf90_put_var(iFileID, iAValVar,  psival(:))
    iStatus = nf90_put_var(iFileID, iBValVar,  alphaval(1:nzeta))
    iStatus = nf90_put_var(iFileID, iCValVar,  chiVal(:))
    iStatus = nf90_put_var(iFileID, idAdRVar,  f(:))
    iStatus = nf90_put_var(iFileID, idBdPVar,  fzet(1:nzeta))

    !! MISC
    iStatus = nf90_put_var(iFileID, iDtVar, DtsNext)
    iStatus = nf90_put_var(iFileID, iTOldVar, TOld)
    iStatus = nf90_put_var(iFileID, icTVar, constTheta)
    iStatus = nf90_put_var(iFileID, icZVar, constZ)
    iStatus = nf90_put_var(iFileID, ipsiInVar, xpsiin)
    iStatus = nf90_put_var(iFileID, ipsiOutVar, xpsiout)
    iStatus = nf90_put_var(iFileID, ikmaxVar, kmax)

    !! Plasmasphere
    iStatus = nf90_put_var(iFileID, iNECRVar, NECR)
    iStatus = nf90_put_var(iFileID, itauVar, tau)

    ! END WRITE MODE
    call ncdf_check(iStatus, NameSub)
    
    ! CLOSE FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_restart

  !==========================================================================
  subroutine read_restart
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRestartIn
    use ModRamGrids,     ONLY: nPa, nT, nR, nS, nE
    use ModRamFunctions, ONLY: RamFileName
    use ModRamTiming,    ONLY: DtsNext, TOld
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP, PParH, PPerH, PParO, PAbn, &
                               PPerO, PParHe, PPerHe, PParE, PPerE, LZ, MU, &
                               NECR, tau, EKEV, species
    use ModScbGrids,     ONLY: npsi, nzeta
    use ModScbParams,    ONLY: constTheta, constZ
    use ModScbVariables, ONLY: x, y, z, alphaVal, psiVal, chiVal, xpsiout, &
                               xpsiin, left, right, kmax, thetaVal, f, fp, &
                               fzet, fzetp, zetaVal, rhoVal, bX, bY, bZ
    !!!! Module Subroutines/Functions
    use ModRamNCDF,      ONLY: ncdf_check, write_ncdf_globatts
    use ModRamFunctions, ONLY: FUNT, FUNI
    !!!! NetCdf Modules
    use netcdf

    use nrtype, ONLY: DP, pi_d, twopi_d
    implicit none
    
    integer :: j, L, iS
    integer :: iFluxVar, iPParTVar, iEGridVar, &
               iPPerTVar, iHVar, iBHVar, iBxVar, iByVar, iBzVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, &
               iTOldVar, iCValVar, icTVar, icZVar, ipsiInVar, ipsiOutVar, &
               ikMaxVar, idAdRVar, idBdPVar, iThetaVar, iRhoVar, iZetaVar, &
               iNECRVar, itauVar

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
    !! MAGNETIC FIELD
    iStatus = nf90_inq_varid(iFileID, 'Bx', iBxVar)
    iStatus = nf90_inq_varid(iFileID, 'By', iByVar)
    iStatus = nf90_inq_varid(iFileID, 'Bz', iBzVar)

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
    iStatus = nf90_inq_varid(iFileID, 'dAlphadRho', idAdRVar)
    iStatus = nf90_inq_varid(iFileID, 'dBetadPhi',  idBdPVar)

    !! MISC
    iStatus = nf90_inq_varid(iFileID, 'DtNext', iDtVar)
    iStatus = nf90_inq_varid(iFileID, 'TOld', iTOldVar)
    iStatus = nf90_inq_varid(iFileID, 'constTheta', icTVar)
    iStatus = nf90_inq_varid(iFileID, 'constZ',     icZVar)
    iStatus = nf90_inq_varid(iFileID, 'xpsiin',     ipsiInVar)
    iStatus = nf90_inq_varid(iFileID, 'xpsiout',    ipsiOutVar)
    iStatus = nf90_inq_varid(iFileID, 'kmax',       ikmaxVar)

    !! Plasmasphere
    iStatus = nf90_inq_varid(iFileID, 'necr', iNECRVar)
    iStatus = nf90_inq_varid(iFileID, 'tau', itauVar)

    ! Check which version we are reading
    iStatus = nf90_inq_varid(iFileID, 'FluxE', iFluxVar)
    if (iStatus == nf90_noerr) then
        iStatus = nf90_inq_varid(iFileID, 'FluxE',  iFluxVar)
        iStatus = nf90_get_var(iFileID, iFluxVar,  F2(1,:,:,:,:))

        iStatus = nf90_inq_varid(iFileID, 'FluxH',  iFluxVar)
        iStatus = nf90_get_var(iFileID, iFluxVar,  F2(2,:,:,:,:))

        iStatus = nf90_inq_varid(iFileID, 'FluxHe',  iFluxVar)
        iStatus = nf90_get_var(iFileID, iFluxVar,  F2(3,:,:,:,:))

        iStatus = nf90_inq_varid(iFileID, 'FluxO',  iFluxVar)
        iStatus = nf90_get_var(iFileID, iFluxVar,  F2(4,:,:,:,:))

        iStatus = nf90_inq_varid(iFileID, 'PParT', iPParTVar)
        iStatus = nf90_get_var(iFileID, iPParTVar, PParT(:,:,:))

        iStatus = nf90_inq_varid(iFileID, 'PPerT', iPPerTVar)
        iStatus = nf90_get_var(iFileID, iPPerTVar, PPerT(:,:,:))
    else 
        do iS = 1, nS
            iStatus = nf90_inq_varid(iFileID, 'EnergyGrid_'//species(iS)%s_name, iEGridVar)
            iStatus = nf90_get_var(iFileID, iEGridVar, EKEV(iS,:))

            iStatus = nf90_inq_varid(iFileID, 'Flux_'//species(iS)%s_name, iFluxVar)
            iStatus = nf90_get_var(iFileID, iFluxVar,  F2(iS,:,:,:,:))

            iStatus = nf90_inq_varid(iFileID, 'PParT_'//species(iS)%s_name, iPParTVar)
            iStatus = nf90_get_var(iFileID, iPParTVar, PParT(iS,:,:))

            iStatus = nf90_inq_varid(iFileID, 'PPerT_'//species(iS)%s_name, iPPerTVar)
            iStatus = nf90_get_var(iFileID, iPPerTVar, PPerT(iS,:,:))
        enddo
    endif

    !! MAGNETIC FIELD
    iStatus = nf90_get_var(iFileID, iBxVar, bX(:,:,:))
    iStatus = nf90_get_var(iFileID, iByVar, bY(:,:,:))
    iStatus = nf90_get_var(iFileID, iBzVar, bZ(:,:,:))

    !! ELECTRIC FIELD
    iStatus = nf90_get_var(iFileID, iVTOLVar, VTOL(2:nR+1,:))
    iStatus = nf90_get_var(iFileID, iVTNVar,  VTN(2:nR+1,:))
    iStatus = nf90_get_var(iFileID, iVTVar,   VT(2:nR+1,:))
    VTOL(1,:) = VTOL(2,:)
    VTN(1,:) = VTN(2,:)
    VT(1,:) = VT(2,:)

    !! hI OUTPUTS
    iStatus = nf90_get_var(iFileID, iHVar,    FNHS(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iBHVar,   BOUNHS(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iIVar,    FNIS(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iBIVar,   BOUNIS(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iBNESVar, BNES(2:nR+1,:))
    iStatus = nf90_get_var(iFileID, iHDNSVar, HDNS(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iEIRVar,  EIR(2:nR+1,:))
    iStatus = nf90_get_Var(iFileID, iEIPVar,  EIP(2:nR+1,:))
    iStatus = nf90_get_var(iFileID, iDBDTVar, dBdt(2:nR+1,:))
    iStatus = nf90_get_var(iFileID, iDIDTVar, dIdt(2:nR+1,:,:))
    iStatus = nf90_get_var(iFileID, iDIBNVar, dIbndt(2:nR+1,:,:))
    DO J=1,NT ! use dipole B at I=1
       BNES(1,J)=0.32/LZ(1)**3/1.e4
       dBdt(1,J) = 0._dp
       EIR(1,J) = 0._dp
       EIP(1,J) = 0._dp
       DO L=1,NPA
          FNHS(1,J,L) = FNHS(2,J,L)
          FNIS(1,J,L) = FNIS(2,J,L)
          BOUNHS(1,J,L)=BOUNHS(2,J,L)
          BOUNIS(1,J,L)=BOUNIS(2,J,L)
          HDNS(1,J,L)=HDNS(2,J,L)
          dIdt(1,J,L)=0._dp
          dIbndt(1,J,L)=0._dp
       ENDDO
    ENDDO

    !! GRID OUTPUTS
    iStatus = nf90_get_var(iFileID, iXVar, x(:,:,1:nzeta))
    iStatus = nf90_get_var(iFileID, iYVar, y(:,:,1:nzeta))
    iStatus = nf90_get_var(iFileID, iZVar, z(:,:,1:nzeta))
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)

    iStatus = nf90_get_var(iFileID, iThetaVar, thetaVal(:))
    iStatus = nf90_get_var(iFileID, iRhoVar, rhoVal(:))
    iStatus = nf90_get_var(iFileID, iZetaVar, zetaVal(:))

    !! ALPHA/BETA
    iStatus = nf90_get_var(iFileID, iAValVar,  psival(:))
    iStatus = nf90_get_var(iFileID, iBValVar,  alphaval(1:nzeta))
    alphaVal(nzeta+1) = alphaVal(2) + twopi_d
    iStatus = nf90_get_var(iFileID, iCValVar,  chiVal(:))
    iStatus = nf90_get_var(iFileID, idAdRVar,  f(:))
    iStatus = nf90_get_var(iFileID, idBdPVar,  fzet(1:nzeta))
    fzet(nZeta+1) = fzet(2)
    fp = 0._dp
    fzetp = 0._dp

    !! MISC
    iStatus = nf90_get_var(iFileID, iDtVar, DTsNext)
    iStatus = nf90_get_var(iFileID, iTOldVar, TOld)
    iStatus = nf90_get_var(iFileID, icTVar, constTheta)
    iStatus = nf90_get_var(iFileID, icZVar, constZ)
    iStatus = nf90_get_var(iFileID, ipsiInVar, xpsiin)
    iStatus = nf90_get_var(iFileID, ipsiOutVar, xpsiout)
    iStatus = nf90_get_var(iFileID, ikMaxVar, kMax)

    !! Plasmasphere
    iStatus = nf90_get_var(iFileID, iNECRVar, NECR)
    iStatus = nf90_get_var(iFileID, itauVar, tau)

    ! CLOSE RESTART FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine read_restart

!==============================================================================
END MODULE ModRamRestart
