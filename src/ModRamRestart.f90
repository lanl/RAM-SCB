!============================================================================
module ModRamRestart
!   A set of variables and tools for general I/O handling, file name
!   creation, etc.
!   Copyright (c) 2016, Los Alamos National Security, LLC
!   All rights reserved.
!******************************************************************************
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRestartOut, PathRestartIn, niter
    use ModRamFunctions, ONLY: RamFileName
    use ModRamTiming,    ONLY: TimeRamElapsed, TimeRamStart, TimeRamNow, DtsNext
    use ModRamGrids,     ONLY: NR, NT, NE, NPA
    use ModRamVariables, ONLY: F2, PParT, PPerT, FNHS, FNIS, BOUNHS, BOUNIS, &
                               BNES, HDNS, EIR, EIP, dBdt, dIdt, dIbndt, VTN, &
                               VTOL, VT, EIR, EIP, FGEOS, PParH, PPerH, PParO, &
                               PPerO, PParHe, PPerHe, PParE, PPerE
    use ModRamScb,       ONLY: indexPA, FLUX3DEQ
    use ModScbGrids,     ONLY: nthe, npsi, nzeta, nzetap
    use ModScbVariables, ONLY: x, y, z, bX, bY, bZ, bf, alfa, psi, alphaVal, psiVal
    !!!! Module Subroutines/Functions
    use ModRamNCDF, ONLY: ncdf_check, write_ncdf_globatts
    !!!! Share Modules
    use ModIOUnit, ONLY: UNITTMP_
    !!!! NetCdf Modules
    use netcdf

  implicit none
  save
  
  contains
  !==========================================================================
  subroutine write_restart

    implicit none
    
    integer :: stat
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iPParTVar, &
               iPPerTVar, iBxVar, iByVar, iBzVar, iBfVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iGEOVar, iPaIndexVar, iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAlphaVar, iBetaVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus
    integer :: nRDim, nTDim, nEDim, nPaDim, nSDim, nThetaDim, nPsiDim, &
               nZetaDim, nRPDim, nZetaPDim, iFlux3DVar

    character(len=2), dimension(4):: NameSpecies = (/'e_','h_','he','o_'/)
    character(len=200)            :: NameFile,CWD

    character(len=*), parameter :: NameSub='write_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    if(DoTest)write(*,'(a,f11.2)') 'RAM-SCB: Writing restarts at t=',&
         TimeRamElapsed

    ! Write ascii portion of restart.
    NameFile=RamFileName(PathRestartOut//'/restart_info','txt',TimeRamNow)
    !NameFile=PathRestartOut//'/restart_info.txt'
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
    NameFile = RamFileName(PathRestartOut//'/restart','nc',TimeRamNow)
    !NameFile = PathRestartOut//'/restart.nc'
    iStatus = nf90_create(trim(NameFile), nf90_clobber, iFileID)
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
    iStatus = nf90_def_var(iFileID, 'FluxE', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxEVar)
    iStatus = nf90_put_att(iFileID, iFluxEVar, 'title', &
                           'This is an example title')

    iStatus = nf90_def_var(iFileID, 'FluxH', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxHVar)
    iStatus = nf90_def_var(iFileID, 'FluxO', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxOVar)
    iStatus = nf90_def_var(iFileID, 'FluxHe', nf90_double, &
                           (/nRdim,nTDim,nEDim,nPaDim/), iFluxHeVar)
    iStatus = nf90_def_var(iFileID, 'FGEOS', nf90_double, &
                           (/nSdim,nTDim,nEDim,nPadim/), iGEOVar) 
    iStatus = nf90_def_var(iFileID, 'Flux3D', nf90_double, &
                           (/nSdim,nPsiDim,nZetaDim,nEDim,nPADim/), iFlux3DVar)

    !! PRESSURES
    iStatus = nf90_def_var(iFileID, 'PParT', nf90_double, &
                           (/nSdim,nRDim,nTDim/), iPParTVar)
    iStatus = nf90_def_var(iFileID, 'PPerT', nf90_double, &
                           (/nSdim,nRDim,nTDim/), iPPerTVar)

    !! MAGNETIC FIELD
    iStatus = nf90_def_var(iFileID, 'Bx', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iBxVar)
    iStatus = nf90_def_var(iFileID, 'By', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iByVar)
    iStatus = nf90_def_var(iFileID, 'Bz', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaDim/), iBzVar)
    iStatus = nf90_def_var(iFileID, 'B', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iBfVar)

    !! ELECTRIC FIELD
    iStatus = nf90_def_var(iFileID, 'VTOL', nf90_double, &
                           (/nRPDim,nTDim/), iVTOLVar)
    iStatus = nf90_def_var(iFileID, 'VTN', nf90_double, &
                           (/nRPDim,nTDim/), iVTNVar)
    iStatus = nf90_def_var(iFileID, 'VT', nf90_double, &
                           (/nRPDim,nTDim/), iVTVar)

    !! hI OUTPUTS
    iStatus = nf90_def_var(iFileID, 'FNHS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iHVar)
    iStatus = nf90_def_var(iFileID, 'BOUNHS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iBHVar)
    iStatus = nf90_def_var(iFileID, 'FNIS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iIVar)
    iStatus = nf90_def_Var(iFileID, 'BOUNIS', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iBIVar)
    iStatus = nf90_def_var(iFileID, 'BNES', nf90_double, &
                           (/nrPDim,nTDim/), iBNESVar)
    iStatus = nf90_def_var(iFileID, 'HDNS', nf90_double, &
                           (/nRPDim,nTDim,nPaDim/), iHDNSVar)
    iStatus = nf90_def_var(iFileID, 'EIR', nf90_double, &
                           (/nrPDim,nTDim/), iEIRVar)
    iStatus = nf90_def_var(iFileID, 'EIP', nf90_double, &
                           (/nrPDim,nTDim/), iEIPVar)
    iStatus = nf90_def_var(iFileID, 'dBdt', nf90_double, &
                           (/nrPDim,nTDim/), iDBDTVar)
    iStatus = nf90_def_var(iFileID, 'dIdt', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iDIDTVar)
    iStatus = nf90_def_var(iFileID, 'dIbndt', nf90_double, &
                           (/nRPDim,nTDim,NPaDim/), iDIBNVar)

    !! GRID OUTPUTS
    iStatus = nf90_def_var(iFileID, 'x', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iXVar)
    iStatus = nf90_def_var(iFileID, 'y', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iYVar)
    iStatus = nf90_def_var(iFileID, 'z', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iZVar)

    !! ALPHA/BETA
    iStatus = nf90_def_var(iFileID, 'alphaVal', nf90_double, &
                           (/nZetaPDim/), iAValVar)
    iStatus = nf90_def_var(iFileID, 'betaVal', nf90_double, &
                           (/nPsiDim/), iBValVar)
    iStatus = nf90_def_var(iFileID, 'alpha', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iAlphaVar)
    iStatus = nf90_def_var(iFileID, 'beta', nf90_double, &
                           (/nThetaDim,nPsiDim,nZetaPDim/), iBetaVar)

    !! MISC
    iStatus = nf90_def_var(iFileID, 'IndexPA', nf90_int, &
                           (/nThetaDim,nPsiDim,nZetaDim,nPaDim/), iPaIndexVar)
    iStatus = nf90_def_var(iFileID, 'DtNext', nf90_double, iDtVar)

    ! END DEFINE MODE
    iStatus = nf90_enddef(iFileID)
    call ncdf_check(iStatus, NameSub)


    ! START WRITE MODE
    !! FLUXES
    iStatus = nf90_put_var(iFileID, iFluxEVar,  F2(1,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHVar,  F2(2,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxHeVar, F2(3,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFluxOVar,  F2(4,:,:,:,:))
    iStatus = nf90_put_var(iFileID, iGEOVar,    FGEOS(:,:,:,:))
    iStatus = nf90_put_var(iFileID, iFlux3DVar, FLUX3DEQ(:,:,:,:,:))

    !! PRESSURES
    iStatus = nf90_put_var(iFileID, iPParTVar, PParT(:,:,:))
    iStatus = nf90_put_var(iFileID, iPPerTVar, PPerT(:,:,:))

    !! MAGNETIC FIELD
    iStatus = nf90_put_var(iFileID, iBxVar, bX(:,:,:))
    iStatus = nf90_put_var(iFileID, iByVar, bY(:,:,:))
    iStatus = nf90_put_var(iFileID, iBzVar, bZ(:,:,:))
    iStatus = nf90_put_var(iFileID, iBfVar, bf(:,:,:))

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

    !! ALPHA/BETA
    iStatus = nf90_put_var(iFileID, iAValVar, alphaval(:))
    iStatus = nf90_put_var(iFileID, iBValVar, psival(:))
    iStatus = nf90_put_var(iFileID, iAlphaVar, alfa(:,:,:))
    iStatus = nf90_put_var(iFileID, iBetaVar, psi(:,:,:))

    !! MISC
    iStatus = nf90_put_var(iFileID, iPaIndexVar, indexPA(:,:,:,:))
    iStatus = nf90_put_var(iFileID, iDtVar, DtsNext)

    ! END WRITE MODE
    call ncdf_check(iStatus, NameSub)
    
    ! CLOSE FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_restart

  !==========================================================================
  subroutine read_restart

    implicit none
    
    integer :: nrIn, ntIn, neIn, npaIn
    integer :: iFluxEVar, iFluxHVar, iFluxHeVar, iFluxOVar, iPParTVar, &
               iPPerTVar, iBxVar, iByVar, iBzVar, iBfVar, iHVar, iBHVar, &
               iIVar, iBIVar, iBNESVar, iHDNSVar, iEIRVar, iEIPVar, &
               iGEOVar, iPaIndexVar, iDtVar, iXVar, iYVar, iZVar, &
               iVTNVar, iAlphaVar, iBetaVar, iAValVar, iBValVar, iVTOLVar, &
               iVTVar, iDBDTVar, iDIDTVar, iDIBNVar, iFileID, iStatus, iFlux3DVar

    character(len=2), dimension(4) :: NameSpecies = (/ 'e_','h_','he','o_' /)
    character(len=100)             :: NameFile, StringLine

    character(len=*), parameter :: NameSub='read_restart'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    call write_prefix
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
    iStatus = nf90_inq_varid(iFileID, 'FGEOS',  iGEOVar)
    iStatus = nf90_inq_varid(iFileID, 'Flux3D', iFlux3DVar)

    !! PRESSURES
    iStatus = nf90_inq_varid(iFileID, 'PParT', iPParTVar)
    iStatus = nf90_inq_varid(iFileID, 'PPerT', iPPerTVar)

    !! MAGNETIC FIELD
    iStatus = nf90_inq_varid(iFileID, 'Bx', iBxVar)
    iStatus = nf90_inq_varid(iFileID, 'By', iByVar)
    iStatus = nf90_inq_varid(iFileID, 'Bz', iBzVar)
    iStatus = nf90_inq_varid(iFileID, 'B',  iBfVar)

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

    !! ALPHA/BETA
    iStatus = nf90_inq_varid(iFileID, 'alphaVal', iAValVar)
    iStatus = nf90_inq_varid(iFileID, 'betaVal', iBValVar)
    iStatus = nf90_inq_varid(iFileID, 'alpha', iAlphaVar)
    iStatus = nf90_inq_varid(iFileID, 'beta', iBetaVar)

    !! MISC
    iStatus = nf90_inq_varid(iFileID, 'IndexPA', iPaIndexVar)
    iStatus = nf90_inq_varid(iFileID, 'DtNext', iDtVar)

    ! READ DATA
    !! FLUXES
    iStatus = nf90_get_var(iFileID, iFluxEVar,  F2(1,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHVar,  F2(2,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxHeVar, F2(3,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFluxOVar,  F2(4,:,:,:,:))
    iStatus = nf90_get_var(iFileID, iGEOVar,    FGEOS(:,:,:,:))
    iStatus = nf90_get_var(iFileID, iFlux3DVar, FLUX3DEQ(:,:,:,:,:))

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
    iStatus = nf90_get_var(iFileID, iBxVar, bX(:,:,:))
    iStatus = nf90_get_var(iFileID, iByVar, bY(:,:,:))
    iStatus = nf90_get_var(iFileID, iBzVar, bZ(:,:,:))
    iStatus = nf90_get_Var(iFileID, iBfVar, bf(:,:,:))

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

    !! ALPHA/BETA
    iStatus = nf90_get_var(iFileID, iAValVar, alphaval(:))
    iStatus = nf90_get_var(iFileID, iBValVar, psival(:))
    iStatus = nf90_get_var(iFileID, iAlphaVar, alfa(:,:,:))
    iStatus = nf90_get_var(iFileID, iBetaVar, psi(:,:,:))

    !! MISC
    iStatus = nf90_get_var(iFileID, iPaIndexVar, indexPA(:,:,:,:))
    iStatus = nf90_get_var(iFileID, iDtVar, DTsNext)

    ! CLOSE RESTART FILE
    iStatus = nf90_close(iFileID)
    call ncdf_check(iStatus, NameSub)

  end subroutine read_restart

!==============================================================================
END MODULE ModRamRestart
