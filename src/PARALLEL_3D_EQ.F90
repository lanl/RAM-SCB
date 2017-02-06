SUBROUTINE PARALLEL_3D_EQ(rhour_local,ST3, iPressure_local, iComp_domain_local, electric_local, event_local)

!******************************************************************************
! For compiling see Makefile
!
! Description: 
!  Calculates the inverse 3-D plasma equilibrium equation in Euler potential
!  coordinates in the closed-field region of the magnetosphere
!
! Numerical scheme for elliptic PDE solution direct inversion:
! (Block Tri-diagonal method + Gauss-Jordan) or SOR
! (SOR faster for most practical apps)
! 
! Method: Iterative flux coordinate scheme 
! [Cheng, 1995; Zaharia et al., 2004]; solves iteratively sets of quasi-2D eqs.
!  
! For each 2-D PDE: Numerical scheme for elliptic PDE solution direct inversion 
! (Block Tri-diagonal method + Gauss-Jordan) or SOR
! (SOR faster for most practical apps)
! 
! Input files: 
! Magnetic field: t****.cdf (usually from tracing Tsyganenko model field, 
! see separate tracing routine)
! Pressure: pressure***.in  ! Pressure input through subroutine pressure.f90
!
! Output files: 
! pressure***.cdf Pressure, anisotropy and other output 
! mag_field***.cdf Magnetic field, coordinates and other output
! hI****.dat Mag-field line integrals/Equatorial B for use with coupling 
! with RAM Output (hI*dat files) through subroutine hRAM.f90
!
! Code Description: 
! Language: Fortran 90 (some F95 constructs as well)
! Author: S. Zaharia, 2003-2009
!
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!******************************************************************************


  USE nrtype
  USE Module1
  USE mpi
  USE ModRamMain, ONLY: PathRamOut, PathScbOut, PathScbIn, Kp, TimeRamNow
  USE ModIoUnit,  ONLY: io_unit_new
  USE ModRamMpi,  ONLY: iProc, nProc, iComm
  Use ModRamIO,   ONLY: RamFileName, UseNewFmt

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rhour_local
  INTEGER, INTENT(IN) :: iPressure_local, iComp_domain_local
  CHARACTER*4, INTENT(IN) :: electric_local
  CHARACTER*4, INTENT(IN) :: ST3  
  INTEGER :: inBlank, iostat, iStart, iStop, iRate, ierr, length 
  REAL(DP) :: time, xTry
  CHARACTER(LEN = MPI_MAX_PROCESSOR_NAME) :: procName
  CHARACTER(LEN = 6) :: event_local

!C	print*, 'P3DEQ: iPressure, iDomain = ', iPressure_local, iComp_domain_local	

  !C CALL MPI_INIT(ierr)  ! SZ, put it in Main
  start_time = MPI_WTIME() ! Start timer
  rank = iProc
  !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  ! Get the total number of processors:
  numProc = nProc
  !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProc, ierr)
  CALL MPI_GET_PROCESSOR_NAME(procName, length, ierr)
  IF (rank == 0) PRINT*, 'Number of procs. = ', numProc; CALL FLUSH(6)

  rHour = rhour_local ! Computation only for one hour, if this is a subroutine for RAM coupling
  event = event_local 

   iDumpRAMFlux = 0 ! Writes RAM flux mapped along 3D field lines in NetCDF format

  isotropy = 0 ! Anisotropic pressure case 

  iInducedE = 0 ! No calculation/writing of induced E-fields	
  iConvE = 1 ! No calculation/writing of convective E-fields (if zero)

  ! For runs with RAM pressures, only the following choices should be used

  iReduceAnisotropy = 0 ! No change in anisotropy 
  ! iReduceAnisotropy = 1 ! Change anisotropy to marginally mirror-stable

  iOuterMethod = 1 ! DO NOT CHANGE FOR RAM-SCB ! Picard iteration  ! Newton has not been benchmarked for anisotropic pressure
  !C iOuterMethod = 2 ! Newton iteration  

  ! Given as arguments now
  iPressureChoice = iPressure_local
  ! iPressureChoice = 5 ! RAM pressures (domain inside 6.6 RE)
  ! iPressureChoice = 6 ! RAM pressures, with extension up to 10 RE via Roeder model
  ! iPressureChoice = 11 ! SWMF driving, RAM pressures inside 6.6 RE and SWMF pressures on boundary/outside

  iInterpMethod = 1 ! DO NOT CHANGE FOR RAM-SCB ! PSPLINE interpolation
 
  !method = 2 ! DO NOT CHANGE FOR RAM-SCB 
  ! SOR computation of PDE; faster than direct method for N larger than about 41 ! method = 3 means no calculation
  ! method NOW SET IN INIT_RAMSCB.F90; DEFAULTS TO 2.

  ! iWantAlphaExtrapolation = 0 ! Do not extrapolate alpha (beta) on the first/last flux surface
  iWantAlphaExtrapolation = 1 ! Extrapolate alpha (beta) on the first/last flux surface

  numit = 200 ! DO NOT CHANGE FOR RAM-SCB Max. # of iterations; the actual number is much lower (code exits on convergence) and depends on convergence speed 

  iConvergenceMethod = 2 ! 1 (global force balance decrease) or 2 (residual between iterations decrease)  

  ! decreaseConv - the amount the global force imbalance decreases throughout iterations (vs. the initial one)
  decreaseConv = 0.5_dp
  ! decreaseConv = 0.33 ! 3X reduction, should be enough for RAM-SCB (very hard to get better w/ coarse grid resolution)
  ! decreaseConv = 0.2   ! 5X reduction, if grid resolution is better
  ! decreaseConv = 0.1   ! 10X reduction, and more needed if state has to be a near-perfect equilibrium for theoretical calculations

  blendInitial = 0.25_dp ! "Blends" solution at iterations (n+1) and n

  decreaseConv = 0.05 ! 20X reduction
  
  decreaseConvAlpha = decreaseConvAlphaMin + (decreaseConvAlphaMax - decreaseConvAlphaMin)*(MIN(Kp,6._dp))**2/36.
  decreaseConvPsi = decreaseConvPsiMin + (decreaseConvPsiMax - decreaseConvPsiMin)*(MIN(Kp,6._dp))**2/36.
  
  print*, 'P3DEQ: Kp, decreaseConv = ', Kp, decreaseConvPsi

  blendAlpha = blendAlphaInit
  blendPsi = blendPsiInit

  thresh = 1.1
  damp = 0.5
  relax = 1.5
  nrelax = 5 ! Try to increase blending factor after nrelax iterations

  iSm = 0 !Even smoother grid respacing, for difficult equilibria
  iSm2 = 1 ! Savitzky-Golay smoothing in pressure.f90 (for RAM-SCB mostly)
  
  nimax = 5000  ! DO NOT CHANGE FOR RAM-SCB ! # of inner iterations in SOR technique; usually 2000 is more than enough for SOR convergence

  iAMR = 1 ! Mesh refinement in magnetic flux, so that one has equidistant magnetic flux surfaces; improves convergence a lot

  iAzimOffset = 2 ! Equidistance sought for most problematic local time
  !C iAzimOffset = 1 ! Equidistance maintained at midnight

  isSORDetailNeeded = 0 ! No details for the inner SOR iterations 
  ! isSORDetailNeeded = 1 ! Details about inner SOR iterations

  ! isEnergDetailNeeded = 0
  isEnergDetailNeeded = 1 ! Dst computation (DPS formula with thermal energy inside domain)

  isFBDetailNeeded = 0 ! Computes global force imbalance

  ! iLossCone = 2 ! More realistic, empty loss cone for RAM computations (M. Liemohn's formalism, Liemohn, 2004); 
  ! can have some errors with very deformed B-field
  iLossCone = 1 ! Filled loss cone
  iOutput = 1 ! Output for Vania's RAM - hI.dat files are written to disk

  iTilt = 0     ! Zero dipole tilt
  ! iTilt = 1   ! Case with non-zero magnetic tilt

  iHourChoice = 0 ! More than 1 equilibrium calculation (for RAM-SCB)

  ! Given as argument now from Main.f
  iCompDomain = iComp_domain_local
  ! iCompDomain = 2 ! TS04 boundary conditions 
  ! iCompDomain = 20 ! SWMF B-field boundary conditions
  ! iCompDomain = 3 ! T89c-based boundary conditions

  IF (electric_local == 'IESC') iElectric = 1 ! Mapping of iono. potentials along SC B-field lines
  IF (electric_local == 'WESC') iElectric = 2 ! Weimer electric iono. potentials mapped along SC B-field  
  IF (electric_local == 'W5SC') iElectric = 3 ! Weimer electric 2005 iono. potentials along SC B-field

  prefixIn     = trim(PathScbIn)  // '/'
  prefixOut    = trim(PathScbOut) // '/'
  prefixRAMOut = trim(PathRamOut) // '/'

  READ(ST3, '(I4)') iST3 ! Index of the equilibrium snapshot

  Only_if_RAM_pressures:      IF (iPressureChoice == 5 .OR. iPressureChoice==6 .OR. iPressureChoice==11) THEN
     CALL MPI_BARRIER(iComm, ierr)
     iCallRoutine = 0 
     WRITE (HourDecimal,'(F5.2)') rHour    ! Internal write to transform rHour into character HourDecimal 

     iDay = FLOOR(rHour/24.)
     !     print*, 'rHour, iDay = ', rHour, iDay
     WRITE(Day, '(I2)') iDay
     IF (iDay < 10) Day(1:1) = '0'
     iHour = MOD(rHour, 24.)
     iHourAbs = FLOOR(rHour) ! Without taking account of the day
     iMin = INT(ABS(rHour-iHourAbs)*60) ! Minutes	 
     print*, 'iDay, iHour, iMin = ', iDay, iHour, iMin	 
     WRITE(HourInteger, '(I2.2)') iHour	       
     IF (iHour < 23) THEN
        WRITE(DayP1,'(I2)') iDay
        WRITE(HourIntegerP1, '(I2)') iHour+1
     ELSE
        WRITE(DayP1,'(I2)') iDay+1
        WRITE(HourIntegerP1, '(I2)') 0
     END IF
     IF (iDay < 10 .AND. iHour<23) DayP1(1:1) = '0'
     IF (iDay < 9 .AND. iHour==23) DayP1(1:1) = '0'
     IF (iHour < 9 .OR. iHour==23) HourIntegerP1(1:1) = '0'

     WRITE(HourDecShort, '(F5.2)') rHour-iDay*24 ! Hour inside the day)
     IF (rHour-iDay*24 < 10.) HourDecShort(1:1) = '0'

     iSWMF = iST3

     ! Define variable decreaseConv
     decreaseConv = 0.5 ! + (1.-0.5)*exp(-0.1*real(iST3))

     PRINT*, 'P3DEQ: iST3, decreaseConv = ', iST3, real(decreaseConv,SP)
     print*, ''
     !C     IF (rank == 0) PRINT*, 'Day, Hour short, Hour = ', Day, HourDecShort, HourDecimal; CALL FLUSH(6)
     if(UseNewFmt)then
	fileNamePressure=trim(prefixRAMOut)//RamFileName('pressure','in',TimeRamNow)
     else
	fileNamePressure = ADJUSTL(TRIM(prefixRAMOut)//'pressure_'//ST3//'.in') ! all pressure files in one directory
     end if
     inBlank = INDEX(fileNamePressure, ' ') - 1
     IF (rank == 0) THEN
        PRINT*, 'Pressure file is ', fileNamePressure(1:inBlank); CALL FLUSH(6)
     END IF

     iMin = MOD(real(iST3*5+2),60.)
     WRITE(MinChar,'(I2.2)') iMin
     !C     print*, 'MinChar = ', MinChar

  END IF Only_if_RAM_pressures


  CALL Computational_domain(iCompDomain)
  CALL input
  itout = 0

  ! Initial guess for alpha(zeta) and psi(rho) 
  CALL alfges
  CALL psiges
  CALL MPI_BARRIER(iComm, ierr)

  !cl      CALL computeEquilibrium  ! Main computation routine
  CALL computeEquilibrium(ST3)  ! Main computation routine

  CALL MPI_BARRIER(iComm, ierr)

  end_time = MPI_WTIME( )          ! stop timer

  !C  IF(rank == 0) THEN
  !C     OPEN(35, file='time_output', form = 'formatted', status = 'unknown', action = 'write')
  !C     WRITE(35,*) '#######################################'
  !C     WRITE(35,*) 'Total cpu time (s) =',end_time - start_time,' x', numProc
  !C     CLOSE(35)
  !C     PRINT*, ' '
  !C     PRINT *, '#######################################'
  !C     PRINT *, 'Total cpu time (s) =',end_time - start_time,' x', numProc
  !C     PRINT*, ' '
  !C     CALL FLUSH(6)
  !C  ENDIF

!C  IF (rank == 0) PRINT*, 'Equilibrium computation finished'; CALL flush(6)
  CALL MPI_BARRIER(iComm, ierr)

  RETURN

END SUBROUTINE PARALLEL_3D_EQ

