!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbMain

  USE nrtype,     ONLY: pi_d, SP, DP ! Definition of single/double precision type
  USE ModRamMain, ONLY: PathRamOut, PathScbOut, PathScbIn

  IMPLICIT NONE
  SAVE

  REAL(DP), PARAMETER :: mu0 = 4._dp*pi_d*1.E-7_dp, &
                         BEarth = 0.31_dp*1.E-4_dp, &
                         REarth = 6.4_dp*1.E6_dp 

  character(len=200) :: prefixIn     = trim(PathScbIn)  // '/'
  character(len=200) :: prefixOut    = trim(PathScbOut) // '/'
  character(len=200) :: prefixRAMOut = trim(PathRamOut) // '/'

!!! Intilization
  integer  :: iCountPressureCall   = 0,     &
              iElectric            = 0
              
  REAL(SP) :: rHour


  CHARACTER(LEN=2)   :: Day, DayP1, HourInteger, HourIntegerP1, MinChar
  CHARACTER(Len=5)   :: HourDecimal, HourDecShort
  CHARACTER(LEN=4)   :: suffixFilePress
  CHARACTER(LEN=500) :: fileNamePressure
  CHARACTER(LEN=200) :: prefixPres
  CHARACTER(LEN=2)   :: Hour

  integer  :: tsygcorrect = 0

  integer  :: iInducedE = 0 ! No calculation/writing of induced E-fields

  integer  :: iConvE = 1 ! No calculation/writing of convective E-fields (if zero)

  ! For runs with RAM pressures, only the following choices should be used

  integer  :: iInterpMethod = 1 ! DO NOT CHANGE FOR RAM-SCB ! PSPLINE interpolation

  integer  :: numit = 200 ! DO NOT CHANGE FOR RAM-SCB Max. # of iterations; the actual number is much lower (code exits on convergence) and depends on convergence speed 

  integer  :: iConvergenceMethod = 2 ! 1 (global force balance decrease) or 2 (residual between iterations decrease)  

  ! decreaseConv - the amount the global force imbalance decreases throughout iterations (vs. the initial one)
  real(dp) :: decreaseConv = 0.5 ! 3X reduction, should be enough for RAM-SCB (very hard to get better w/ coarse grid resolution)
 !real(dp) :: decreaseConv = 0.2 ! 5X reduction, if grid resolution is better
 !real(dp) :: decreaseConv = 0.1 ! 10X reduction, and more needed if state has to be a near-perfect equilibrium for theoretical calculations

  real(dp) :: blendInitial = 0.25_dp ! "Blends" solution at iterations (n+1) and n

  real(dp) :: thresh = 1.1

  real(dp) :: damp = 0.90

  real(dp) :: relax = 1.1

  integer  :: nrelax = 100 ! Try to increase blending factor after nrelax iterations

  integer  :: iSm  = 0 ! Even smoother grid respacing, for difficult equilibria

  integer  :: nimax = 5001  ! DO NOT CHANGE FOR RAM-SCB ! # of inner iterations in SOR technique; usually 2000 is more than enough for SOR convergence

END MODULE ModScbMain
