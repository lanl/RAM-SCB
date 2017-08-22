!==============================================================================
module ModRamMain
!    This module replaces the common blocks in RAM.
!    DTW, EDIT: 2009-04-21: updated for latest version of RAM.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModTimeConvert, ONLY: TimeType
  use nrtype,         ONLY: DP, SP
  use ModRamGrids,    ONLY: NS, NR, NT, NE, NPA

  implicit none
  save

! RAM version.
  real, parameter :: CodeVersion = 3.0

! Kind for double precision, set here for portability.
  integer, parameter :: Real8_ = DP
  integer, parameter :: Real4_ = SP

!!!!! PATHS
  ! File paths for restarts
  character(len=*), parameter :: PathRestartIn = 'IM/restartIN'
  character(len=*), parameter :: PathRestartOut= 'IM/restartOUT'

  ! File paths to fix hardcoding of path issues.
  character(len=*), parameter :: PathRamOut = 'IM/output/ram/'
  character(len=*), parameter :: PathScbOut = 'IM/output/scb/'
  character(len=*), parameter :: PathSwmfOut= 'IM/output_swmf/'
  character(len=*), parameter :: PathRamIn  = 'IM/input_ram/'
  character(len=*), parameter :: PathScbIn  = 'IM/input_scb/'
!!!!!

!!!!! ITERATIONS
  integer :: iCal     ! Tracks iteration number, kind of.
  integer :: nIter=0  ! Tracks iteration, including restart.
!!!!!

!!!!! MAIN RAM VARIABLES (Pressures, Fluxes, and hI variables)
!  real(kind=Real8_), dimension(NS,NR,NT,NE,NPA) :: F2, FLUX
!  real(kind=Real8_), dimension(NR,NT) :: PPerH, PParH, PPerO, PParO, PPerHe, &
!                                         PParHe, PPerE, PParE, PAllSum, PParSum
!  real(kind=Real8_), dimension(NS,NR,NT) :: PPerT, PParT
!  real(kind=Real8_), dimension(NR+1,NT,NPA) :: dIdt, dIbndt, HDNS, FNHS, FNIS, &
!                                               BOUNHS, BOUNIS
!  real(kind=Real8_), dimension(NR+1,NT) :: BNES, dBdt
!!!!!

!\
! PARAMS Block.
!/
integer :: S

!character(len=3) :: ST0, ST1
!character(len=2) :: ST2
!character(len=19):: ST7

end Module ModRamMain
!==============================================================================
