!==============================================================================
module ModPlane
!    This module replaces the variable declarations in plane.h and plane_scb.f
!    CAJ, EDIT: 2012-12-03: original version
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================

  use ModRamPl_Ne, ONLY: NL, NLT, NR, PI, RE, RE_cm, RE_km, KB, NT, Real8_, Real4_


  implicit none
  save

  !\
  ! Grid size parameters
  !/
  
  integer, parameter :: CARTESIAN = 1
  
  !\
  ! Grid parameters
  !/
  
  real(kind=Real8_), parameter :: &
       DL=8.5/(NL-1), &
       DLT=24./NLT, &
       DTR=PI/180., &
       RTD=180./PI
  
  !\
  ! DateTime structure
  !/
  
  common /DateTime/ year, month, day, sec
  
  integer :: year, month, day
  real(kind=Real4_) :: sec
       
  !\
  ! Array variables that are better in modules.
  !/
  real(kind=Real8_) :: nsw(0:NLT,NL)
  
  real(kind=Real4_) :: L(0:NL+1)
  real(kind=Real4_) :: mlt(0:NLT)
  real(kind=Real4_) :: tau0(NL,0:NLT)
  real(kind=Real4_) :: n0(0:NL+1)
  real(kind=Real4_) :: uL(0:NL+1,0:NLT)
  real(kind=Real4_) :: ut(0:NLT,0:NL+1)
  real(kind=Real4_) :: vol(NL)

end Module ModPlane
!==============================================================================
