!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamConst

  use ModRamMain,  ONLY: Real8_
  use ModNumConst, ONLY: cPi

  implicit none

  real(kind=Real8_) :: ME    = 7.9E15,    &   ! Magnetic moment of the earth [T*m3]
                       RE    = 6.371E6,   &   ! Earth's radius [m]
                       RE_cm = 6.371E8,   &   ! Earth's radius [cm]
                       RE_km = 6.371E3,   &   ! Earth's radius [km]
                       HMIN  = 2E5,       &   ! Altitude of dense atmosphere [m]
                       MP    = 1.673E-27, &   ! Mass of H+ [kg]
                       Q     = 1.602E-19, &   ! Elementary charge [C]
                       CS    = 2.998E8,   &   ! Speed of light [m/s]
                       PI    = cPi,       &   ! pi used to be 3.141592654
                       KB    = 1.38065E-16, &  ! Boltzmann constant [cm^2 g s^-2 K^-1]
                       b0dip = 30570.0        ! Doesn't really matter what this value is

  integer, parameter :: monthday(12) = (/0,30,59,90,120,151,181,212,243,273,304,334/)
  integer, parameter :: leapday(12)  = (/0,30,60,91,121,152,182,213,244,274,305,335/)

END MODULE ModRamConst

