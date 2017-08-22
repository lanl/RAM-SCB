MODULE ModRamConst

  use ModRamMain,  ONLY: Real8_
  use ModNumConst, ONLY: cPi

  real(kind=Real8_) :: ME   = 7.9E15,    &   ! Magnetic moment of the earth [T*m3]
                       RE   = 6.371E6,   &   ! Earth's radius [m]
                       HMIN = 2E5,       &   ! Altitude of dense atmosphere [m]
                       MP   = 1.673E-27, &   ! Mass of H+ [kg]
                       Q    = 1.602E-19, &   ! Elementary charge [C]
                       CS   = 2.998E8,   &   ! Speed of light [m/s]
                       PI   = cPi            ! pi used to be 3.141592654

  real(kind=Real8_), dimension(4) :: M1 = (/5.4462E-4,1.,4.,16./) ! Mass number of e-, H+, He+ and O+

  integer, parameter :: monthday(12) = (/0,30,59,90,120,151,181,212,243,273,304,334/)
  integer, parameter :: leapday(12)  = (/0,30,60,91,121,152,182,213,244,274,305,335/)


END MODULE ModRamConst

