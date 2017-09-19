
Module ModConductance
  use ModKind
  implicit none
  save

  integer :: longmx, latmx, ndx
  real(real8_)    :: steplat

  integer, parameter :: latmdx = 40
  integer, parameter :: lonmdx = 30
  integer, parameter :: mndx   = 10
  integer, parameter :: nConductanceSolutions = 4

  real(real8_), dimension(mndx) :: halmin, pedmin, avk50, efx50

  real(real8_), dimension(0:lonmdx,0:latmdx,mndx) :: &
        halar, pedar, avkar, efxar

  real(real8_), dimension(nConductanceSolutions) :: ConductanceBackground

  real(real8_), dimension(4) :: bkgc, bkgcer
  real(real8_)               :: flxmax, dbymax, dbzmax

  integer, parameter :: pedersen_ = 1
  integer, parameter :: hall_     = 2
  integer, parameter :: eflux_    = 3
  integer, parameter :: avee_     = 4

end Module ModConductance
