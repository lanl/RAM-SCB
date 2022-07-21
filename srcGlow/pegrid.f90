! Subroutine PEGRID sets up proton energy grid

! This software is part of the GLOW model.  Use is governed by the open source
! academic research license agreement contained in the file Glowlicense.txt.
! For more information see the file Glow.txt.

! Stan Solomon, 2/2018

! Inputs:
!   npbins   number of bins in the proton energy grid
! Outputs:
!   pener    energy at center of each bin, eV
!   pdel     width of each bin, eV

! Energy grid is set up to extend from 37 eV to 803 keV at bin centers(for npbins=27).

subroutine pegrid (pener, pdel, npbins)

  implicit none

  integer,intent(in) :: npbins
  real,intent(out) :: pener(npbins), pdel(npbins)

  integer :: n
  real,parameter :: e0=30.0

  do n=1,npbins
    pener(n) = e0 * exp (0.384*float(n))
  enddo

  pdel(1) = pener(1) - e0
  do n=2,npbins
    pdel(n) = pener(n) - pener(n-1)
  enddo

  do n=1,npbins
    pener(n) = pener(n) - pdel(n)/2.0
  enddo

  return

end subroutine pegrid
