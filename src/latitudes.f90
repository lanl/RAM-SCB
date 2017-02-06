!*************************************************
SUBROUTINE latitudes(min_lat, max_lat)
  !*************************************************
  ! const_Z is doing AMR in zeta

  USE nrtype, ONLY : DP, PI_D
  USE Module_points

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: min_lat, max_lat
  INTEGER :: nc, j
  REAL(DP), PARAMETER :: pow = 1.0_dp 
  REAL(DP), EXTERNAL :: ftor
  REAL(DP) :: xzero3, re, bi, psibc
  REAL(DP) :: xzero, bnormal, &
       psiin, psiout, aa, bb, psitot, xpsiin1, xpsiout1, &
       xpsitot, psis, xpl, cost, cost2, thet

  !NAMELIST/eqdat/nthe,npsi,nzeta,xpsiin,xpsiout

  re = radiusStart ! If the latitudes are NOT from Earth's surface, but from the surface of the sphere of radius radiusStart
  xzero = 6.6

  xzero3 = xzero**3

  xpsiin = re / ((COS(pi_d*min_lat/180._dp))**2)
  xpsiout = re / ((COS(pi_d*max_lat/180._dp))**2)

  !IF (rank == 0) WRITE(*, '(A, 5X, 2F8.2)') 'xpsiin, xpsiout = ', xpsiin, xpsiout

  !cc.. bnormal is Earth's dipole magnetic field at xzero in equator (in nT)
  !cc.. for xzero = 6.6 R_E, bnormal = 107.83 nT
  bnormal = 0.31_dp / xzero3 * 1.E5_dp
 
  psiin = -xzero3/xpsiin
  psiout = -xzero3/xpsiout
  psitot = psiout-psiin
  aa = -xzero3
  bb = 0. 

  xpsiout1=xpsiout
  xpsiin1=xpsiin

  xpsitot = xpsiout1 - xpsiin1

  !c  Need to make sure that psival is a monotonically increasing function of j
  !c  define psival grids that correspond to dipole psivals for j=1 and j=nps
  !c  define psival grids that correspond to equal equatorial distance grids in the midnight sector

  !c...  the (i,j,k) coordinate system is related to (theta,psis,zeta) flux coordinate
  !c...   by theta = (i-1) * pi / (nthe-1), 0 < theta < pi
  !c...   by zeta = phi + delta(i,j,k), 0 < zeta < 2*pi
  !c...   by psis = (j-1) / (npsi-1), 0 < psis < 1.0
  !c...      psiin < psival < psiout

  !C  if (rank == 0) then
  !C     WRITE(*,*) 'psiin,psiout,xpsiin,xpsiout,xpsiin1,xpsiout1'
  !C     WRITE(*,*) psiin,psiout,xpsiin,xpsiout,xpsiin1,xpsiout1
  !C     WRITE(*,9991)
  !C  end if
9991 FORMAT('j,psival(j)'//)

  DO  j = 1, npsi
     psis = REAL(j-1, DP) / REAL(npsi-1, DP)
     ! xpl = xpsiin1 + xpsitot * 0.1 * (psis + 9.*SIN(0.5*pi_d*psis))
     xpl = xpsiin1 + xpsitot * psis**pow

     !     xpl = xpsiin + 0.5_dp*xpsitot * (psis + sin(0.5_dp*pi_d*psis))
     psival(j) = aa / xpl + bb * (xpl-re)  ! Actual Euler potential value for this surface
  END DO
  !cc corresponding polar angle on ! the tracing start surface; Earth's surface
 

  re = radiusStart           ! added, to compute latitudes for tracing from radiusStart correctly

  !C  if (rank == 0) WRITE(*,*) 'nc,latitude,psival'

  nc = 0
  DO j = 1, npsi
     psibc = psival(j)
     cost2 = (re * psibc) / (-xzero3)
     cost = SQRT(cost2)
     thet = ACOS(cost)
     nc = nc + 1
     rlatval(nc) = thet * 90. / (0.5 * pi_d)
     !if (rank == 0) WRITE(*,9992) j, rlatval(nc), psival(nc)
9992 FORMAT(I5, 2G15.5)
  END DO

  ! With the new Euler potential choice, the idea is that we can choose any latitudes on the midnight 
  ! (or equivalent line), but then for other local times they are found using the F(beta) function

  RETURN

END SUBROUTINE latitudes

!*****************************
