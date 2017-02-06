SUBROUTINE barycoord_mine_2D (a, b, c, e, alpha, degen)
  !****************************************************************************
  !
  !! BARYCOORD_MINE computes barycentric coordinates of a point in 3D.
  !
  !  Purpose: 
  !
  !    Compute barycentric coordinates of 3D point with respect
  !    to four vertices of a tetrahedron.
  !
  !  Author: Sorin Zaharia
  !
  !  Parameters:
  !
  !    Input, A(1:3),B(1:3),C(1:3),D(1:3) - 4 vertices of tetrahedron.
  !
  !    Input, E(1:3) - fifth point for which barycentric coordinates found
  !
  !    Output, ALPHA(1:3) - scaled barycentric coordinates (if DEGEN = .FALSE.)
  !    such that E = (ALPHA(1)*A + ALPHA(2)*B + ALPHA(3)*C)/DET 
  !    
  !    Output, DEGEN - .TRUE. if A,B,C
  !
  !    Copyright (c) 2016, Los Alamos National Security, LLC
  !    All rights reserved.
  !****************************************************************************


  USE nrtype, ONLY : DP
  !
  IMPLICIT NONE

  !
  REAL(DP) ::  a(2)
  REAL(DP) ::  alpha(3)
  REAL(DP) ::  b(2)
  REAL(DP) ::  c(2)
  REAL(DP) ::  d(3)
  LOGICAL  ::  degen
  REAL(DP) ::  det
  REAL(DP) ::  e(2)
  REAL(DP) ::  tol
  REAL(DP) ::  d1, d2, d3
  !
  tol = 1.E02_dp * EPSILON (tol)
  degen = .FALSE.

  det = a(1)*b(2) + a(2)*c(1) + b(1)*c(2) - b(2)*c(1) - c(2)*a(1) - a(2)*b(1)

  IF (ABS(det) <= 0.01d0*tol*MAX(maxval(abs(a)),maxval(abs(b)),maxval(abs(c)))) THEN
     degen = .TRUE.
     ! PRINT*, 'barycoord_mine_2D: MAX(amax,bmax,cmax), det: ', maxval(abs(a)), maxval(abs(b)), maxval(abs(c)),  det
     RETURN
  END IF

 d1 = e(1)*b(2) + e(2)*c(1) + b(1)*c(2) - b(2)*c(1) - c(2)*e(1) - e(2)*b(1)
 d2 = a(1)*e(2) + a(2)*c(1) + e(1)*c(2) - e(2)*c(1) - a(1)*c(2) - e(1)*a(2)
 d3 = a(1)*b(2) + a(2)*e(1) + b(1)*e(2) - e(1)*b(2) - a(1)*e(2) - a(2)*b(1)


  ! print*, 'barycoord_mine: sum(alpha) = ', sum(alpha)

alpha(1) = d1/det
alpha(2) = d2/det
alpha(3) = d3/det

  RETURN
END SUBROUTINE barycoord_mine_2D
