SUBROUTINE baryth ( a, b, c, d, e, alpha, degen )
  !******************************************************************************
  !
  !! BARYTH computes barycentric coordinates of a point in 3D.
  !
  !
  !  Purpose: 
  !
  !    Compute barycentric coordinates of 3D point with respect
  !    to four vertices of a tetrahedron.
  !
  !  Author:
  !
  !    Barry Joe, 
  !    Department of Computing Science, 
  !    University of Alberta,
  !    Edmonton, Alberta, Canada  T6G 2H1
  !    Email: barry@cs.ualberta.ca
  !
  !  Parameters:
  !
  !    Input, A(1:3),B(1:3),C(1:3),D(1:3) - 4 vertices of tetrahedron.
  !
  !    Input, E(1:3) - fifth point for which barycentric coordinates found
  !
  !    Output, ALPHA(1:4) - scaled barycentric coordinates (if DEGEN = .FALSE.)
  !    such that E = (ALPHA(1)*A + ALPHA(2)*B + ALPHA(3)*C +
  !    ALPHA(4)*D)/DET where DET = 6 * (volume of tetra ABCD);
  !    an ALPHA(I) may be set to 0 after tolerance test to
  !    indicate that E is coplanar with a face, so sum of
  !    ALPHA(I)/DET may not be 1; if the actual barycentric
  !    coordinates rather than just their signs are needed,
  !    modify this routine to divide ALPHA(I) by DET.
  !
  !    Output, DEGEN - .TRUE. if A,B,C,D are coplanar.

  USE nrtype, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) ::  a(3)
  REAL(DP) ::  alpha(4)
  REAL(DP) ::  amax
  REAL(DP) ::  b(3)
  REAL(DP) ::  bmax
  REAL(DP) ::  c(3)
  REAL(DP) ::  cmax
  REAL(DP) ::  cp1
  REAL(DP) ::  cp2
  REAL(DP) ::  cp3
  REAL(DP) ::  d(3)
  REAL(DP) ::  da(3)
  REAL(DP) ::  db(3)
  REAL(DP) ::  dc(3)
  REAL(DP) ::  de(3)
  LOGICAL  ::  degen
  REAL(DP) ::  det
  REAL(DP) ::  dmax
  REAL(DP) ::  e(3)
  REAL(DP) ::  ea(3)
  REAL(DP) ::  eb(3)
  REAL(DP) ::  ec(3)
  REAL(DP) ::  emax
  INTEGER  :: i
  REAL(DP) ::  tol
  !
  
  alpha = 0.0_dp
  
  tol = 1.0D+01 * EPSILON (tol)
  degen = .FALSE.

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  amax = MAX(ABS(a(1)),ABS(a(2)),ABS(a(3)))
  bmax = MAX(ABS(b(1)),ABS(b(2)),ABS(b(3)))
  cmax = MAX(ABS(c(1)),ABS(c(2)),ABS(c(3)))
  dmax = MAX(ABS(d(1)),ABS(d(2)),ABS(d(3)))
  cp1 = db(2)*dc(3) - db(3)*dc(2)
  cp2 = db(3)*dc(1) - db(1)*dc(3)
  cp3 = db(1)*dc(2) - db(2)*dc(1)
  det = da(1)*cp1 + da(2)*cp2 + da(3)*cp3

  IF (ABS(det) <= 0.01d0*tol*MAX(amax,bmax,cmax,dmax)) THEN
     degen = .TRUE.
     ! PRINT*, 'baryth: MAX(amax,bmax,cmax,dmax), det: ', amax, bmax, cmax, dmax, det
     RETURN
  END IF

  de(1:3) = e(1:3) - d(1:3)
  ea(1:3) = a(1:3) - e(1:3)
  eb(1:3) = b(1:3) - e(1:3)
  ec(1:3) = c(1:3) - e(1:3)

  alpha(1) = de(1)*cp1 + de(2)*cp2 + de(3)*cp3
  cp1 = da(2)*de(3) - da(3)*de(2)
  cp2 = da(3)*de(1) - da(1)*de(3)
  cp3 = da(1)*de(2) - da(2)*de(1)
  alpha(2) = dc(1)*cp1 + dc(2)*cp2 + dc(3)*cp3
  alpha(3) = db(1)*cp1 + db(2)*cp2 + db(3)*cp3
  alpha(4) = ea(1)*(eb(2)*ec(3) - eb(3)*ec(2)) + ea(2)*(eb(3)*ec(1) &
       - eb(1)*ec(3)) + ea(3)*(eb(1)*ec(2) - eb(2)*ec(1))

  IF (det < 0.0d0) THEN
     alpha(1) = -alpha(1)
     alpha(2) = -alpha(2)
     alpha(4) = -alpha(4)
  ELSE
     alpha(3) = -alpha(3)
  END IF

  emax = MAX(ABS(e(1)),ABS(e(2)),ABS(e(3)))

  IF ( ABS(alpha(1)) <= tol*MAX(bmax,cmax,dmax,emax)) THEN
     alpha(1) = 0.0d0
  END IF

  IF (ABS(alpha(2)) <= tol*MAX(amax,cmax,dmax,emax)) THEN
     alpha(2) = 0.0d0
  END IF

  IF (ABS(alpha(3)) <= tol*MAX(amax,bmax,dmax,emax)) THEN
     alpha(3) = 0.0d0
  END IF

  IF (ABS(alpha(4)) <= tol*MAX(amax,bmax,cmax,emax)) THEN
     alpha(4) = 0.0d0
  END IF

  alpha = - alpha/det  ! to get the real barycentric coordinates

  RETURN
END SUBROUTINE baryth
