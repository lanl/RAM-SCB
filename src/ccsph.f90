SUBROUTINE ccsph(intest, a, b, c, d, e, center, radsq, in)
  !
  !******************************************************************************
  !
  !! CCSPH finds the circumsphere through the vertices of a tetrahedron.
  !
  !
  !  Purpose: 
  !
  !    Find center and square of radius of circumsphere through
  !    four vertices of a tetrahedron, and possibly determine whether
  !    a fifth 3D point is inside sphere.
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
  !    Input, INTEST - .TRUE. iff test for fifth point in sphere to be made.
  !
  !    Input, A(1:3), B(1:3), C(1:3), D(1:3) - 4 vertices of tetrahedron.
  !
  !    Input, E(1:3) - fifth point; referenced if INTEST is .TRUE.
  !
  !    Output, CENTER(1:3) - center of sphere; undefined if A,B,C,D coplanar.
  !
  !    Output, RADSQ - square of radius of sphere; -1 if A,B,C,D coplanar.
  !
  !    Output, IN - contains following value if INTEST is .TRUE.:
  !    2 if A,B,C,D coplanar; 1 if E inside sphere;
  !    0 if E on sphere; -1 if E outside sphere
  !
  
  use nrtype, ONLY : DP
  
  IMPLICIT NONE
  !
  REAL(DP) a(3)
  REAL(DP) b(3)
  REAL(DP) c(3)
  REAL(DP) center(3)
  REAL(DP) cmax
  REAL(DP) cp1
  REAL(DP) cp2
  REAL(DP) cp3
  REAL(DP) d(3)
  REAL(DP) det
  REAL(DP) dsq
  REAL(DP) e(3)
  INTEGER i
  INTEGER in
  LOGICAL intest
  REAL(DP) radsq
  REAL(DP) tol

  REAL(DP) da(3),db(3),dc(3),rhs(3)
  !
  tol = 100.0D+00 * EPSILON (tol)
  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  rhs(1) = 0.5d0*(da(1)**2 + da(2)**2 + da(3)**2)
  rhs(2) = 0.5d0*(db(1)**2 + db(2)**2 + db(3)**2)
  rhs(3) = 0.5d0*(dc(1)**2 + dc(2)**2 + dc(3)**2)

  cmax = MAX ( &
       ABS(a(1)), ABS(a(2)), ABS(a(3)), &
       ABS(b(1)), ABS(b(2)), ABS(b(3)), &
       ABS(c(1)), ABS(c(2)), ABS(c(3)), &
       ABS(d(1)), ABS(d(2)), ABS(d(3)) )

  cp1 = db(2)*dc(3) - dc(2)*db(3)
  cp2 = dc(2)*da(3) - da(2)*dc(3)
  cp3 = da(2)*db(3) - db(2)*da(3)
  det = da(1)*cp1 + db(1)*cp2 + dc(1)*cp3

  IF (ABS(det) <= 0.01d0*tol*cmax) THEN
     radsq = -1.0d0
     in = 2
     RETURN
  END IF

  center(1) = (rhs(1)*cp1 + rhs(2)*cp2 + rhs(3)*cp3)/det
  cp1 = db(1)*rhs(3) - dc(1)*rhs(2)
  cp2 = dc(1)*rhs(1) - da(1)*rhs(3)
  cp3 = da(1)*rhs(2) - db(1)*rhs(1)
  center(2) = (da(3)*cp1 + db(3)*cp2 + dc(3)*cp3)/det
  center(3) = -(da(2)*cp1 + db(2)*cp2 + dc(2)*cp3)/det
  radsq = center(1)**2 + center(2)**2 + center(3)**2

  center(1:3) = center(1:3) + d(1:3)

  IF (intest) THEN
     dsq = SUM ( ( e(1:3) - center(1:3) )**2 )
     IF ( dsq > (1.0d0 + tol) * radsq ) THEN
        in = -1
     ELSE IF ( dsq < (1.0d0 - tol) * radsq ) THEN
        in = 1
     ELSE
        in = 0
     END IF
  END IF

  RETURN
END SUBROUTINE ccsph
