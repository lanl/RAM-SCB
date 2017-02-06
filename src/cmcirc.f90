FUNCTION cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )
  !
  !******************************************************************************
  !
  !! CMCIRC determines if a point is in the circumcircle of three points.
  !
  !
  !  Purpose: 
  !
  !    Determine whether (X0,Y0) is in the circumcircle through
  !    the three points (X1,Y1), (X2,Y2), (X3,Y3).
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
  !    Input, X0, Y0, X1, Y1, X2, Y2, X3, Y3 - vertex coordinates.
  !
  !    Output, CMCIRC -  
  !     2 if three vertices collinear,
  !     1 if (X0,Y0) inside circle,
  !     0 if (X0,Y0) on circle,
  !    -1 if (X0,Y0) outside circle
  !
  
use nrtype, ONLY : DP

IMPLICIT NONE
  !
  REAL (DP) :: a11
  REAL (DP) :: a12
  REAL (DP) ::  a21
  REAL (DP) ::  a22
  REAL (DP) ::  b1
  REAL (DP) ::  b2
  INTEGER cmcirc
  REAL (DP) ::  det
  REAL (DP) ::  diff
  REAL (DP) ::  rsq
  REAL (DP) ::  tol
  REAL (DP) ::  tolabs
  REAL (DP) ::  x0
  REAL (DP) ::  x1
  REAL (DP) ::  x2
  REAL (DP) ::  x3
  REAL (DP) ::  xc
  REAL (DP) ::  y0
  REAL (DP) ::  y1
  REAL (DP) ::  y2
  REAL (DP) ::  y3
  REAL (DP) ::   yc
  !
  tol = 100.0D+00 * EPSILON ( tol )
  cmcirc = 2
  a11 = x2 - x1
  a12 = y2 - y1
  a21 = x3 - x1
  a22 = y3 - y1
  tolabs = tol*MAX(ABS(a11),ABS(a12),ABS(a21),ABS(a22))
  det = a11*a22 - a21*a12
  IF (ABS(det) <= tolabs) RETURN
  b1 = a11**2 + a12**2
  b2 = a21**2 + a22**2
  det = det + det
  xc = (b1*a22 - b2*a12)/det
  yc = (b2*a11 - b1*a21)/det
  rsq = xc**2 + yc**2
  diff = ((x0 - x1 - xc)**2 + (y0 - y1 - yc)**2) - rsq
  tolabs = tol*rsq

  IF (diff < -tolabs) THEN
     cmcirc = 1
  ELSE IF (diff > tolabs) THEN
     cmcirc = -1
  ELSE
     cmcirc = 0
  END IF

  RETURN
END FUNCTION cmcirc
