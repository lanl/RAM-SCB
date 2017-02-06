SUBROUTINE quad(x1,x2,x3,y1,y2,y3,y,ans,ierr,deriv)
  ! cray  lcm(quad)
  d1 = (x1-x2)*(x1-x3)
  d2 = (x2-x1)*(x2-x3)
  d3 = (x3-x1)*(x3-x2)
  IF(d1.EQ.0. .OR. d2.EQ.0. .OR. d3.EQ.0.) go to 350
  aco = y1/d1 + y2/d2 + y3/d3
  bco = -(y1*(x2+x3)/d1 + y2*(x1+x3)/d2 + y3*(x1+x2)/d3)
  cco = y1*x2*x3/d1 + y2*x1*x3/d2 + y3*x1*x2/d3 - y
  IF(aco .EQ. 0.) go to 150
  disc = (bco**2 - 4.*aco*cco)
  IF(disc.LT.0) go to 100
  disc = SQRT(disc)
  root1 = (-bco + disc)/(2.*aco)
  root2 = (-bco - disc)/(2.*aco)
  IF(root1.LT.x1 .OR. root1.GT.x3) go to 250
  IF(root2.LT.x1 .OR. root2.GT.x3) go to 260
240 IF(ABS(root1-x2) .LT. ABS(root2-x2)) go to 265
  go to 255
250 CONTINUE
  IF(root2.LT.x1 .OR. root2.GT.x3) go to 240
255 ans = root2
  go to 270
260 CONTINUE
  IF(root1.LT.x1 .OR. root1.GT.x3) go to 240
265 ans = root1
270 CONTINUE
  denom = 2.*aco*ans + bco
  IF(denom.EQ.0.) go to 100
  deriv = 1./denom
  ierr = 0
  RETURN
100 CONTINUE
  ierr = 1
  RETURN
150 CONTINUE
  ans = -cco/bco
  deriv = 1./bco
  RETURN
350 CONTINUE
  ans = 0.
  deriv = 0.
  RETURN
END SUBROUTINE quad
