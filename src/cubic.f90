SUBROUTINE cubic(x1,x2,x3,x4,y1,y2,y3,y4,y,ans,dans,is,ierr)
  ! ray  lcm(cubic)
  ! cubic polynomial interpolation
  ! this routine returns interpolated value and first derivative

  d1 = (y1-y2)*(y1-y3)*(y1-y4)
  d2 = (y2-y1)*(y2-y3)*(y2-y4)
  d3 = (y3-y1)*(y3-y2)*(y3-y4)
  d4 = (y4-y1)*(y4-y2)*(y4-y3)
  ans = x1*(y-y2)*(y-y3)*(y-y4)/d1 + x2*(y-y1)*(y-y3)*(y-y4)/d2  &
       + x3*(y-y1)*(y-y2)*(y-y4)/d3 + x4*(y-y1)*(y-y2)*(y-y3)/d4
  
  IF ( is == 0 )  RETURN
  IF(is /= 2) THEN
     dans = x1/d1*((y-y3)*(y-y4)+(y-y2)*(y-y4)+(y-y2)*(y-y3))  &
          + x2/d2*((y-y3)*(y-y4)+(y-y1)*(y-y4)+(y-y1)*(y-y3))  &
          + x3/d3*((y-y2)*(y-y4)+(y-y1)*(y-y4)+(y-y1)*(y-y2))  &
          + x4/d4*((y-y2)*(y-y3)+(y-y1)*(y-y3)+(y-y1)*(y-y2))
     RETURN
  END IF
  
  dans = x1/d1*((y-y3)+(y-y4)+(y-y2)+(y-y4)+(y-y2)+(y-y3))  &
       + x2/d2*((y-y3)+(y-y4)+(y-y1)+(y-y4)+(y-y1)+(y-y3))  &
       + x3/d3*((y-y2)+(y-y4)+(y-y1)+(y-y4)+(y-y1)+(y-y2))  &
       + x4/d4*((y-y2)+(y-y3)+(y-y1)+(y-y3)+(y-y1)+(y-y2))
  RETURN
  
  ! 20 CONTINUE
  anum = (y2-y3)*(x1-x2)+(y1-y2)*(x3-x2)
  adem = (y1-y3)*(y1-y2)**2
  a = anum/adem
  b = (x3-x1)/(y3-y1) - 2.*y2*a
  dans = 2.*a*y + b
  RETURN
END SUBROUTINE cubic
