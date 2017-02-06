SUBROUTINE quad2(x1,x2,x3,y1,y2,y3,y,ans,ierr)
!ray  lcm(quad2)
ierr = 0
d1 = (y1-y2)*(y1-y3)
d2 = (y2-y3)*(y2-y1)
d3 = (y3-y1)*(y3-y2)
ans = x1*(y-y2)*(y-y3)/d1 + x2*(y-y3)*(y-y1)/d2  &
    + x3*(y-y1)*(y-y2)/d3

RETURN
END SUBROUTINE quad2
