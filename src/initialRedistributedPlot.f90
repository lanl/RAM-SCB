SUBROUTINE initialRedistributedPlot

  USE Module1
  USE ModIoUnit, ONLY: UNITTMP_

  ! First 3 files are for flux plotting 

  OPEN(UNITTMP_, file='plotfile1_initial', status = 'replace', action = 'write')  
  !Noon-midnight plane, day-side
  DO i=1,nThetaEquator
     DO j=1,npsi
        WRITE(UNITTMP_,*) x(i,j,2), z(i,j,2)
     END DO
  END DO
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file='plotfile1neg_initial', status = 'replace', action = 'write')  
  !Noon-midnight plane, night-side
  DO i = 1,nThetaEquator
     DO j = 1,npsi
        WRITE(UNITTMP_,*) x(i, j, nZetaMidnight), z(i, j, nZetaMidnight)  
     END DO
  END DO
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file='plotfile2_initial', status = 'replace', action = 'write')
  DO j=1,npsi
     DO k=2,nZetaMidnight
        WRITE(UNITTMP_,*) x(nThetaEquator,j,k), y(nThetaEquator,j,k)  !Equatorial plane
     END DO
  END DO
  CLOSE(UNITTMP_1)


  ! Next file is for field line plotting

  OPEN(UNITTMP_, file = 'plotline_Tsyg_redistributed', status = 'replace', action = 'write')
  DO k = 2, nZetaMidnight
     DO i = 1, nThetaEquator
        WRITE(UNITTMP_,999) x(i,npsi-1,k), y(i,npsi-1,k), z(i,npsi-1,k)
     END DO
  END DO
999 FORMAT(3f10.5)
  CLOSE(UNITTMP_)


END SUBROUTINE initialRedistributedPlot
