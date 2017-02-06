SUBROUTINE MagneticData

  USE Module1
  USE ModIoUnit, ONLY: UNITTMP_

  OPEN(UNITTMP_, file='plotfile1', status = 'replace', action = 'write')  
  ! Noon-midnight plane, day-side
  DO i=1,nThetaEquator
     DO j=1,npsi
        WRITE(UNITTMP_,*) x(i,j,2), z(i,j,2)
     END DO
  END DO
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file='plotfile1neg', status = 'replace', action = 'write')  
  !Noon-midnight plane, night-side
  DO i=1,nThetaEquator
     DO j=1,npsi
        WRITE(UNITTMP_,*) x(i,j,nZetaMidnight), z(i,j,nZetaMidnight)  
     END DO
  END DO
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file='plotfile2', status = 'replace', action = 'write')
  DO j=1,npsi
     DO k=2,nZetaMidnight
        WRITE(UNITTMP_,*) x(nThetaEquator,j,k), y(nThetaEquator,j,k)  ! Equatorial plane
     END DO
  END DO
  CLOSE(UNITTMP_)


  OPEN(UNITTMP_, file = 'plotline', status = 'replace', action = 'write')
  DO k = 2, nZetaMidnight
     DO i = 1, nThetaEquator
        WRITE(UNITTMP_,999) x(i,npsi,k), y(i,npsi,k), z(i,npsi,k)
     END DO
  END DO
999 FORMAT(3f10.5)
  CLOSE(UNITTMP_)

  OPEN(UNITTMP_, file = 'plotline_inner', status = 'replace', action = 'write')
  DO k = 2, nZetaMidnight
     DO i = 1, nThetaEquator
        WRITE(UNITTMP_,999) x(i,40,k), y(i,40,k), z(i,40,k)
     END DO
  END DO
  CLOSE(UNITTMP_)

 OPEN(UNITTMP_, file = 'plotline2', status = 'replace', action = 'write')
  DO k = nZetaMidnight+1, nzeta
     DO i = 1, nThetaEquator
        WRITE(UNITTMP_,999) x(i,npsi-1,k), y(i,npsi-1,k), z(i,npsi-1,k)
     END DO
  END DO

  CLOSE(UNITTMP_)

  RETURN

END SUBROUTINE MagneticData
