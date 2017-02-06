SUBROUTINE findx

  ! Obtains equidistant x-s on the midnight equator

  USE nrtype
  USE Module1

  IMPLICIT NONE

  REAL(DP), PARAMETER :: aCS = 1
  REAL(DP) :: deltaX, deltaX2, sum, delX(npsi), w(npsi)
  INTEGER :: j

  deltaX = ABS(x(nThetaEquator, npsi, nZetaMidnight) - &
       x(nThetaEquator, 1, nZetaMidnight)) / (npsi-1)

  currentSheetCase:  IF (iPressureChoice == 10) THEN 
     Isotropic_case:   IF (isotropy == 1) THEN
        CALL pressure
        sum = 0._dp
        DO j = 2, npsi
           sum = sum + 1./ABS(dpdPsi(nThetaEquator,j,nZetaMidnight))
        END DO

        xEqMidNew(1) = x(nThetaEquator,1,nZetaMidnight)

        DO j = 2, npsi
           w(j) = 1./sum * 1./ABS(dpdPsi(nThetaEquator,j,nZetaMidnight))
           delX(j) = w(j) * deltaX * REAL(npsi-1)
           xEqMidNew(j) = xEqMidNew(j-1) - delX(j)
        END DO
     END IF Isotropic_case

     Anisotropic_case:     IF (isotropy == 0) THEN
        CALL computeBandJacob_initial ! to have bf field, needed for computing anisotropic P
        CALL pressure
        sum = 0._dp
        DO j = 2, npsi
           sum = sum + 1./ABS(dPperdPsi(nThetaEquator,j,nZetaMidnight))
        END DO

        ! print*, dPperdPsi(nThetaEquator,:,nZetaMidnight)

        xEqMidNew(1) = x(nThetaEquator,1,nZetaMidnight)

        DO j = 2, npsi
           w(j) = 1./sum * 1./ABS(dPperdPsi(nThetaEquator,j,nZetaMidnight))
           delX(j) = w(j) * deltaX * REAL(npsi-1)
           xEqMidNew(j) = xEqMidNew(j-1) - delX(j)
        END DO
     END IF Anisotropic_case

     !     print*, w(:)*deltaX

  ELSE   ! No current sheet - keep equidistant flux surfaces
     DO j = 2, npsi
        xEqMidNew(j) = x(nThetaEquator, 1, nZetaMidnight) - &
             & REAL(j-1) * deltaX
     END DO
     xEqMidNew(1) = x(nThetaEquator, 1, nZetaMidnight)

  END IF currentSheetCase



  RETURN
END SUBROUTINE findx
