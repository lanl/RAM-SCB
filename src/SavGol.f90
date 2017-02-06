FUNCTION SavGol(pres)
  USE nrtype, ONLY : DP

  IMPLICIT NONE

  REAL(DP), intent(IN) :: pres(:,:)
  REAL(DP), DIMENSION(SIZE(pres,1),SIZE(pres,2)) :: SavGol, pres1
  REAL(DP) :: BSav(5,5)

  INTEGER :: ier, j, k, nrad, nphi

  nrad = SIZE(pres,1)
  nphi = SIZE(pres,2)


  BSav = RESHAPE((/31, 9, -3, -5, 3, 9, 13, 12, 6, -5, -3, 12, 17, 12, -3, -5, 6, 12, 13, 9, 3, -5, -3, 9, 31/), (/5,5/))
  BSav = BSav/35.

  DO k = 1, nphi
     DO j = 1, nrad
        IF (j > 2 .AND. j < nrad-1) THEN
           pres1(j,k) = DOT_PRODUCT(BSav(3,:), pres(j-2:j+2,k))
        ELSE IF (j == 1) THEN
           pres1(j,k) = DOT_PRODUCT(BSav(1,:), pres(1:5,k))
        ELSE IF (j == 2) THEN
           pres1(j,k) = DOT_PRODUCT(BSav(2,:), pres(1:5,k))
        ELSE IF (j == nrad-1) THEN
           pres1(j,k) = DOT_PRODUCT(BSav(4,:), pres(nrad-4:nrad,k))
        ELSE IF (j == nrad) THEN
           pres1(j,k) = DOT_PRODUCT(BSav(5,:), pres(nrad-4:nrad,k))
        ELSE
           STOP 'SavGol problem.'
        END IF
     END DO
  END DO

  DO j = 1, nrad
     DO k = 1, nphi
        IF (k>2 .AND. k < nphi-1) THEN
           SavGol(j,k) = DOT_PRODUCT(BSav(3,:), pres1(j,k-2:k+2))
        ELSE IF (k == 1 .OR. k == nphi) THEN
           SavGol(j,k) = DOT_PRODUCT(BSav(3,:), (/pres1(j,nphi-2),pres1(j,nphi-1),pres1(j,1),pres1(j,2),pres1(j,3)/))
        ELSE IF (k == 2) THEN
           SavGol(j,k) = DOT_PRODUCT(BSav(3,:), (/pres1(j,nphi-1),pres1(j,1),pres1(j,2),pres1(j,3),pres1(j,4)/))
        ELSE IF (k == nphi-1) THEN
           SavGol(j,k) = DOT_PRODUCT(BSav(3,:), (/pres1(j,nphi-3),pres1(j,nphi-2),pres1(j,nphi-1),pres1(j,1),pres1(j,2)/))
        ELSE
           STOP 'SavGol problem'
        END IF
     END DO
  END DO

  WHERE (SavGol < 0) SavGol = MINVAL(pres)

  RETURN

END FUNCTION SavGol
