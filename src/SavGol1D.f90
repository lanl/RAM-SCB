SUBROUTINE SavGol1D(flux, NT, iorder)
  ! Savitzky-Golay smoothing (1-D) using quadratic polynomial and 5 pts

  IMPLICIT NONE

  REAL, INTENT(IN OUT) :: flux(NT)
  INTEGER, INTENT(IN) :: NT, iorder

  REAL :: flux1(NT)
  REAL :: BSav(iorder, iorder)

  INTEGER :: ier, iPass, j, k, nphi

  nphi = NT
!C  PRINT*, 'SavGol1D: iorder, nphi, size(flux) = ', iorder, nphi, SIZE(flux)
  IF (iorder == 5) THEN
     BSav = RESHAPE((/31, 9, -3, -5, 3, 9, 13, 12, 6, -5, -3, 12, 17, 12, -3, -5, 6, 12, 13, 9, 3, -5, -3, 9, 31/), (/5,5/))
     BSav = BSav/35.
  ELSE IF (iorder == 7) THEN
     BSav = RESHAPE((/32, 5, 1, -2, -2, -1, 5, &
          15, 4, 3, 3, 1, 0, -3, &
          3, 3, 4, 6, 3, 1, -6, &
          -4, 2, 4, 7, 4, 2, -4, &
          -6, 1, 3, 6, 4, 3, 3, &
          -3, 0, 1, 3, 3, 4, 15, &
          5, -1, -2, -2, 1, 5, 32/), (/7,7/))
     ! Reshapes in Fortran normal order (storage by columns)
     ! Here all are centered and there is a common weight of 21 (Table I of [Gorry, 1990])
     ! Pass in the k (phi) direction
     BSav = Bsav/21.
  ELSE 
     STOP 'SavGol1D: problem, invalid choice'
  END IF

  !C  print*, 'SavGol1D: BSav(3,:), SUM(BSav(3,:)) = ', BSav(3,:), SUM(BSav(3,:))
  !C IF (iorder == 7) PRINT*, 'SavGol1D: BSav(4,:), SUM(BSav(4,:)) = ', BSav(4,:), SUM(BSav(4,:))

  DO iPass = 1, 3 ! 3 passes

     IF (iorder == 5) THEN
        DO k = 1, nphi
           IF (k > 2 .AND. k < nphi-1) THEN
              flux1(k) = DOT_PRODUCT(BSav(3,:), flux(k-2:k+2))
           ELSE IF (k == 1 .OR. k == nphi) THEN
              flux1(k) = DOT_PRODUCT(BSav(3,:), (/flux(nphi-2), &
                   flux(nphi-1),flux(1),flux(2),flux(3)/))
           ELSE IF (k == 2) THEN
              flux1(k) = DOT_PRODUCT(BSav(3,:), (/flux(nphi-1),flux(1),flux(2),flux(3),flux(4)/))
           ELSE IF (k == nphi-1) THEN
              flux1(k) = DOT_PRODUCT(BSav(3,:), (/flux(nphi-3),&
                   flux(nphi-2),flux(nphi-1),flux(1),flux(2)/))
           ELSE
              STOP 'SavGol1D problem'
           END IF
        END DO
     ELSE
        DO k = 1, nphi
           IF (k > 3 .AND. k < nphi-2) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), flux(k-3:k+3))
           ELSE IF (k == 1 .OR. k == nphi) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), (/flux(nphi-3),flux(nphi-2), &
                   flux(nphi-1),flux(1),flux(2),flux(3),flux(4)/))
           ELSE IF (k == 2) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), (/flux(nphi-2),flux(nphi-1),flux(1),flux(2),&
                   flux(3),flux(4),flux(5)/))
           ELSE IF (k == 3) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), (/flux(nphi-1),flux(1),flux(2),flux(3),&
                   flux(4),flux(5),flux(6)/))
           ELSE IF (k == nphi-2) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), (/flux(nphi-5),flux(nphi-4),&
                   flux(nphi-3),flux(nphi-2),flux(1),flux(2),flux(3)/))
           ELSE IF (k == nphi-1) THEN
              flux1(k) = DOT_PRODUCT(BSav(4,:), (/flux(nphi-4),flux(nphi-3),&
                   flux(nphi-2),flux(nphi-1),flux(1),flux(2),flux(3)/))
           ELSE
              STOP 'SavGol1D problem'
           END IF
        END DO
     END IF

     WHERE (flux1 < 0) flux1 = MINVAL(flux)

     flux = flux1

  END DO

  RETURN

END SUBROUTINE SavGol1D
