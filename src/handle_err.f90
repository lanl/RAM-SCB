SUBROUTINE HANDLE_ERR(STATUS)

  USE netcdf

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: STATUS

  IF(status /= nf90_noerr) THEN
     PRINT *, TRIM(nf90_strerror(status))
     STOP "Stopped by handle_err."
  END IF
END SUBROUTINE HANDLE_ERR

