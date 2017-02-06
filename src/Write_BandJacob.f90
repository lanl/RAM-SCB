SUBROUTINE Write_BandJacob

  USE Module1
  USE ezcdf

  IMPLICIT NONE

  INTEGER, DIMENSION(3) :: dimlens = (/1, 1, 1/)
  INTEGER :: ierr, ncdfId

  CALL cdf_open(ncdfId, 'mag_field.cdf', 'w')
  dimlens(1) = npsi
  dimlens(2) = nzeta
  dimlens(3) = 1
  CALL cdf_define(ncdfId, 'xEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'yEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'BEq', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'JacobianEq', dimlens, 'R8')

  dimlens(1) = nthe
  dimlens(2) = npsi
  dimlens(3) = 1
  CALL cdf_define(ncdfId, 'xMid', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'zMid', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'BMid', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'JacobianMid', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'xNoo', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'zNoo', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'BNoo', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'JacobianNoo', dimlens, 'R8')

  dimlens(1) = nthe
  dimlens(2) = npsi
  dimlens(3) = nzeta
  CALL cdf_define(ncdfId, 'xFull', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'yFull', dimlens, 'R8')
  CALL cdf_define(ncdfId, 'zFull', dimlens, 'R8')

  IF (ALLOCATED(bX)) THEN
     CALL cdf_define(ncdfId, 'bXFull', dimlens, 'R8')
     CALL cdf_define(ncdfId, 'bYFull', dimlens, 'R8')
     CALL cdf_define(ncdfId, 'bZFull', dimlens, 'R8')
     CALL cdf_define(ncdfId, 'bFull', dimlens, 'R8')
  END IF

  CALL cdf_write(ncdfId, 'xEq', x(nThetaEquator,:,1:nzeta))
  CALL cdf_write(ncdfId, 'yEq', y(nThetaEquator,:,1:nzeta))
  CALL cdf_write(ncdfId, 'BEq', bnormal*bf(nThetaEquator,1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'JacobianEq', jacobian(nThetaEquator,1:npsi,1:nzeta))

  CALL cdf_write(ncdfId, 'xMid', x(:,:,nZetaMidnight))
  CALL cdf_write(ncdfId, 'zMid', z(:,:,nZetaMidnight))
  CALL cdf_write(ncdfId, 'BMid', bnormal*bf(1:nthe,1:npsi,nZetaMidnight))
  CALL cdf_write(ncdfId, 'JacobianMid', jacobian(1:nthe,1:npsi,nZetaMidnight))
  CALL cdf_write(ncdfId, 'xNoo', x(:,:,2))
  CALL cdf_write(ncdfId, 'zNoo', z(:,:,2))
  CALL cdf_write(ncdfId, 'BNoo', bnormal*bf(1:nthe,1:npsi,2))
  CALL cdf_write(ncdfId, 'JacobianNoo', jacobian(1:nthe,1:npsi,2))

  CALL cdf_write(ncdfId, 'xFull', x(1:nthe,1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'yFull', y(1:nthe,1:npsi,1:nzeta))
  CALL cdf_write(ncdfId, 'zFull', z(1:nthe,1:npsi,1:nzeta))
  
  IF (ALLOCATED(bX)) THEN
     CALL cdf_write(ncdfId, 'bFull', bnormal*bf(1:nthe,1:npsi,1:nzeta))
     CALL cdf_write(ncdfId, 'bXFull', bnormal*bX(1:nthe,1:npsi,1:nzeta))
     CALL cdf_write(ncdfId, 'bYFull', bnormal*bY(1:nthe,1:npsi,1:nzeta))
     CALL cdf_write(ncdfId, 'bZFull', bnormal*bZ(1:nthe,1:npsi,1:nzeta))
  END IF

  CALL cdf_close(ncdfId)

  RETURN
END SUBROUTINE Write_BandJacob
