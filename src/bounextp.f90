SUBROUTINE bounextp

USE Module1

DO  k = 1,nzetp
  DO  j = 2,npsim
    CALL extap(bj(4,j,k),bj(3,j,k),bj(2,j,k),bj(1,j,k))
    CALL extap(bj(nthe-3,j,k),bj(nthe-2,j,k),bj(nthem,j,k) ,bj(nthe,j,k))
    CALL extap(phij(4,j,k),phij(3,j,k),phij(2,j,k),phij(1,j,k))
    CALL extap(phij(nthe-3,j,k),phij(nthe-2,j,k),phij(nthem,j,k)  &
        ,phij(nthe,j,k))

    CALL extap(bf(4,j,k),bf(3,j,k),bf(2,j,k),bf(1,j,k))
    CALL extap(bf(nthe-3,j,k),bf(nthe-2,j,k),bf(nthem,j,k) ,bf(nthe,j,k))

    bsq(1,j,k)=bf(1,j,k)**2
    bsq(nthe,j,k)=bf(nthe,j,k)**2
  
END DO
  DO  i=1,nthe
     CALL extap(bj(i,4,k),bj(i,3,k),bj(i,2,k),bj(i,1,k))
     CALL extap(bj(i,npsi-3,k),bj(i,npsi-2,k),bj(i,npsi-1,k) ,bj(i,npsi,k))
     CALL extap(phij(i,4,k),phij(i,3,k),phij(i,2,k),phij(i,1,k))
     CALL extap(phij(i,npsi-3,k),phij(i,npsi-2,k),phij(i,npsi-1,k)  &
          ,phij(i,npsi,k))
     CALL extap(bf(i,4,k),bf(i,3,k),bf(i,2,k),bf(i,1,k))
     CALL extap(bf(i,npsi-3,k),bf(i,npsi-2,k),bf(i,npsi-1,k) ,bf(i,npsi,k))
     IF (bf(i,npsi,k) < 0._dp) bf(i,npsi,k) = bf(i,npsi-1,k)
     bsq(i,1,k) = bf(i,1,k)**2
     bsq(i,npsi,k) = bf(i,npsi,k)**2
  END DO
END DO

! k0=nZetaMidnight
! DO  k=1,nzeta
!  DO  j=1,npsi
!    DO  i=1,nthe
!c  anisotropic pressure functions
!      tratio(i,j,k)=1.0/(1.-SQRT(bsq(nThetaEquator,j,k0)/bsq(i,j,k))*at0(j))
!      pressure3D(i,j,k)=p(j)*tratio(i,j,k)
!      pressure3D(i,j,k)=pressure3D(i,j,k)*tratio(i,j,k)
!      sigma(i,j,k)=1.0+(pressure3D(i,j,k)-pressure3D(i,j,k))/bsq(i,j,k)
!      tau(i,j,k)=1. -2.*(pressure3D(i,j,k)-pressure3D(i,j,k))/bsq(i,j,k)*tratio(i,j,k)
!      dpper(i,j,k)=(pp(j)  &
!          +p(j)*(tratio(i,j,k)-1.)*dpsibsq(nThetaEquator,j,k)/bsq(nThetaEquator,j,k))  &
!          *tratio(i,j,k)**2

!c  isotropic pressure functions
!      tratio(i,j,k)=1.0
!      pressure3D(i,j,k) = p(j)
!      pressure3D(i,j,k) = p(j)
!      sigma(i,j,k)=1.0
!      tau(i,j,k)=1.0

! New pressure functions
!       tratio(i,j,k)=1.0
!       pressure3D(i,j,k) =  Pressure_rad(sqrt(x(i,j,k)**2 + y(i,j,k)**2))
!       pressure3D(i,j,k)=pressure3D(i,j,k)
!       sigma(i,j,k)=1.0
!       tau(i,j,k)=1.0

!    END DO
!  END DO
! END DO

RETURN

END SUBROUTINE bounextp
