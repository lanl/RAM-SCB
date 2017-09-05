SUBROUTINE tridag_ser(a,b,c,r,u)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
        REAL(DP), DIMENSION(:), INTENT(OUT) :: u
        REAL(DP), DIMENSION(size(b)) :: gam
        INTEGER(I4B) :: n,j
        REAL(DP) :: bet
        n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
        bet=b(1)
        if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
        u(1)=r(1)/bet
        do j=2,n
                gam(j)=c(j-1)/bet
                bet=b(j)-a(j-1)*gam(j)
                if (bet == 0.0) then
                  write(*,*) j, b(j), a(j-1), gam(j)
                  call nrerror('tridag_ser: Error at code stage 2')
                end if
                u(j)=(r(j)-a(j-1)*u(j-1))/bet
        end do
        do j=n-1,1,-1
                u(j)=u(j)-gam(j+1)*u(j+1)
        end do
        END SUBROUTINE tridag_ser
!===
        RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        USE nrmod, ONLY : tridag_ser
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
        REAL(DP), DIMENSION(:), INTENT(OUT) :: u
        INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
        INTEGER(I4B) :: n,n2,nm,nx
        REAL(DP), DIMENSION(size(b)/2) :: y,q,piva
        REAL(DP), DIMENSION(size(b)/2-1) :: x,z
        REAL(DP), DIMENSION(size(a)/2) :: pivc
        n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
        if (n < NPAR_TRIDAG) then
                call tridag_ser(a,b,c,r,u)
        else
                if (maxval(abs(b(1:n))) == 0.0) &
                        call nrerror('tridag_par: possible singular matrix')
                n2=size(y)
                nm=size(pivc)
                nx=size(x)
                piva = a(1:n-1:2)/b(1:n-1:2)
                pivc = c(2:n-1:2)/b(3:n:2)
                y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
                q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
                if (nm < n2) then
                        y(n2) = b(n)-piva(n2)*c(n-1)
                        q(n2) = r(n)-piva(n2)*r(n-1)
                end if
                x = -piva(2:n2)*a(2:n-2:2)
                z = -pivc(1:nx)*c(3:n-1:2)
                call tridag_par(x,y,z,q,u(2:n:2))
                u(1) = (r(1)-c(1)*u(2))/b(1)
                u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
                        -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
                if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
        end if
        END SUBROUTINE tridag_par

!==============================================================================
SUBROUTINE trapzd(func,a,b,s,n)
  USE nrtype; USE nrutil, ONLY : arth
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP), INTENT(INOUT) :: s
  INTEGER(I4B), INTENT(IN) :: n
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  REAL(DP) :: del,fsum
  INTEGER(I4B) :: it
  IF (n == 1) THEN
     s = 0.5_dp*(b-a)*SUM(func( (/ a, b /) ))
  ELSE
     it = 2**(n-2)
     del=(b-a)/REAL(it, dp)
     fsum = SUM(func(arth(a+0.5_dp*del,del,it)))
     s = 0.5_dp*(s+del*fsum)
  END IF

  RETURN
END SUBROUTINE trapzd

!==============================================================================
        FUNCTION gcf_s(a,x,gln)
        USE nrtype; USE nrutil, ONLY : nrerror
        USE nrmod, ONLY : gammln
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: a,x
        REAL(DP), OPTIONAL, INTENT(OUT) :: gln
        REAL(DP) :: gcf_s
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
        INTEGER(I4B) :: i
        REAL(DP) :: an,b,c,d,del,h
        if (x == 0.0) then
                gcf_s=1.0
                RETURN
        end if
        b=x+1.0_DP-a
        c=1.0_DP/FPMIN
        d=1.0_DP/b
        h=d
        do i=1,ITMAX
                an=-i*(i-a)
                b=b+2.0_DP
                d=an*d+b
                if (abs(d) < FPMIN) d=FPMIN
                c=b+an/c
                if (abs(c) < FPMIN) c=FPMIN
                d=1.0_DP/d
                del=d*c
                h=h*del
                if (abs(del-1.0_DP) <= EPS) exit
        end do
        if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
        if (present(gln)) then
                gln=gammln(a)
                gcf_s=exp(-x+a*log(x)-gln)*h
        else
                gcf_s=exp(-x+a*log(x)-gammln(a))*h
        end if
        END FUNCTION gcf_s

!===
        FUNCTION gcf_v(a,x,gln)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        USE nrmod, ONLY : gammln
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
        REAL(DP), DIMENSION(size(a)) :: gcf_v
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
        INTEGER(I4B) :: i
        REAL(DP), DIMENSION(size(a)) :: an,b,c,d,del,h
        LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
        i=assert_eq(size(a),size(x),'gcf_v')
        zero=(x == 0.0)
        where (zero)
                gcf_v=1.0
        elsewhere
                b=x+1.0_DP-a
                c=1.0_DP/FPMIN
                d=1.0_DP/b
                h=d
        end where
        converged=zero
        do i=1,ITMAX
                where (.not. converged)
                        an=-i*(i-a)
                        b=b+2.0_DP
                        d=an*d+b
                        d=merge(FPMIN,d, abs(d)<FPMIN )
                        c=b+an/c
                        c=merge(FPMIN,c, abs(c)<FPMIN )
                        d=1.0_DP/d
                        del=d*c
                        h=h*del
                        converged = (abs(del-1.0_DP)<=EPS)
                end where
                if (all(converged)) exit
        end do
        if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
        if (present(gln)) then
                if (size(gln) < size(a)) call &
                        nrerror('gser: Not enough space for gln')
                gln=gammln(a)
                where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
        else
                where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
        end if
        END FUNCTION gcf_v

!==============================================================================
        FUNCTION gser_s(a,x,gln)
        USE nrtype; USE nrutil, ONLY : nrerror
        USE nrmod, ONLY : gammln
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: a,x
        REAL(DP), OPTIONAL, INTENT(OUT) :: gln
        REAL(DP) :: gser_s
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(DP), PARAMETER :: EPS=epsilon(x)
        INTEGER(I4B) :: n
        REAL(DP) :: ap,del,summ
        if (x == 0.0) then
                gser_s=0.0
                RETURN
        end if
        ap=a
        summ=1.0_DP/a
        del=summ
        do n=1,ITMAX
                ap=ap+1.0_DP
                del=del*x/ap
                summ=summ+del
                if (abs(del) < abs(summ)*EPS) exit
        end do
        if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
        if (present(gln)) then
                gln=gammln(a)
                gser_s=summ*exp(-x+a*log(x)-gln)
        else
                gser_s=summ*exp(-x+a*log(x)-gammln(a))
        end if
        END FUNCTION gser_s
!===
        FUNCTION gser_v(a,x,gln)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        USE nrmod, ONLY : gammln
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
        REAL(DP), DIMENSION(size(a)) :: gser_v
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(DP), PARAMETER :: EPS=epsilon(x)
        INTEGER(I4B) :: n
        REAL(DP), DIMENSION(size(a)) :: ap,del,summ
        LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
        n=assert_eq(size(a),size(x),'gser_v')
        zero=(x == 0.0)
        where (zero) gser_v=0.0
        ap=a
        summ=1.0_DP/a
        del=summ
        converged=zero
        do n=1,ITMAX
                where (.not. converged)
                        ap=ap+1.0_DP
                        del=del*x/ap
                        summ=summ+del
                        converged = (abs(del) < abs(summ)*EPS)
                end where
                if (all(converged)) exit
        end do
        if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
        if (present(gln)) then
                if (size(gln) < size(a)) call &
                        nrerror('gser: Not enough space for gln')
                gln=gammln(a)
                where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
        else
                where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
        end if
        END FUNCTION gser_v

!==============================================================================
        SUBROUTINE gaussj(a,b)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
        LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
        REAL(DP) :: pivinv
        REAL(DP), DIMENSION(size(a,1)) :: dumc
        INTEGER(I4B), TARGET :: irc(2)
        INTEGER(I4B) :: i,l,n
        INTEGER(I4B), POINTER :: irow,icol
        n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
        irow => irc(1)
        icol => irc(2)
        ipiv=0
        do i=1,n
                lpiv = (ipiv == 0)
                irc=maxloc(abs(a),outerand(lpiv,lpiv))
                ipiv(icol)=ipiv(icol)+1
                if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
                if (irow /= icol) then
                        call swap(a(irow,:),a(icol,:))
                        call swap(b(irow,:),b(icol,:))
                end if
                indxr(i)=irow
                indxc(i)=icol
                if (a(icol,icol) == 0.0) &
                        call nrerror('gaussj: singular matrix (2)')
                pivinv=1.0_DP/a(icol,icol)
                a(icol,icol)=1.0
                a(icol,:)=a(icol,:)*pivinv
                b(icol,:)=b(icol,:)*pivinv
                dumc=a(:,icol)
                a(:,icol)=0.0
                a(icol,icol)=pivinv
                a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
                b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
                a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
                b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
        end do
        do l=n,1,-1
                call swap(a(:,indxr(l)),a(:,indxc(l)))
        end do
        END SUBROUTINE gaussj

!==============================================================================
        FUNCTION gammp_s(a,x)
        USE nrtype; USE nrutil, ONLY : assert
        USE nrmod, ONLY : gcf,gser
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: a,x
        REAL(DP) :: gammp_s
        call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
        if (x<a+1.0_DP) then
                gammp_s=gser(a,x)
        else
                gammp_s=1.0_DP-gcf(a,x)
        end if
        END FUNCTION gammp_s
!===
        FUNCTION gammp_v(a,x)
        USE nrtype; USE nrutil, ONLY : assert,assert_eq
        USE nrmod, ONLY : gcf,gser
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(DP), DIMENSION(size(x)) :: gammp_v
        LOGICAL(LGT), DIMENSION(size(x)) :: mask
        INTEGER(I4B) :: ndum
        ndum=assert_eq(size(a),size(x),'gammp_v')
        call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
        mask = (x<a+1.0_DP)
        gammp_v=merge(gser(a,merge(x,0.0_DP,mask)), &
                1.0_DP-gcf(a,merge(x,0.0_DP,.not. mask)),mask)
        END FUNCTION gammp_v

!==============================================================================
        FUNCTION gammln_s(xx)
        USE nrtype; USE nrutil, ONLY : arth,assert
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: xx
        REAL(DP) :: gammln_s
        REAL(DP) :: tmp,x
        REAL(DP) :: stp = 2.5066282746310005_dp
        REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                -86.50532032941677_dp,24.01409824083091_dp,&
                -1.231739572450155_dp,0.1208650973866179e-2_dp,&
                -0.5395239384953e-5_dp/)
        call assert(xx > 0.0, 'gammln_s arg')
        x=xx
        tmp=x+5.5_dp
        tmp=(x+0.5_dp)*log(tmp)-tmp
        gammln_s=tmp+log(stp*(1.000000000190015_dp+&
                sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
        END FUNCTION gammln_s
!===
        FUNCTION gammln_v(xx)
        USE nrtype; USE nrutil, ONLY: assert
        IMPLICIT NONE
        INTEGER(I4B) :: i
        REAL(DP), DIMENSION(:), INTENT(IN) :: xx
        REAL(DP), DIMENSION(size(xx)) :: gammln_v
        REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
        REAL(DP) :: stp = 2.5066282746310005_dp
        REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                -86.50532032941677_dp,24.01409824083091_dp,&
                -1.231739572450155_dp,0.1208650973866179e-2_dp,&
                -0.5395239384953e-5_dp/)
        if (size(xx) == 0) RETURN
        call assert(all(xx > 0.0), 'gammln_v arg')
        x=xx
        tmp=x+5.5_dp
        tmp=(x+0.5_dp)*log(tmp)-tmp
        ser=1.000000000190015_dp
        y=x
        do i=1,size(coef)
                y=y+1.0_dp
                ser=ser+coef(i)/y
        end do
        gammln_v=tmp+log(stp*ser/x)
        END FUNCTION gammln_v

!==============================================================================
SUBROUTINE polint(xa,ya,x,y,dy)
  USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(OUT) :: y, dy
  INTEGER(I4B) :: m,n,ns
  REAL(DP), DIMENSION(SIZE(xa)) :: c,d,den,ho
  n = assert_eq(SIZE(xa),SIZE(ya),'polint')
  c=ya
  d=ya
  ho=xa-x
  ns=iminloc(ABS(x-xa))
  y=ya(ns)
  ns=ns-1
  DO m=1,n-1
     den(1:n-m)=ho(1:n-m)-ho(1+m:n)
     IF (ANY(den(1:n-m) == 0.0)) THEN
        PRINT*, 'xa = ', xa
        PRINT*, 'ya = ', ya
        PRINT*, 'x = ', x
        PRINT*, 'y = ', y
        PAUSE 'Problem'
        CALL nrerror('polint: calculation failure')
     END IF
     den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
     d(1:n-m)=ho(1+m:n)*den(1:n-m)
     c(1:n-m)=ho(1:n-m)*den(1:n-m)
     IF (2*ns < n-m) THEN
        dy=c(ns+1)
     ELSE
        dy=d(ns)
        ns=ns-1
     END IF
     y=y+dy
  END DO

return
END SUBROUTINE polint

!==============================================================================
FUNCTION qtrap(func,a,b)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nrmod, ONLY : trapzd
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP) :: qtrap
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: JMAX=20
  REAL(DP), PARAMETER :: EPS=1.0e-6_DP, UNLIKELY=-1.0e30_DP
  REAL(DP) :: olds
  INTEGER(I4B) :: j
  olds=UNLIKELY
  DO j=1, JMAX
     CALL trapzd(func,a,b,qtrap,j)
     ! print*, 'qtrap: j, qtrap = ', j, qtrap
     IF (ABS(qtrap-olds) < EPS*ABS(olds)) RETURN
     IF (qtrap == 0.0 .AND. olds == 0.0 .AND. j > 6) RETURN
     olds=qtrap
  END DO
  CALL nrerror('qtrap: too many steps')
END FUNCTION qtrap

!==============================================================================
FUNCTION qromb(func,a,b)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nrmod, ONLY : polint, trapzd
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a, b
  REAL(DP) :: qromb
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: JMAX=50, JMAXP=JMAX+1, K=2, KM=K-1 ! With K = 2 this is equivalent to Simpson rule
  REAL(DP), PARAMETER :: EPS = 1.0e-3_dp
  REAL(DP), DIMENSION(JMAXP) :: h, s
  REAL(DP) :: dqromb
  INTEGER(I4B) :: j
  h(1) = 1.0_dp
  DO j = 1, JMAX
     CALL trapzd(func,a,b,s(j),j)
     ! print*, 'qromb: j, s(j) = ', j, s(j)
     IF (j >= K) THEN
        CALL polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
        IF (ABS(dqromb) <= EPS*ABS(qromb)) RETURN
     END IF
     s(j+1)=s(j)
     h(j+1)=0.25_dp*h(j)
  END DO
  CALL nrerror('qromb: too many steps')
END FUNCTION qromb

!==============================================================================
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

!==============================================================================
FUNCTION locate(xx,x)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  INTEGER(I4B) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=SIZE(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  DO
     IF (ju-jl <= 1) EXIT
     jm=(ju+jl)/2
     IF (ascnd .EQV. (x >= xx(jm))) THEN
        jl=jm
     ELSE
        ju=jm
     END IF
  END DO
  IF (x == xx(1)) THEN
     locate=1
  ELSE IF (x == xx(n)) THEN
     locate=n-1
  ELSE
     locate=jl
  END IF
END FUNCTION locate

