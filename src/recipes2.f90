!=============================================================================
subroutine LINTP(XX,YY,N,X,Y,IER)

  use ModRamMain, ONLY: Real8_
  implicit none
  
  ! Arguments:
  integer, intent(in)  :: N
  integer, intent(out) :: ier
  real(kind=Real8_), intent(in), dimension(N) :: XX, YY
  real(kind=Real8_), intent(in)               :: X
  real(kind=Real8_), intent(out)              :: Y

  ! Internal variables
  integer           :: JL, JM, JU, J
  real(kind=Real8_) :: D
  !-----------------------------------------------------------------------

  IER = 0

  ! Initialize upper and lower boundaries.
  JL=1
  JU=N

  ! if not dne compute a midpoint

  do while (JU-JL.GT.1)
     JM=(JU+JL)/2
     ! now replace lower or upper limit
     IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
        JL=JM
     ELSE
        JU=JM
     ENDIF
     ! try again
  end do

  J=JL
  D=XX(J+1)-XX(J)
  Y=(YY(J)*(XX(J+1)-X)+YY(J+1)*(X-XX(J)))/D

  RETURN
end subroutine LINTP

!=============================================================================
subroutine LINTP2(X1A,X2A,YA,M,N,X1,X2,Y,IER)

  use ModRamMain, ONLY: Real8_
  implicit none

  ! Arguments:
  integer, intent(in) :: m, n
  integer, intent(out):: ier
  real(kind=Real8_), intent(in), dimension(n)  :: x2a
  real(kind=Real8_), intent(in), dimension(m)  :: x1a
  real(kind=Real8_), intent(in), dimension(m,n):: ya
  real(kind=Real8_), intent(in)  :: x1, x2
  real(kind=Real8_), intent(out) :: y

  ! Local variables:
  integer :: j, k
  real(kind=Real8_) :: ytmp(N), yytmp(M)

  !-----------------------------------------------------------------------
  IER = 0

  ! do M evaluations of row constructed using 1-dimensional evaluator LINTP
  do J=1,M
     do K=1,N
        YTMP(K)=YA(J,K)
     end do
     call lintp(X2A,YTMP,N,X2,YYTMP(J),IER)
     IER = 10 * IER
     if (IER.EQ.10) return
  end do

  ! evaluate it
  CALL LINTP(X1A,YYTMP,M,X1,Y,IER)
  IER = IER * 10
  return
end subroutine LINTP2

!=============================================================================
! Need to pass interpolation routine indices of neigbhors, X1ind & X2ind
subroutine LINTP2_ind(X1A,X2A,YA,M,N,X1,X2,X1ind,X2ind,Y,IER)

  use ModRamMain, ONLY: Real8_
  implicit none

  ! Arguments:
  integer, intent(in) :: m, n
  integer, intent(out):: ier
  real(kind=Real8_), intent(in), dimension(n)  :: x2a
  real(kind=Real8_), intent(in), dimension(m)  :: x1a
  real(kind=Real8_), intent(in), dimension(m,n):: ya
  real(kind=Real8_), intent(in)  :: x1, x2
  integer, intent(in), dimension(2) :: X1ind,X2ind
  real(kind=Real8_), intent(out) :: y

  !-----------------------------------------------------------------------
  IER = 0

  Y = YA(X1ind(1),X2ind(1))*(X1A(X1ind(2))-X1)*(X2A(X2ind(2))-X2) + &
       YA(X1ind(2),X2ind(1))*(X1-X1A(X1ind(1)))*(X2A(X2ind(2))-X2) + &
       YA(X1ind(1),X2ind(2))*(X1A(X1ind(2))-X1)*(X2-X2A(X2ind(1))) + &
       YA(X1ind(2),X2ind(2))*(X1-X1A(X1ind(1)))*(X2-X2A(X2ind(1)))

  Y = Y/(X1A(X1ind(2))-X1A(X1ind(1)))/(X2A(X2ind(2))-X2A(X2ind(1)))
  
  return
end subroutine LINTP2_IND

!=============================================================================
! Need to pass interpolation routine indices of neigbhors, X1ind & X2ind AND coefficients Xnncoef
subroutine LINTP2_ind_coef(X1A,X2A,YA,M,N,X1,X2,X1ind,X2ind,Xnncoef,Y,IER)

  use ModRamMain, ONLY: Real8_
  implicit none

  ! Arguments:
  integer, intent(in) :: m, n
  integer, intent(out):: ier
  real(kind=Real8_), intent(in), dimension(n)  :: x2a
  real(kind=Real8_), intent(in), dimension(m)  :: x1a
  real(kind=Real8_), intent(in), dimension(m,n):: ya
  real(kind=Real8_), intent(in), dimension(2,2):: Xnncoef
  real(kind=Real8_), intent(in)  :: x1, x2
  integer, intent(in), dimension(2) :: X1ind,X2ind
  real(kind=Real8_), intent(out) :: y

  !-----------------------------------------------------------------------
  IER = 0

  Y = YA(X1ind(1),X2ind(1))*Xnncoef(1,1) + YA(X1ind(2),X2ind(1))*Xnncoef(2,1) + &
       YA(X1ind(1),X2ind(2))*Xnncoef(1,2) + YA(X1ind(2),X2ind(2))*Xnncoef(2,2)
  
  return
end subroutine LINTP2_IND_COEF

!=============================================================================
! Addition from V.J. Sep, 1996 - 3 Dimensional linear interpolation
subroutine LINTP3(X1A,X2A,X3A,YA,M,N,L,X1,X2,X3,Y,IER)

  use ModRamMain, ONLY: Real8_
  implicit none
  
  ! Arguments:
  integer, intent(in) :: m, n, l
  integer, intent(out):: ier
  real(kind=Real8_), intent(in), dimension(m)     :: x1a
  real(kind=Real8_), intent(in), dimension(n)     :: x2a
  real(kind=Real8_), intent(in), dimension(l)     :: x3a
  real(kind=Real8_), intent(in), dimension(m,n,l) :: ya
  real(kind=Real8_), intent(in) :: x1, x2, x3
  real(kind=Real8_), intent(out):: y

  ! Local variables:
  integer :: ii, jj, kk
  real(kind=Real8_) :: ytmp(l), yytmp(m,n)
  !-----------------------------------------------------------------------
  IER = 0
  do II=1,M
     do JJ=1,N
        do KK=1,L
           YTMP(KK)=YA(II,JJ,KK)
        end do
        CALL LINTP(X3A,YTMP,L,X3,YYTMP(II,JJ),IER)
     end do
  end do


  CALL LINTP2(X1A,X2A,YYTMP,M,N,X1,X2,Y,IER)
  IER = IER * 10
  RETURN
end subroutine LINTP3

!=============================================================================
function GAMMLN(XX,IER)

  use ModRamMain, ONLY: Real4_, Real8_
  implicit none

  ! Arguments
  real(kind=Real8_) :: gammln
  real(kind=Real8_), intent(in)     :: XX
  integer, intent(out) :: IER

  ! Local vars:
  integer :: j
  real(kind=Real8_) :: x, tmp, ser, &
       half = 0.5, &
       one  = 1.0, &
       fpf  = 5.5, &
       stp  = 2.50662827465E0
  real(kind=Real8_) :: cof(6) = (/76.18009173D0, -86.50532033D0, &
       24.01409822D0, -1.231739516D0, .120858003D-2, -.536382D-5/)
  !-----------------------------------------------------------------------
  IER = 0
  X=XX-ONE
  TMP=X+FPF
  TMP=(X+HALF)*LOG(TMP)-TMP
  SER=ONE
  do J=1,6
     X=X+ONE
     SER=SER+COF(J)/X
  end do

  GAMMLN=TMP+LOG(STP*SER)
  RETURN
end function GAMMLN

!=============================================================================
subroutine GSER(GAMSER,A,X,GLN,IER)
  
  use ModRamMain, ONLY: Real8_
  implicit none
  
  ! Arguments:
  real(kind=Real8_), intent(in) :: a, x
  real(kind=Real8_), intent(out):: gamser, gln
  integer, intent(out) :: IER

  ! Local vars:
  integer, parameter :: ItMax = 100 ! Max iterations.
  integer :: n
  real(kind=Real8_), parameter :: eps = 3.0E-7! A very small number.
  real(kind=Real8_) :: sum, ap, del
  real(kind=Real8_) :: gammln
  !-----------------------------------------------------------------------
  gamser = 0
  ier = 0
  GLN=GAMMLN(A,IER)
  IF(X.LE.0.)THEN
     IF(X.LT.0.) write(*,*) & !PAUSE is bad.
          'WARNING: ARGUMENT X SHOULD NOT BE NEGATIVE IN GSER!'
     RETURN
  ENDIF
  
  AP=A
  SUM=1./A
  DEL=SUM
  do N=1,ITMAX
     AP=AP+1.
     DEL=DEL*X/AP
     SUM=SUM+DEL
     if(abs(DEL) .lt. abs(SUM)*EPS) exit ! Break from loop as necessary.
  end do
 
  ! If max iterations are reached:
  if (n .eq. ItMax) then
     IER = 1
     return
  endif
  ! Else:
  GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
  RETURN
end subroutine GSER

!=============================================================================
subroutine GCF(GAMMCF,A,X,GLN,IER)
  
  use ModRamMain, ONLY: Real8_
  implicit none

  ! Arguments
  real(kind=Real8_), intent(in) :: a, x
  real(kind=Real8_), intent(out):: gln, gammcf
  integer,           intent(out):: ier

  ! Local vars:
  integer,           parameter :: ItMax = 100 ! Max iterations.
  real(kind=Real8_), parameter :: eps = 3.0E-7! A very small number.

  integer :: n
  real(kind=Real8_) :: A0, A1, AN, ANA, ANF, B0, B1, fac, G, gold
  real(kind=Real8_) :: gammln ! external function

  !-----------------------------------------------------------------------
  GLN=GAMMLN(A,IER)
  IER = 0
  
  ! previous value to check for convergence
  GOLD=0.
  ! setting up to evaluate continuous fraction
  A0=1.
  A1=X
  B0=0.
  B1=1.

  ! renormalized factor preventing overflow
  FAC=1.
  do N=1,ITMAX
     AN=FLOAT(N)
     ANA=AN-A
     
     ! one step of the recurrence
     A0=(A1+A0*ANA)*FAC
     B0=(B1+B0*ANA)*FAC
     
     ! next step
     ANF=AN*FAC
     A1=X*A0+ANF*A1
     B1=X*B0+ANF*B1
     
     ! time to renormalize ?
     IF(A1.NE.0.)THEN
        FAC=1./A1
        G=B1*FAC
        
        ! Exit loop if converged
        IF(ABS((G-GOLD)/G).LT.EPS) exit
        GOLD=G
     ENDIF
  end do
     write(*,*) 'THE VALUE OF N IS', n
  ! If error:
  if(n .ge. ItMax) then
     IER = 1
     RETURN
  end if

  ! If converged:
  GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
  RETURN

end subroutine GCF
  
!=============================================================================
subroutine TRIDAG(A,B,C,R,U,N,IER)
  
  use ModRamMain, ONLY: Real8_

  implicit none

  ! Arguments
  integer, intent(in)  :: N
  integer, intent(out) :: IER
  real(kind=Real8_), intent(in),  dimension(N) :: A, B, C, R
  real(kind=Real8_), intent(out), dimension(N) :: U
 
  ! Local vars
  integer, parameter :: nMax = 100
  integer :: j
  real(kind=Real8_), dimension(N) ::  GAM
  real(kind=Real8_) :: BET
  !-----------------------------------------------------------------------
  ! problem can be simplified to N-1
  IF(B(1).EQ.0.)THEN
     IER = 1
     RETURN
  ENDIF
  IER = 0
  BET=B(1)
  U(1)=R(1)/BET
 
  ! decomposition and forward substitution
  do J=2,N
     GAM(J)=C(J-1)/BET
     BET=B(J)-A(J)*GAM(J)
     ! if algotithm fails...
     IF(BET.EQ.0.)THEN
        IER = 2
        RETURN
     ENDIF
     U(J)=(R(J)-A(J)*U(J-1))/BET
  end do

  ! back substitution
  do J=N-1,1,-1
     U(J)=U(J)-GAM(J+1)*U(J+1)
  end do
  RETURN
end subroutine TRIDAG

!=============================================================================
subroutine bessel2(jn,arg1,bs)
  ! This version of bessel2 returns  bs(1),n=jn-1;bs(2),n=jn;bs(3),n=jn+1.

  use ModRamMain, ONLY: Real8_
  implicit none
  ! Arguments:
  integer, intent(in) :: jn
  real(kind=Real8_),    intent(in) :: arg1
  real(kind=Real8_), intent(out) :: bs(3)  ! Can you trust a variable named B.S.?

  ! Local vars:
  integer :: ier, n, jj(3), in, ij, ii, ie
  real(kind=Real8_) :: arg, bb, efac, bessj, bessj0, bessj1
  !-----------------------------------------------------------------------
  arg=arg1
  if(arg.lt.0)arg=-arg
  jj(1)=jn-1
  jj(2)=jn
  jj(3)=jn+1

  do in=1,3
     ij=abs(jj(in))
     select case (ij)
     case(0)
        bs(in)=bessj0(arg)
     case(1)
        bs(in)=bessj1(arg)
     case default
        bs(in)=bessj(ij,arg)
     end select
  end do

  if(arg1.lt.0.)then
     do in=1,3
        ij=jj(in)
        if(ij.lt.0)ij=-ij
        ie=ij
        if(ie.eq.0) exit 
        efac=1.0
        do  ii=1,ie
           efac=efac*(-1)
        end do
        bs(in)=efac*bs(in)
     end do
  endif

  ! check jj for negative values and calculate negative harmonics
  do in=1,3
     ij=jj(in)
     if(ij.ge.0) cycle
     ie=-ij
     efac=1.0
     do  ii=1,ie
        efac=efac*(-1)
     end do
     bs(in)=efac*bs(in)
  end do
  return
end subroutine bessel2

!=============================================================================
function bessj(n,x)

  use ModRamMain, ONLY: Real8_
  implicit none

  ! Arguments
  integer, intent(in) :: n
  real(kind=Real8_), intent(in)    :: x
  real(kind=Real8_) :: bessj!(n,int(x))

  ! Local vars:
  real(kind=Real8_),    parameter :: bigno=1.0e10, bigni=1.0e-10
  integer, parameter :: iacc=40
  integer :: j, m
  real(kind=Real8_) :: tox, bjm, bj, bjp, sum, jsum
  real(kind=Real8_) :: bessj0, bessj1 ! External functions

  !-----------------------------------------------------------------------
  if(n.lt.2) then
     write(*,*) 'bad argument n in bessj-exiting routine '
     return
  end if

  tox=2./x
  if(x.gt.float(n))then
     bjm=bessj0(x)
     bj=bessj1(x)
     do j=1,n-1
        bjp=j*tox*bj-bjm
        bjm=bj
        bj=bjp
     end do
     bessj=bj
  else
     m=2*((n+int(sqrt(float(iacc*n))))/2)
     bessj=0.
     jsum=0
     sum=0.
     bjp=0.
     bj=1.
     do j=m,1,-1
        bjm=j*tox*bj-bjp
        bjp=bj
        bj=bjm
        if(abs(bj).gt.bigno)then
           bj=bj*bigni
           bjp=bjp*bigni
           bessj=bessj*bigni
           sum=sum*bigni
        endif
        if(jsum.ne.0)sum=sum+bj
        jsum=1-jsum
        if(j.eq.n)bessj=bjp
     end do
     sum=2.*sum-bj
     bessj=bessj/sum
  endif
  return
end function bessj
!=============================================================================
function bessj0(x)
  
  use ModRamMain, ONLY: Real8_
  implicit none
  
  ! Arguments:
  real(kind=Real8_), intent(in) :: x
  real(kind=Real8_) :: bessj0!(int(x))
  
  !Local Vars 
  real(kind=Real8_) :: y, z, xx, ax
  real(kind=Real8_), parameter :: &
       p1 = 1.e0,            q1 = -.1562499995e-1,&
       p2 = -.1098628627e-2, q2 = .1430488765e-3,&
       p3 = .2734510407e-4,  q3 = -.6911147651e-5,&
       p4 = -.2073370639e-5, q4 = .7621095161e-6,& 
       p5= .2093887211e-6,   q5 = -.934945152e-7
  real(kind=Real8_), parameter :: &
       r1 = 57568490574.e0,  s1 = 57568490411.e0,&
       r2 = -13362590354.e0, s2 = 1029532985.e0,&
       r3 = 651619640.7e0,   s3 = 9494680.718e0,&
       r4 = -11214424.18e0,  s4 = 59272.64853e0,&
       r5 = 77392.33017e0,   s5 = 267.8532712e0,&
       r6 = -184.9052456e0,  s6 = 1.e0
  !-----------------------------------------------------------------------
  if(abs(x).lt.8.)then
     y=x**2
     bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) / &
          (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
  else
     ax=abs(x)
     z=8./ax
     y=z**2
     xx=ax-.785398164
     bessj0=sqrt(.636619772/ax)*(cos(xx) * &
          (p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx) * &
          (q1+y*(q2+y*(q3+y*(q4+y*q5)))))
  endif
  return
end function bessj0

!=============================================================================
function bessj1(x)
  
  use ModRamMain, ONLY: Real4_, Real8_
  
  implicit none

  real(kind=Real8_), intent(in) :: x
  real(kind=Real8_) :: bessj1!(int(x))
  real(kind=Real8_) :: xx, ax, y, z
  real(kind=Real8_), parameter :: cOne = 1.0
  real(kind=Real8_), parameter :: &
       r1 = 72362614232.E0, s1 = 144725228442.d0, &
       r2 = -7895059235.d0, s2 = 2300535178.d0, &
       r3 = 242396853.1d0,  s3 = 18583304.74d0, &
       r4 = -2972611.439d0, s4 = 99447.43394d0, &
       r5 = 15704.48260d0,  s5 = 376.9991397d0, &
       r6 = -30.16036606d0, s6 = 1.d0
  real(kind=Real8_), parameter :: &
       p1 = 1.d0,            q1 = .04687499995d0, &
       p2 = .183105d-2,      q2 = -.2002690873d-3, &
       p3 = -.3516396496d-4, q3 = .8449199096d-5, &
       p4 = .2457520174d-5,  q4 = -.88228987d-6, &
       p5 = -.240337019d-6,  q5 = .105787412d-6
  !-----------------------------------------------------------------------
  if(abs(x).lt.8.)then
     y=x**2
     bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) / &
          (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
  else
     ax=abs(x)
     z=8./ax
     y=z**2
     xx=ax-2.356194491
     bessj1=sqrt(.636619772/ax)*(cos(xx) * &
          (p1+y*(p2+y*(p3+y*(p4+y*p5)))) - &
          z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))) * &
          sign(cOne,x) ! Why not abs???
  endif
  return
end function bessj1

!=============================================================================
