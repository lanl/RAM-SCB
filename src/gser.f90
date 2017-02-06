	FUNCTION gser_s(a,x,gln)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : gammln
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


	FUNCTION gser_v(a,x,gln)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
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
