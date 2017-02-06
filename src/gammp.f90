	FUNCTION gammp_s(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
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


	FUNCTION gammp_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
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
