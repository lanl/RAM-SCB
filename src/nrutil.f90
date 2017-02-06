MODULE nrutil
        USE nrtype
        IMPLICIT NONE
        INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
        INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
        INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
        INTEGER(I4B), PARAMETER :: NPAR_POLY=8
        INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
        INTERFACE array_copy
                MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
        END INTERFACE
        INTERFACE swap
                MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
                        swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
                        masked_swap_rs,masked_swap_rv,masked_swap_rm
        END INTERFACE
        INTERFACE reallocate
                MODULE PROCEDURE reallocate_rv,reallocate_rm,&
                        reallocate_iv,reallocate_im,reallocate_hv
        END INTERFACE
        INTERFACE imaxloc
                MODULE PROCEDURE imaxloc_r,imaxloc_i
        END INTERFACE
        INTERFACE assert
                MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
        END INTERFACE
        INTERFACE assert_eq
                MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
        END INTERFACE
        INTERFACE arth
                MODULE PROCEDURE arth_r, arth_d, arth_i
        END INTERFACE
        INTERFACE geop
                MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
        END INTERFACE
        INTERFACE poly
                MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
                        poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
        END INTERFACE
        INTERFACE outerprod
                MODULE PROCEDURE outerprod_r,outerprod_d
        END INTERFACE
        INTERFACE outerdiff
                MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
        END INTERFACE
        INTERFACE scatter_add
                MODULE PROCEDURE scatter_add_r,scatter_add_d
        END INTERFACE
        INTERFACE scatter_max
                MODULE PROCEDURE scatter_max_r,scatter_max_d
        END INTERFACE
        INTERFACE diagadd
                MODULE PROCEDURE diagadd_rv,diagadd_r
        END INTERFACE
        INTERFACE diagmult
                MODULE PROCEDURE diagmult_rv,diagmult_r
        END INTERFACE
        INTERFACE get_diag
                MODULE PROCEDURE get_diag_rv, get_diag_dv
        END INTERFACE
        INTERFACE put_diag
                MODULE PROCEDURE put_diag_rv, put_diag_r
        END INTERFACE
CONTAINS
!BL
        SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
        REAL(SP), DIMENSION(:), INTENT(IN) :: src
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied=MIN(SIZE(src),SIZE(dest))
        n_not_copied=SIZE(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        END SUBROUTINE array_copy_r
!BL
        SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
        REAL(DP), DIMENSION(:), INTENT(IN) :: src
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied=MIN(SIZE(src),SIZE(dest))
        n_not_copied=SIZE(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        END SUBROUTINE array_copy_d
!BL
        SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied=MIN(SIZE(src),SIZE(dest))
        n_not_copied=SIZE(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
        END SUBROUTINE array_copy_i
!BL
!BL
        SUBROUTINE swap_i(a,b)
        INTEGER(I4B), INTENT(INOUT) :: a,b
        INTEGER(I4B) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_i
!BL
        SUBROUTINE swap_r(a,b)
        REAL(SP), INTENT(INOUT) :: a,b
        REAL(SP) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_r
!BL
        SUBROUTINE swap_rv(a,b)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(DP), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_rv
!BL
        SUBROUTINE swap_c(a,b)
        COMPLEX(SPC), INTENT(INOUT) :: a,b
        COMPLEX(SPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_c
!BL
        SUBROUTINE swap_cv(a,b)
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cv
!BL
        SUBROUTINE swap_cm(a,b)
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a,1),SIZE(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_cm
!BL
        SUBROUTINE swap_z(a,b)
        COMPLEX(DPC), INTENT(INOUT) :: a,b
        COMPLEX(DPC) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_z
!BL
        SUBROUTINE swap_zv(a,b)
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zv
!BL
        SUBROUTINE swap_zm(a,b)
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a,1),SIZE(a,2)) :: dum
        dum=a
        a=b
        b=dum
        END SUBROUTINE swap_zm
!BL
        SUBROUTINE masked_swap_rs(a,b,mask)
        REAL(SP), INTENT(INOUT) :: a,b
        LOGICAL(LGT), INTENT(IN) :: mask
        REAL(SP) :: swp
        IF (mask) THEN
                swp=a
                a=b
                b=swp
        END IF
        END SUBROUTINE masked_swap_rs
!BL
        SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(SIZE(a)) :: swp
        WHERE (mask)
                swp=a
                a=b
                b=swp
        END WHERE
        END SUBROUTINE masked_swap_rv
!BL
        SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(SIZE(a,1),SIZE(a,2)) :: swp
        WHERE (mask)
                swp=a
                a=b
                b=swp
        END WHERE
        END SUBROUTINE masked_swap_rm
!BL
!BL
        FUNCTION reallocate_rv(p,n)
        REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        ALLOCATE(reallocate_rv(n),stat=ierr)
        IF (ierr /= 0) CALL &
                nrerror('reallocate_rv: problem in attempt to allocate memory')
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold=SIZE(p)
        reallocate_rv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
        END FUNCTION reallocate_rv
!BL
        FUNCTION reallocate_iv(p,n)
        INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        ALLOCATE(reallocate_iv(n),stat=ierr)
        IF (ierr /= 0) CALL &
                nrerror('reallocate_iv: problem in attempt to allocate memory')
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold=SIZE(p)
        reallocate_iv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
        END FUNCTION reallocate_iv
!BL
        FUNCTION reallocate_hv(p,n)
        CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        ALLOCATE(reallocate_hv(n),stat=ierr)
        IF (ierr /= 0) CALL &
                nrerror('reallocate_hv: problem in attempt to allocate memory')
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold=SIZE(p)
        reallocate_hv(1:MIN(nold,n))=p(1:MIN(nold,n))
        DEALLOCATE(p)
        END FUNCTION reallocate_hv
!BL
        FUNCTION reallocate_rm(p,n,m)
        REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        ALLOCATE(reallocate_rm(n,m),stat=ierr)
        IF (ierr /= 0) CALL &
                nrerror('reallocate_rm: problem in attempt to allocate memory')
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold=SIZE(p,1)
        mold=SIZE(p,2)
        reallocate_rm(1:MIN(nold,n),1:MIN(mold,m))=&
                p(1:MIN(nold,n),1:MIN(mold,m))
        DEALLOCATE(p)
        END FUNCTION reallocate_rm
!BL
        FUNCTION reallocate_im(p,n,m)
        INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        ALLOCATE(reallocate_im(n,m),stat=ierr)
        IF (ierr /= 0) CALL &
                nrerror('reallocate_im: problem in attempt to allocate memory')
        IF (.NOT. ASSOCIATED(p)) RETURN
        nold=SIZE(p,1)
        mold=SIZE(p,2)
        reallocate_im(1:MIN(nold,n),1:MIN(mold,m))=&
                p(1:MIN(nold,n),1:MIN(mold,m))
        DEALLOCATE(p)
        END FUNCTION reallocate_im
!BL
        FUNCTION ifirstloc(mask)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        INTEGER(I4B) :: ifirstloc
        INTEGER(I4B), DIMENSION(1) :: loc
        loc=MAXLOC(MERGE(1,0,mask))
        ifirstloc=loc(1)
        IF (.NOT. mask(ifirstloc)) ifirstloc=SIZE(mask)+1
        END FUNCTION ifirstloc
!BL
        FUNCTION imaxloc_r(arr)
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B) :: imaxloc_r
        INTEGER(I4B), DIMENSION(1) :: imax
        imax=MAXLOC(arr(:))
        imaxloc_r=imax(1)
        END FUNCTION imaxloc_r
!BL
        FUNCTION imaxloc_i(iarr)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
        INTEGER(I4B), DIMENSION(1) :: imax
        INTEGER(I4B) :: imaxloc_i
        imax=MAXLOC(iarr(:))
        imaxloc_i=imax(1)
        END FUNCTION imaxloc_i
!BL
        FUNCTION iminloc(arr)
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr   
        ! changed type here from SP to DP to be able to use it 
        INTEGER(I4B), DIMENSION(1) :: imin
        INTEGER(I4B) :: iminloc
        imin=MINLOC(arr(:))
        iminloc=imin(1)
        END FUNCTION iminloc
!BL
        SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        IF (.NOT. n1) THEN
                WRITE (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert1'
        END IF
        END SUBROUTINE assert1
!BL
        SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        IF (.NOT. (n1 .AND. n2)) THEN
                WRITE (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert2'
        END IF
        END SUBROUTINE assert2
!BL
        SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        IF (.NOT. (n1 .AND. n2 .AND. n3)) THEN
                WRITE (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert3'
        END IF
        END SUBROUTINE assert3
!BL
        SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        IF (.NOT. (n1 .AND. n2 .AND. n3 .AND. n4)) THEN
                WRITE (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert4'
        END IF
        END SUBROUTINE assert4
!BL
        SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        IF (.NOT. ALL(n)) THEN
                WRITE (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert_v'
        END IF
        END SUBROUTINE assert_v
!BL
        FUNCTION assert_eq2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2
        INTEGER :: assert_eq2
        IF (n1 == n2) THEN
                assert_eq2=n1
        ELSE
                WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq2'
        END IF
        END FUNCTION assert_eq2
!BL
        FUNCTION assert_eq3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3
        INTEGER :: assert_eq3
        IF (n1 == n2 .AND. n2 == n3) THEN
                assert_eq3=n1
        ELSE
                WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq3'
        END IF
        END FUNCTION assert_eq3
!BL
        FUNCTION assert_eq4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1, n2, n3, n4
        INTEGER :: assert_eq4
        IF (n1 == n2 .AND. n2 == n3 .AND. n3 == n4) THEN
                assert_eq4 = n1
        ELSE
                WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq4'
        END IF
        END FUNCTION assert_eq4
!BL
        FUNCTION assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER :: assert_eqn
        IF (ALL(nn(2:) == nn(1))) THEN
                assert_eqn=nn(1)
        ELSE
                WRITE (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eqn'
        END IF
        END FUNCTION assert_eqn
!BL
        SUBROUTINE nrerror(string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        WRITE (*,*) 'nrerror: ',string
        STOP 'program terminated by nrerror'
        END SUBROUTINE nrerror
!BL
        FUNCTION arth_r(first,increment,n)
        REAL(SP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: arth_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        IF (n > 0) arth_r(1)=first
        IF (n <= NPAR_ARTH) THEN
                DO k=2,n
                        arth_r(k)=arth_r(k-1)+increment
                END DO
        ELSE
                DO k=2,NPAR2_ARTH
                        arth_r(k)=arth_r(k-1)+increment
                END DO
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        arth_r(k+1:MIN(k2,n))=temp+arth_r(1:MIN(k,n-k))
                        temp=temp+temp
                        k=k2
                END DO
        END IF
        END FUNCTION arth_r
!BL
        FUNCTION arth_d(first,increment,n)
        REAL(DP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: arth_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        IF (n > 0) arth_d(1)=first
        IF (n <= NPAR_ARTH) THEN
                DO k=2,n
                        arth_d(k)=arth_d(k-1)+increment
                END DO
        ELSE
                DO k=2,NPAR2_ARTH
                        arth_d(k)=arth_d(k-1)+increment
                END DO
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        arth_d(k+1:MIN(k2,n))=temp+arth_d(1:MIN(k,n-k))
                        temp=temp+temp
                        k=k2
                END DO
        END IF
        END FUNCTION arth_d
!BL
        FUNCTION arth_i(first,increment,n)
        INTEGER(I4B), INTENT(IN) :: first,increment,n
        INTEGER(I4B), DIMENSION(n) :: arth_i
        INTEGER(I4B) :: k,k2,temp
        IF (n > 0) arth_i(1)=first
        IF (n <= NPAR_ARTH) THEN
                DO k=2,n
                        arth_i(k)=arth_i(k-1)+increment
                END DO
        ELSE
                DO k=2,NPAR2_ARTH
                        arth_i(k)=arth_i(k-1)+increment
                END DO
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        arth_i(k+1:MIN(k2,n))=temp+arth_i(1:MIN(k,n-k))
                        temp=temp+temp
                        k=k2
                END DO
        END IF
        END FUNCTION arth_i
!BL
!BL
        FUNCTION geop_r(first,factor,n)
        REAL(SP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: geop_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        IF (n > 0) geop_r(1)=first
        IF (n <= NPAR_GEOP) THEN
                DO k=2,n
                        geop_r(k)=geop_r(k-1)*factor
                END DO
        ELSE
                DO k=2,NPAR2_GEOP
                        geop_r(k)=geop_r(k-1)*factor
                END DO
                temp=factor**NPAR2_GEOP
                k=NPAR2_GEOP
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        geop_r(k+1:MIN(k2,n))=temp*geop_r(1:MIN(k,n-k))
                        temp=temp*temp
                        k=k2
                END DO
        END IF
        END FUNCTION geop_r
!BL
        FUNCTION geop_d(first,factor,n)
        REAL(DP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: geop_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        IF (n > 0) geop_d(1)=first
        IF (n <= NPAR_GEOP) THEN
                DO k=2,n
                        geop_d(k)=geop_d(k-1)*factor
                END DO
        ELSE
                DO k=2,NPAR2_GEOP
                        geop_d(k)=geop_d(k-1)*factor
                END DO
                temp=factor**NPAR2_GEOP
                k=NPAR2_GEOP
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        geop_d(k+1:MIN(k2,n))=temp*geop_d(1:MIN(k,n-k))
                        temp=temp*temp
                        k=k2
                END DO
        END IF
        END FUNCTION geop_d
!BL
        FUNCTION geop_i(first,factor,n)
        INTEGER(I4B), INTENT(IN) :: first,factor,n
        INTEGER(I4B), DIMENSION(n) :: geop_i
        INTEGER(I4B) :: k,k2,temp
        IF (n > 0) geop_i(1)=first
        IF (n <= NPAR_GEOP) THEN
                DO k=2,n
                        geop_i(k)=geop_i(k-1)*factor
                END DO
        ELSE
                DO k=2,NPAR2_GEOP
                        geop_i(k)=geop_i(k-1)*factor
                END DO
                temp=factor**NPAR2_GEOP
                k=NPAR2_GEOP
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        geop_i(k+1:MIN(k2,n))=temp*geop_i(1:MIN(k,n-k))
                        temp=temp*temp
                        k=k2
                END DO
        END IF
        END FUNCTION geop_i
!BL
        FUNCTION geop_c(first,factor,n)
        COMPLEX(SP), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        COMPLEX(SP), DIMENSION(n) :: geop_c
        INTEGER(I4B) :: k,k2
        COMPLEX(SP) :: temp
        IF (n > 0) geop_c(1)=first
        IF (n <= NPAR_GEOP) THEN
                DO k=2,n
                        geop_c(k)=geop_c(k-1)*factor
                END DO
        ELSE
                DO k=2,NPAR2_GEOP
                        geop_c(k)=geop_c(k-1)*factor
                END DO
                temp=factor**NPAR2_GEOP
                k=NPAR2_GEOP
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        geop_c(k+1:MIN(k2,n))=temp*geop_c(1:MIN(k,n-k))
                        temp=temp*temp
                        k=k2
                END DO
        END IF
        END FUNCTION geop_c
!BL
        FUNCTION geop_dv(first,factor,n)
        REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(SIZE(first),n) :: geop_dv
        INTEGER(I4B) :: k,k2
        REAL(DP), DIMENSION(SIZE(first)) :: temp
        IF (n > 0) geop_dv(:,1)=first(:)
        IF (n <= NPAR_GEOP) THEN
                DO k=2,n
                        geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
                END DO
        ELSE
                DO k=2,NPAR2_GEOP
                        geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
                END DO
                temp=factor**NPAR2_GEOP
                k=NPAR2_GEOP
                DO
                        IF (k >= n) EXIT
                        k2=k+k
                        geop_dv(:,k+1:MIN(k2,n))=geop_dv(:,1:MIN(k,n-k))*&
                                SPREAD(temp,2,SIZE(geop_dv(:,1:MIN(k,n-k)),2))
                        temp=temp*temp
                        k=k2
                END DO
        END IF
        END FUNCTION geop_dv
!BL
!BL
        FUNCTION poly_rr(x,coeffs)
        REAL(SP), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(SP) :: poly_rr
        REAL(SP) :: pow
        REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=SIZE(coeffs)
        IF (n <= 0) THEN
                poly_rr=0.0_sp
        ELSE IF (n < NPAR_POLY) THEN
                poly_rr=coeffs(n)
                DO i=n-1,1,-1
                        poly_rr=x*poly_rr+coeffs(i)
                END DO
        ELSE
                ALLOCATE(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                DO
                        vec(n+1)=0.0_sp
                        nn=ISHFT(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        IF (nn == 1) EXIT
                        pow=pow*pow
                        n=nn
                END DO
                poly_rr=vec(1)
                DEALLOCATE(vec)
        END IF
        END FUNCTION poly_rr
!BL
        FUNCTION poly_dd(x,coeffs)
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(DP) :: poly_dd
        REAL(DP) :: pow
        REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=SIZE(coeffs)
        IF (n <= 0) THEN
                poly_dd=0.0_dp
        ELSE IF (n < NPAR_POLY) THEN
                poly_dd=coeffs(n)
                DO i=n-1,1,-1
                        poly_dd=x*poly_dd+coeffs(i)
                END DO
        ELSE
                ALLOCATE(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                DO
                        vec(n+1)=0.0_dp
                        nn=ISHFT(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        IF (nn == 1) EXIT
                        pow=pow*pow
                        n=nn
                END DO
                poly_dd=vec(1)
                DEALLOCATE(vec)
        END IF
        END FUNCTION poly_dd
!BL
        FUNCTION poly_rc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_rc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=SIZE(coeffs)
        IF (n <= 0) THEN
                poly_rc=0.0_sp
        ELSE IF (n < NPAR_POLY) THEN
                poly_rc=coeffs(n)
                DO i=n-1,1,-1
                        poly_rc=x*poly_rc+coeffs(i)
                END DO
        ELSE
                ALLOCATE(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                DO
                        vec(n+1)=0.0_sp
                        nn=ISHFT(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        IF (nn == 1) EXIT
                        pow=pow*pow
                        n=nn
                END DO
                poly_rc=vec(1)
                DEALLOCATE(vec)
        END IF
        END FUNCTION poly_rc
!BL
        FUNCTION poly_cc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_cc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=SIZE(coeffs)
        IF (n <= 0) THEN
                poly_cc=0.0_sp
        ELSE IF (n < NPAR_POLY) THEN
                poly_cc=coeffs(n)
                DO i=n-1,1,-1
                        poly_cc=x*poly_cc+coeffs(i)
                END DO
        ELSE
                ALLOCATE(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                DO
                        vec(n+1)=0.0_sp
                        nn=ISHFT(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        IF (nn == 1) EXIT
                        pow=pow*pow
                        n=nn
                END DO
                poly_cc=vec(1)
                DEALLOCATE(vec)
        END IF
        END FUNCTION poly_cc
!BL
        FUNCTION poly_rrv(x,coeffs)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(SP), DIMENSION(SIZE(x)) :: poly_rrv
        INTEGER(I4B) :: i,n,m
        m=SIZE(coeffs)
        n=SIZE(x)
        IF (m <= 0) THEN
                poly_rrv=0.0_sp
        ELSE IF (m < n .OR. m < NPAR_POLY) THEN
                poly_rrv=coeffs(m)
                DO i=m-1,1,-1
                        poly_rrv=x*poly_rrv+coeffs(i)
                END DO
        ELSE
                DO i=1,n
                        poly_rrv(i)=poly_rr(x(i),coeffs)
                END DO
        END IF
        END FUNCTION poly_rrv
!BL
        FUNCTION poly_ddv(x,coeffs)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(DP), DIMENSION(SIZE(x)) :: poly_ddv
        INTEGER(I4B) :: i,n,m
        m=SIZE(coeffs)
        n=SIZE(x)
        IF (m <= 0) THEN
                poly_ddv=0.0_dp
        ELSE IF (m < n .OR. m < NPAR_POLY) THEN
                poly_ddv=coeffs(m)
                DO i=m-1,1,-1
                        poly_ddv=x*poly_ddv+coeffs(i)
                END DO
        ELSE
                DO i=1,n
                        poly_ddv(i)=poly_dd(x(i),coeffs)
                END DO
        END IF
        END FUNCTION poly_ddv
!BL
        FUNCTION poly_msk_rrv(x,coeffs,mask)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(SIZE(x)) :: poly_msk_rrv
        poly_msk_rrv=UNPACK(poly_rrv(PACK(x,mask),coeffs),mask,0.0_sp)
        END FUNCTION poly_msk_rrv
!BL
        FUNCTION poly_msk_ddv(x,coeffs,mask)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(DP), DIMENSION(SIZE(x)) :: poly_msk_ddv
        poly_msk_ddv=UNPACK(poly_ddv(PACK(x,mask),coeffs),mask,0.0_dp)
        END FUNCTION poly_msk_ddv
!BL
!BL
        FUNCTION zroots_unity(n,nn)
        INTEGER(I4B), INTENT(IN) :: n,nn
        COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
        INTEGER(I4B) :: k
        REAL(SP) :: theta
        zroots_unity(1)=1.0
        theta=TWOPI/n
        k=1
        DO
                IF (k >= nn) EXIT
                zroots_unity(k+1)=CMPLX(COS(k*theta),SIN(k*theta),SPC)
                zroots_unity(k+2:MIN(2*k,nn))=zroots_unity(k+1)*&
                        zroots_unity(2:MIN(k,nn-k))
                k=2*k
        END DO
        END FUNCTION zroots_unity
!BL
        FUNCTION outerprod_r(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_r
        outerprod_r = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerprod_r
!BL
        FUNCTION outerprod_d(a,b)
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_d
        outerprod_d = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerprod_d
!BL
        FUNCTION outerdiv(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerdiv
        outerdiv = SPREAD(a,dim=2,ncopies=SIZE(b)) / &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerdiv
!BL
        FUNCTION outersum(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outersum
        outersum = SPREAD(a,dim=2,ncopies=SIZE(b)) + &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outersum
!BL
        FUNCTION outerdiff_r(a,b)
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerdiff_r
        outerdiff_r = SPREAD(a,dim=2,ncopies=SIZE(b)) - &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerdiff_r
!BL
        FUNCTION outerdiff_d(a,b)
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(SIZE(a),SIZE(b)) :: outerdiff_d
        outerdiff_d = SPREAD(a,dim=2,ncopies=SIZE(b)) - &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerdiff_d
!BL
        FUNCTION outerdiff_i(a,b)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
        INTEGER(I4B), DIMENSION(SIZE(a),SIZE(b)) :: outerdiff_i
        outerdiff_i = SPREAD(a,dim=2,ncopies=SIZE(b)) - &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerdiff_i
!BL
        FUNCTION outerand(a,b)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
        LOGICAL(LGT), DIMENSION(SIZE(a),SIZE(b)) :: outerand
        outerand = SPREAD(a,dim=2,ncopies=SIZE(b)) .AND. &
                SPREAD(b,dim=1,ncopies=SIZE(a))
        END FUNCTION outerand
!BL
        SUBROUTINE scatter_add_r(dest,source,dest_index)
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(SP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_add_r')
        m=SIZE(dest)
        DO j=1,n
                i=dest_index(j)
                IF (i > 0 .AND. i <= m) dest(i)=dest(i)+source(j)
        END DO
        END SUBROUTINE scatter_add_r
        SUBROUTINE scatter_add_d(dest,source,dest_index)
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_add_d')
        m=SIZE(dest)
        DO j=1,n
                i=dest_index(j)
                IF (i > 0 .AND. i <= m) dest(i)=dest(i)+source(j)
        END DO
        END SUBROUTINE scatter_add_d
        SUBROUTINE scatter_max_r(dest,source,dest_index)
        REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(SP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_max_r')
        m=SIZE(dest)
        DO j=1,n
                i=dest_index(j)
                IF (i > 0 .AND. i <= m) dest(i)=MAX(dest(i),source(j))
        END DO
        END SUBROUTINE scatter_max_r
        SUBROUTINE scatter_max_d(dest,source,dest_index)
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_max_d')
        m=SIZE(dest)
        DO j=1,n
                i=dest_index(j)
                IF (i > 0 .AND. i <= m) dest(i)=MAX(dest(i),source(j))
        END DO
        END SUBROUTINE scatter_max_d
!BL
        SUBROUTINE diagadd_rv(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), DIMENSION(:), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = assert_eq2(SIZE(diag),MIN(SIZE(mat,1),SIZE(mat,2)),'diagadd_rv')
        DO j=1,n
                mat(j,j)=mat(j,j)+diag(j)
        END DO
        END SUBROUTINE diagadd_rv
!BL
        SUBROUTINE diagadd_r(mat,diag)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(SP), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = MIN(SIZE(mat,1),SIZE(mat,2))
        DO j=1,n
                mat(j,j)=mat(j,j)+diag
        END DO
        END SUBROUTINE diagadd_r
!BL
        SUBROUTINE diagmult_rv(mat,diag)   ! changed types to DP
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), DIMENSION(:), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = assert_eq2(SIZE(diag),MIN(SIZE(mat,1),SIZE(mat,2)),'diagmult_rv')
        DO j=1,n
                mat(j,j)=mat(j,j)*diag(j)
        END DO
        END SUBROUTINE diagmult_rv
!BL
        SUBROUTINE diagmult_r(mat,diag)    ! changed types to DP
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), INTENT(IN) :: diag
        INTEGER(I4B) :: j,n
        n = MIN(SIZE(mat,1),SIZE(mat,2))
        DO j=1,n
                mat(j,j)=mat(j,j)*diag
        END DO
        END SUBROUTINE diagmult_r
!BL
        FUNCTION get_diag_rv(mat)
        REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(SP), DIMENSION(SIZE(mat,1)) :: get_diag_rv
        INTEGER(I4B) :: j
        j=assert_eq2(SIZE(mat,1),SIZE(mat,2),'get_diag_rv')
        DO j=1,SIZE(mat,1)
                get_diag_rv(j)=mat(j,j)
        END DO
        END FUNCTION get_diag_rv
!BL
        FUNCTION get_diag_dv(mat)
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(DP), DIMENSION(SIZE(mat,1)) :: get_diag_dv
        INTEGER(I4B) :: j
        j=assert_eq2(SIZE(mat,1),SIZE(mat,2),'get_diag_dv')
        DO j=1,SIZE(mat,1)
                get_diag_dv(j)=mat(j,j)
        END DO
        END FUNCTION get_diag_dv
!BL
        SUBROUTINE put_diag_rv(diagv,mat)
        REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        INTEGER(I4B) :: j,n
        n=assert_eq2(SIZE(diagv),MIN(SIZE(mat,1),SIZE(mat,2)),'put_diag_rv')
        DO j=1,n
                mat(j,j)=diagv(j)
        END DO
        END SUBROUTINE put_diag_rv
!BL
        SUBROUTINE put_diag_r(scal,mat)
        REAL(SP), INTENT(IN) :: scal
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
        INTEGER(I4B) :: j,n
        n = MIN(SIZE(mat,1),SIZE(mat,2))
        DO j=1,n
                mat(j,j)=scal
        END DO
        END SUBROUTINE put_diag_r
!BL
        SUBROUTINE unit_matrix(mat)
        REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
        INTEGER(I4B) :: i,n
        n=MIN(SIZE(mat,1),SIZE(mat,2))
        mat(:,:)=0.0_sp
        DO i=1,n
                mat(i,i)=1.0_sp
        END DO
        END SUBROUTINE unit_matrix
!BL
        FUNCTION upper_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
        INTEGER(I4B) :: n
        n=0
        IF (PRESENT(extra)) n=extra
        upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
        END FUNCTION upper_triangle
!BL
        FUNCTION lower_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
        INTEGER(I4B) :: n
        n=0
        IF (PRESENT(extra)) n=extra
        lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
        END FUNCTION lower_triangle
!BL
        FUNCTION vabs(v)  ! Changed type to DP
        REAL(DP), DIMENSION(:), INTENT(IN) :: v
        REAL(DP) :: vabs
        vabs = SQRT(DOT_PRODUCT(v,v))
        END FUNCTION vabs
!BL
END MODULE nrutil