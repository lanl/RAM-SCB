!******************************************************************************
! Copyright notice
! This program was prepared by Los Alamos National Security, LLC at Los Alamos 
! National Laboratory (LANL) under contract No. DE-AC52-06NA25396 with the U.S. 
! Department of Energy (DOE).
! Permission is granted to the public to copy and use this software without 
! charge, provided that this Notice and any statement of authorship are 
! reproduced on all copies.
! Neither the U.S. Government nor LANS makes any warranty, express or implied, 
! or assumes any liability or responsibility for the use of this software.
! Author: Enrico Camporeale, June 2013
! Converted to f90, summer 2015, by Chris Jeffery
! in't Hout paper (2009) Eq 1.3
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!******************************************************************************

subroutine diff_2d(f,J1,J2,D1,D2,D3,dx,dy,dt,Nx,Ny,Titer)

  use ModRamMain, ONLY: Real8_
  implicit none

  integer, intent(in) :: Nx,Ny,Titer
  real(kind=Real8_), intent(in) :: dx,dy,dt
  real(kind=Real8_), intent(in), dimension(Nx) :: J1
  real(kind=Real8_), intent(in), dimension(Ny) :: J2
  real(kind=Real8_), intent(inout), dimension(Ny,Nx) :: f,D1,D2,D3

  integer :: i,j,t
  real(kind=Real8_) :: bet,cxxfac,cxyfac,cyyfac,ryyfac,rxyfac,rxyfac2

  real(kind=Real8_), dimension(Nx) :: gamx,J1summ1,J1sump1
  real(kind=Real8_), dimension(Ny) :: gamy,J2summ1,J2sump1
  real(kind=Real8_), dimension(Ny,Nx) :: ax,bx,cx,D1sumip1,D1sumim1,J1D3,J2D3,D2sumjp1,D2sumjm1
  real(kind=Real8_), dimension(Nx,Ny) :: ay,by,cy
  real(kind=Real8_), dimension(Ny,Nx) :: f_1, f_2, f_3, r, r1, r2

!	double bx[Ny][Nx],ax[Ny][Nx], cx[Ny][Nx]
!	double by[Nx][Ny],ay[Nx][Ny], cy[Nx][Ny]
	
  do j=1,Ny,1
     cx(j,1)=0
     bx(j,1)=1
     ax(j,1)=0
     ax(j,Nx)=0
     bx(j,Nx)=1
     cx(j,Nx)=0
  end do

  do i=1,Nx,1
     cy(i,1)=0 
     by(i,1)=1 
     ay(i,1)=0
     ay(i,Ny)=0
     by(i,Ny)=1
     cy(i,Ny)=0
  end do

  do i=2,Nx-1,1
     J1summ1(i) = J1(i)+J1(i-1)
     J1sump1(i) = J1(i)+J1(i+1)
  end do

  do i=2,Ny-1,1
     J2summ1(i) = J2(i)+J2(i-1)
     J2sump1(i) = J2(i)+J2(i+1)
  end do

  r(:,:) = 0
  f_1(:,:) = 0

  cxxfac = dt/dx/dx/8.
  cxyfac = dt/dx/dy/8.
  cyyfac = dt/dy/dy/8.
  ryyfac = dt/dy/dy/4.
  rxyfac = dt/dx/dy/4.
  rxyfac2 = dt/dx/dy/2.

! initialize a,b,c
  do j=2,Ny-1,1
     
     do i=2,Nx-1,1
        
        D1sumip1(j,i) = D1(j,i)+D1(j,i+1)
        D1sumim1(j,i) = D1(j,i)+D1(j,i-1)
        J2D3(j,i) = J2(j+1)*D3(j+1,i) - J2(j-1)*D3(j-1,i)
        J1D3(j,i) = J1(i+1)*D3(j,i+1) - J1(i-1)*D3(j,i-1)

        cx(j,i) = -cxxfac/J1(i)*D1sumip1(j,i)*J1sump1(i)
        cx(j,i) = cx(j,i) - cxyfac/J2(j)*J2D3(j,i)
        
        bx(j,i)= 1 + cxxfac/J1(i)*(D1sumip1(j,i)*J1sump1(i) + D1sumim1(j,i)*J1summ1(i))
        ax(j,i) = - cxxfac/J1(i)*D1sumim1(j,i)*J1summ1(i)
        ax(j,i) = ax(j,i) - cxyfac/J2(j)*(-J2D3(j,i))

        D2sumjm1(j,i) = D2(j,i)+D2(j-1,i)
        D2sumjp1(j,i) = D2(j,i)+D2(j+1,i)

        ay(i,j) = - cyyfac/J2(j)*D2sumjm1(j,i)*J2summ1(j) 
        ay(i,j) = ay(i,j) - cxyfac/J1(i)*(-J1D3(j,i))
			
        by(i,j)= 1 + cyyfac/J2(j)*(D2sumjp1(j,i)*J2sump1(j) + D2sumjm1(j,i)*J2summ1(j))
			
        cy(i,j) = -cyyfac/J2(j)*D2sumjp1(j,i)*J2sump1(j)
        cy(i,j) = cy(i,j) - cxyfac/J1(i)*(J1D3(j,i))
			
     end do
  end do

  do t=1,Titer,1 ! loop in time

! first tridiagonal system: f1	
     do j=2,Ny-1,1

        do i=2,Nx-1,1
						
           r1(j,i) =f(j,i)*(1 - cxxfac/J1(i)*(D1sumip1(j,i)*J1sump1(i) + &
                D1sumim1(j,i)*J1summ1(i)) - ryyfac/J2(j)*(D2sumjp1(j,i)*J2sump1(j) + &
                D2sumjm1(j,i)*J2summ1(j)) )
! Lx
           r1(j,i) = r1(j,i) + cxxfac*f(j,i+1)/J1(i)*D1sumip1(j,i)*J1sump1(i) 
           r1(j,i) = r1(j,i) + cxxfac*f(j,i-1)/J1(i)*D1sumim1(j,i)*J1summ1(i) 
           r1(j,i) = r1(j,i) + cxyfac*f(j,i+1)/J2(j)*J2D3(j,i)
           r1(j,i) = r1(j,i) + cxyfac*f(j,i-1)/J2(j)*(-J2D3(j,i))
! Ly
           r1(j,i) = r1(j,i) + ryyfac*f(j+1,i)/J2(j)*D2sumjp1(j,i)*J2sump1(j)  ! superdiag f(i,j+1)
           r1(j,i) = r1(j,i) + ryyfac*f(j-1,i)/J2(j)*D2sumjm1(j,i)*J2summ1(j)	!  subdiag  f(i,j-1)
           r1(j,i) = r1(j,i) + rxyfac*f(j+1,i)/J1(i)*J1D3(j,i)
           r1(j,i) = r1(j,i) + rxyfac*f(j-1,i)/J1(i)*(-J1D3(j,i))
! Lxy			
           r(j,i) = r1(j,i) + rxyfac2*f(j-1,i-1)*D3(j,i) ! f(i-1,j-1)
           r(j,i) = r(j,i) - rxyfac2*f(j+1,i-1)*D3(j,i)  !f(i-1,j+1)
           r(j,i) = r(j,i) - rxyfac2*f(j-1,i+1)*D3(j,i)   ! f(i+1,j-1)
           r(j,i) = r(j,i) + rxyfac2*f(j+1,i+1)*D3(j,i)   ! f(i+1,j+1)
			
        end do

        ! Outer BC
        r(j,1) = f(j,1)		
        f_1(j,1) = f(j,1)
        f_2(j,1) = f(j,1)
        f_2(j,Nx) = f(j,Nx)
        bet=1.0

        do i=2,Nx,1

           gamx(i) = cx(j,i-1)/bet
           bet=bx(j,i)-ax(j,i)*gamx(i)
!	if (bet== 0.0) {printf("Error in tridiagonal solver !\n")return}
		
           f_1(j,i)=(r(j,i)-ax(j,i)*f_1(j,i-1))/bet

        end do

! Outer BC
        i=Nx
        f_1(j,i)=f(j,i)

        do i=Nx-1,1,-1
           f_1(j,i) = f_1(j,i) - gamx(i+1)*f_1(j,i+1)
        end do

     end do
		
     ! second tridiagonal system: f2
     do i=2,Nx-1,1
        do j=2,Ny-1,1			
! -0.5*Ly					
           r2(j,i)  =  f(j,i)*cyyfac/J2(j)*(D2sumjp1(j,i)*J2sump1(j) + &
                D2sumjm1(j,i)*J2summ1(j))  !diagonal f(i,j)
              
           r2(j,i) = r2(j,i) - cyyfac*f(j+1,i)/J2(j)*D2sumjp1(j,i)*J2sump1(j)  ! superdiag f(i,j+1)
           r2(j,i) = r2(j,i) - cyyfac*f(j-1,i)/J2(j)*D2sumjm1(j,i)*J2summ1(j)	!  subdiag  f(i,j-1)
           r2(j,i) = r2(j,i) - cxyfac*f(j+1,i)/J1(i)*J1D3(j,i)
           r2(j,i) = r2(j,i) - cxyfac*f(j-1,i)/J1(i)*(-J1D3(j,i))
           
           r(j,i) = r2(j,i)+f_1(j,i)
				
        end do

        ! Outer BC
        r(1,i) = f(1,i)
        r(Ny,i)=f(Ny,i)
        r2(Ny,i)=f(Ny,i)
        f_2(1,i) = f(1,i)
        f_3(Ny,i) = f(Ny,i)
        f_3(1,i) = f(1,i)
        bet=1.0

        do j=2,Ny,1
           gamy(j) = cy(i,j-1)/bet
           bet=by(i,j)-ay(i,j)*gamy(j)
	! if (bet== 0.0) {printf("Error in tridiagonal solver !\n")return}
              
           f_2(j,i)=(r(j,i)-ay(i,j)*f_2(j-1,i))/bet

        end do

        do j=Ny-1,1,-1
           f_2(j,i) = f_2(j,i) - gamy(j+1)*f_2(j+1,i)
        end do

     end do

     ! third tridiagonal: f3	

     do j=2,Ny-1,1
        do i=2,Nx-1,1
			
              ! Lxy			
           r(j,i) = r1(j,i) + rxyfac*(f(j-1,i-1)+f_2(j-1,i-1))*D3(j,i) ! f(i-1,j-1)
           r(j,i) = r(j,i) - rxyfac*(f(j+1,i-1)+f_2(j+1,i-1))*D3(j,i)  !f(i-1,j+1)
           r(j,i) = r(j,i) - rxyfac*(f(j-1,i+1)+f_2(j-1,i+1))*D3(j,i)   ! f(i+1,j-1)
           r(j,i) = r(j,i) + rxyfac*(f(j+1,i+1)+f_2(j+1,i+1))*D3(j,i)   ! f(i+1,j+1)

        end do
			
        f_3(j,1) = f(j,1)
        bet=1.0

        do i=2,Nx,1
           gamx(i) = cx(j,i-1)/bet
           bet=bx(j,i)-ax(j,i)*gamx(i)
!	if (bet== 0.0) {printf("Error in tridiagonal solver !\n")return}
			
           f_3(j,i)=(r(j,i)-ax(j,i)*f_3(j,i-1))/bet

        end do

! Outer BC
        f_3(j,Nx)=f(j,Nx)

        do i=Nx-1,1,-1
           f_3(j,i) = f_3(j,i) - gamx(i+1)*f_3(j,i+1)
        end do
        
     end do

     ! fourth tridiagonal system: f2
     do i=2,Nx-1,1
!		r2(Ny*i) = 0
        bet=1.0!by(0)

        do j=2,Ny,1
           gamy(j) = cy(i,j-1)/bet
           bet=by(i,j)-ay(i,j)*gamy(j)
           f(j,i)=((r2(j,i)+f_3(j,i))-ay(i,j)*f(j-1,i))/bet
           
        end do

        do j=Ny-1,1,-1
           f(j,i) = f(j,i) - gamy(j+1)*f(j+1,i)
        end do
        
     end do
      	
     ! apply bc	
     i=1
     do j=1,Ny,1
        f(j,i) = f(j,i+1)
     end do
     
     i=Nx
     do j=1,Ny,1
        f(j,i) = f(j,i-1)
     end do
	
     j=1
     do i=1,Nx,1
        f(j,i) = f(j+1,i)
     end do
	
     j=Ny
     do i=1,Nx,1
        f(j,i) = f(j-1,i)
     end do

  end do

  return
end subroutine diff_2d
