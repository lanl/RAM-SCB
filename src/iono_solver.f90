!^CFG COPYRIGHT UM

subroutine ionosphere_solver(iBlock, Jr, &
     SigmaThTh, SigmaThPs, SigmaPsPs, &
     dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
     dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
     Theta, Psi, dTheta, dPsi,Phi_C)

  !
  ! This subroutine solves for the ionospheric potential PHI
  ! using the field aligned currents Jr and the conductivity tensor Sigma
  ! and its various derivatives as input.
  !
  ! The idea is that the Jr current is balanced by the divergence of the
  ! horizontal currents. The horizontal currents are Sigma . E, 
  ! the height integrated conductivity Sigma is 2x2 antisymmetric matrix
  ! acting on the Theta and Phi components of the electric field.
  ! The electric field is assumed to be a potential field: E = Grad Phi.
  !
  ! We are solving 
  !
  !    Div( Sigma . Grad Phi ) = - Jr
  ! 
  ! in spherical coordinates. 
  !
  ! This leads to a penta-diagonal linear equation:
  !  
  !  C_A*Phi(i,j) + C_B*Phi(i-1,j) + C_C*Phi(i+1,j) +
  !                 C_D*Phi(i,j-1) + C_E*Phi(i,j+1) = RHS
  !
  ! To avoid division by zero at the poles, the equation is 
  ! multiplied by (sin(Theta)*Radius)**2, thus the 
  ! RHS = Jr * (sin(Theta)*Radius)**2
  !
  ! The linear elliptic PDE for the
  ! electric field potential PHI is defined on the domain 
  ! 0 < Theta < ThetaMax  for the northern hemisphere, and
  ! ThetaMin < Theta < PI for the southern hemisphere), and
  ! 0 < Psi < 2 PI.
  !
  ! The following boundary conditions are applied:
  !
  !      PHI(ThetaMax,Psi) = PHI(ThetaMin,Psi) = 0,
  !
  !      PHI(Theta,0) = PHI(Theta,2*PI).
  ! 
  ! There is no boundary at the poles, but to avoid numerical 
  ! difficulties, the cell value at the pole is replaced with the 
  ! average of the first neighbors:
  !
  !      PHI(0,Psi)  = average( PHI(dTheta, Psi) )
  !
  !      PHI(PI,Psi) = average( PHI(PI-dTheta, Psi) )
  !
  ! where dTheta is the grid resolution at the poles.
  !/

  use ModIonosphere
  use IE_ModMain,      ONLY: DoCoupleUaCurrent, LatBoundaryGm, LatBoundary, &
                             NameSolver, UsePreconditioner, UseInitialGuess, &
                             Tolerance, MaxIteration, time_array
  use IE_ModIo,        ONLY: write_prefix, iUnitOut
  use ModIoUnit,       ONLY: UnitTmp_
  use ModLinearSolver, ONLY: gmres, bicgstab, prehepta, Uhepta, Lhepta
  use ModRamMain,      ONLY: PathSceOut
  implicit none

  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1

  ! INPUT ARGUMENTS:
  integer, intent(in) :: iBlock
  real(real8_), intent(in), dimension(nTheta,nPsi) ::  &
       Jr, &
       SigmaThTh, SigmaThPs, SigmaPsPs, &
       dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
       dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
       Theta, Psi
  real(real8_), intent(in) :: dTheta(nTheta), dPsi(nPsi)
  !  real(real8_), intent(in) :: LatBoundaryIM(nPsi)

  ! OUTPUT ARGUMENTS:
  real(real8_), intent(out) :: Phi_C(nTheta,nPsi) ! Potential

  ! Local variables
  real(real8_) :: dTheta2(nTheta), dPsi2(nPsi)
  real(real8_) :: SinTheta_I(nTheta), CosTheta_I(nTheta)
  real(real8_) :: sn, cs, sn2
  real(real8_) :: TermTheta2, TermTheta1, TermPsi2, TermPsi1

  real(real8_) :: PhiOld_CB(nTheta, nPsi, 2) = 0.0  ! Previous solution
  real(real8_) :: lat_weimer(0:nTheta+1,nPsi), mlt_weimer(0:ntheta+1,nPsi), mlt_tmp(0:ntheta+1,nPsi)
  real(real8_), dimension(nTheta, nPsi) :: Jr_sol, PhiOld_CB_sol, SigmaThTh_sol, &
       SigmaPsPs_sol, dSigmaThTh_dTheta_sol, dSigmaThPs_dPsi_sol, dSigmaThPs_dTheta_sol, &
       dSigmaPsPs_dPsi_sol, Theta_sol, Psi_sol, Phi_sol
  integer :: i, j, k, iI, nTotal, iMin, iMax
  ! 1D vectors for the linear solver
  real(real8_), dimension(:), allocatable :: x, y, rhs, b, Bnd_I
  real(real8_) :: dPsi_sol,  dTheta_sol(nPsi), dtheta_tmp
  real(real8_) :: PhiMin, PhiMax

  logical :: DoTest, DoTestMe

  real(real8_)    :: Residual
  integer :: nIteration, iError
  character(len=100) :: NameFile
  external :: matvec_ionosphere

  SAVE
  character(len=*), parameter:: NameSub = 'ionosphere_solver'
  !-------------------------------------------------------------------------

  if(DoTest)write(*,*)'iono_solve starting'

  ! Count the points above the latitude boundary
  nThetaUsed = count(abs(cHalfPi-Theta(1:nTheta,1)) > LatBoundary)
  
  ! The pole is a boundary point
  if(north)then
     iMin = nTheta - floor(maxval(HighLatBoundary)); iMax = nThetaUsed+1  
  else
     iMin = nTheta - nThetaUsed; iMax = ceiling(maxval(HighLatBoundary))
  end if
  if (north) then
     do j=1, nPsi-1
        iHighBnd(j) = nTheta - floor(HighLatBoundary(j)) 
     end do
  else
     do j=1, nPsi-1
        iHighBnd(j) = ceiling(HighLatBoundary(j)) 
     end do
  end if

  nThetaSolver = iMax - iMin + 1

  ! get the weimer potential for the high polar cap region
  lat_weimer(1:nTheta,:) = 90. - Theta(1:nTheta,:)*cRadToDeg
  mlt_weimer(1:nTheta,:) = Psi(1:nTheta,:)*12.0/cPi
  where(mlt_weimer <= 12)
     mlt_weimer = mlt_weimer + 12 
  elsewhere
     mlt_weimer = mlt_weimer - 12  
  end where

  ! symmetric for southern/northern hemisphere: lat_weimer is positive.
  call iono_potential_weimer(abs(lat_weimer(1:nTheta,1:nPsi-1)), mlt_weimer(1:nTheta,1:nPsi-1),&
       PhiIono_weimer(1:nTheta,1:nPsi-1))
  PhiIono_weimer(:,nPsi) = PhiIono_weimer(:,1)

!  write(NameFile, '(a,i4.4,2i2.2,"-",3i2.2,"-",i3.3,"_b",i1,a)')&
!       trim(PathSceOut)//'Weimer_potential_', time_array(1:7),iBlock,'.txt'
!  open(unittmp_, file=NameFile)
!  do j=1,nPsi
!     do i=1,nTheta
!        write(unittmp_, '(3e13.5, i5)')lat_weimer(i,j),mlt_weimer(i,j),PhiIono_weimer(i,j)
!     end do
!  end do
!
!  PhiIono_weimer = 10.0 !!! testing

  write(*,*)'max min(weimer):',maxval(PhiIono_weimer), minval(PhiIono_weimer)

  ! Calculate (Delta Psi)^2 (note that dPsi = 2 Delta Psi)
  dPsi2 = (dPsi/2.0)**2

  ! Calculate (Delta Theta)^2  (dTheta = 2 Delta Theta for i=1 and nTheta)
  dTheta2        = (dTheta/2.0)**2
  dTheta2(1)     = 4*dTheta2(1)
  dTheta2(nTheta)= 4*dTheta2(nTheta)

  SinTheta_I = sin(Theta(:,1))
  CosTheta_I = cos(Theta(:,1))

  if(DoTest .and. north) &
       write(*,*)'north:nThetaUsed, LatBoundary, Lat(iMax)=',&
       nThetaUsed,LatBoundary*cRadToDeg,&
       abs(cHalfPi-Theta(iMax,1))*cRadToDeg
  if(DoTest .and. .not. north) &
       write(*,*)'south:nThetaUsed, LatBoundary, Lat(iMin)=',&
       nThetaUsed,LatBoundary*cRadToDeg,&
       abs(cHalfPi-Theta(iMin,1))*cRadToDeg

  nX = nPsiUsed*nThetaSolver

  allocate( x(nX), y(nX), rhs(nX), b(nX), Bnd_I(nX), &
       d_I(nX), e_I(nX), e1_I(nX), f_I(nX), f1_I(nX))

  do j = 1, nPsiUsed
     do i= iMin, iMax
        ! set the matrix to identity for regions with Weimer potential
        if (north .and. i <= iHighBnd(j) .or. .not. north .and. i >= iHighBnd(j))then
           C_A(i,j) = 1.
           C_B(i,j) = 0.
           C_C(i,j) = 0.
           C_D(i,j) = 0.
           C_E(i,j) = 0.
        else
           sn  = SinTheta_I(i)
           cs  = CosTheta_I(i)
           sn2 = sn**2
           
           ! Central difference coefficients for second and first derivatives
           TermTheta2 = SigmaThTh(i,j)*sn2/dTheta2(i)
           TermTheta1 = ( SigmaThTh(i,j)*sn*cs   &
                + dSigmaThTh_dTheta(i,j)*sn2 &
                - dSigmaThPs_dPsi(i,j)*sn) / dTheta(i)
           TermPsi2 = SigmaPsPs(i,j) / dPsi2(j)            
           TermPsi1 = (dSigmaThPs_dTheta(i,j)*sn + dSigmaPsPs_dPsi(i,j)) / dPsi(j)
           
           ! Form the complete matrix
           C_A(i,j) = -2.0 * (TermTheta2 + TermPsi2)
           C_B(i,j) =         TermTheta2 - TermTheta1
           C_C(i,j) =         TermTheta2 + TermTheta1
           C_D(i,j) =         TermPsi2   - TermPsi1
           C_E(i,j) =         TermPsi2   + TermPsi1        

!!!         ! Form the complete matrix, testing with the simple laplace matrix
!           C_A(i,j) =-4.0
!           C_B(i,j) = 1.0
!           C_C(i,j) = 1.0
!           C_D(i,j) = 1.0
!           C_E(i,j) = 1.0
        end if
     end do
  enddo

  ! Fill in the diagonal vectors
  iI = 0
  do j = 1, nPsiUsed; do i=iMin,iMax
     iI = iI + 1
     ! The right-hand-side is Jr * Radius^2 sin^2(Theta).
     b(iI)    = Jr(i,j)*(Radius*SinTheta_I(i))**2
     d_I(iI)  = C_A(i,j)
     e_I(iI)  = C_B(i,j)
     f_I(iI)  = C_C(i,j)
     e1_I(iI) = C_D(i,j)
     f1_I(iI) = C_E(i,j)

     if(i == iMin)     e_I(iI)  = 0.0
     if(i == iMax)     f_I(iI) = 0.0
     if(j == 1)        e1_I(iI) = 0.0
     if(j == nPsiUsed) f1_I(iI) = 0.0
  end do; end do

  if(DoTest) write(*,*)'north,sum(b,abs(b),d,e,f,e1,f1)=',&
       north,sum(b),sum(abs(b)),sum(d_I),sum(e_I),sum(f_I),&
       sum(e1_I),sum(f1_I)

  ! Save original RHS for checking
  !b = 0 !!! test for zero RHS
  Rhs = b

  ! move nonlinear part of operator to RHS
  x = 0.0
  UseWeimer = .True.
  DoPrecond = .False.
  call matvec_ionosphere(x, Bnd_I, nX)
  b = b - Bnd_I
  UseWeimer = .False.
!  write(*,*)'C_A(iMin:iMax,1):', C_A(iMin:iMax,1)
!  write(*,*)'b(iMin:iMax)     ', b(1:iMax-iMin+1)
!  write(*,*)'Bnd_I(iMin:iMax) ', Bnd_I(1:iMax-iMin+1)

  DoPrecond = UsePreconditioner
  if(DoPrecond)then
     ! A -> LU
     call prehepta(nX,1,nThetaSolver,nX,real(-0.5,kind=8),d_I,e_I,f_I,e1_I,f1_I)
     ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs
     if(DoTest)write(*,*)'after precond: north,sum(b,abs(b),x,d,e,f,e1,f1)=',&
          north,sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
          sum(e1_I),sum(f1_I)

     ! rhs'=U^{-1}.L^{-1}.rhs
     call Lhepta(       nX,1,nThetaSolver,nX,b,d_I,e_I,e1_I)
     call Uhepta(.true.,nX,1,nThetaSolver,nX,b,f_I,f1_I)
  end if

  if(DoTest)write(*,*)'after precond: north,sum(b,abs(b),x,d,e,f,e1,f1)=',&
       north,sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
       sum(e1_I),sum(f1_I)


  ! Solve A'.x = rhs'
  Residual    = Tolerance
  if(UseInitialGuess) then
     iI = 0
     do j = 1, nPsiUsed; do i=iMin,iMax
        iI = iI + 1
        x(iI)    = PhiOld_CB(i,j,iBlock)
     end do; end do
  else
     x = 0.0
  end if

  select case(NameSolver)
  case('gmres')
     nIteration = MaxIteration
     call gmres(matvec_ionosphere, b, x, UseInitialGuess, nX, MaxIteration, &
          Residual, 'abs', nIteration, iError, DoTestMe)
  case('bicgstab')
     nIteration = 3*MaxIteration
     call bicgstab(matvec_ionosphere, b, x, UseInitialGuess, nX, &
          Residual, 'abs', nIteration, iError, DoTestMe)
  case default
     call CON_stop(NameSub//': unknown NameSolver='//NameSolver)
  end select

  if(DoTest .or. (iError /= 0 .and. iError /=3) ) &
       write(*,*)'iono_solve: north, iter, resid, iError=',&
       north, nIteration, Residual, iError
  if(iError /= 0 .and. iError /=3)then
     write(*,*)'IE_ERROR in iono_solve: solver failed !!!'
     if(iError < 0) &
          call CON_stop('IE_ERROR in iono_solve: residual did not decrease')
  end if

  !Check solution:
  if(DoTest)then
     DoPrecond = .false.
     UseWeimer = .true.
     call check_solution(x, y, rhs, nX)
  end if

  ! Phi_C is the solution within the solved region
  Phi_C(:,:) = 0.0
  iI = 0
  do j=1, nPsiUsed; do i = iMin, iMax
     iI = iI + 1
     if (north)then
        if (i>=iHighBnd(j)) Phi_C(i,j) = x(iI)     
     else
        if (i<=iHighBnd(j)) Phi_C(i,j) = x(iI)
     end if
  end do; end do
  ! Apply periodic boundary condition in Psi direction
  Phi_C(:,nPsi) = Phi_C(:,1)
  PhiMax = maxval(Phi_C)
  PhiMin = minval(Phi_C)

  ! Save the solution for next time
  PhiOld_CB(:,:,iBlock) = Phi_C
  if(north)then
     do j=1, nPsi-1
        PhiOld_CB(1:iHighBnd(j), j, iBlock) = PhiIono_Weimer(1:iHighBnd(j),j)
     end do
  else
     do j=1, nPsi-1
        PhiOld_CB(iHighBnd(j):nTheta,j, iBlock) = PhiIono_Weimer(iHighBnd(j):nTheta,j)
     end do
  end if
  ! apply periodic boundary condition in Psi direction
  PhiOld_CB(:,nPsi,iBlock) = PhiOld_CB(:,1,iBlock)

  if (north) then 
     ! Apply average condition at north pole
     !Phi_C(1,:) = sum(Phi_C(2,1:nPsiUsed))/nPsiUsed
     cpcp_north = (PhiMax - PhiMin)/1000.0
     call write_prefix; write(iUnitOut,*) &
          "iono_solver: Northern Cross Polar Cap Potential=",&
          cpcp_north," kV"
  else
     ! Apply average condition at south pole
     !Phi_C(nTheta,:) = sum(Phi_C(nTheta-1,1:nPsiUsed))/nPsiUsed
     cpcp_south = (PhiMax - PhiMin)/1000.0
     call write_prefix; write(iUnitOut,*) &
          "iono_solver: Southern Cross Polar Cap Potential=",&
          cpcp_south," kV"
  endif
  if(allocated(b)) deallocate(x, y, b, rhs,Bnd_I, d_I, e_I, f_I, e1_I, f1_I)

  write(NameFile, '(a,i4.4,2i2.2,"_",3i2.2,"_",i3.3,"_b",i1,a)')&
       trim(PathSceOut)//'Output_potential_', time_array(1:7),iBlock,'.txt'
  open(unittmp_, file=NameFile)
  if (north)then
     write(unittmp_,'(a)')'Northern Hemisphere'
     write(unittmp_,'(i4.4,2i2.2,"_",3i2.2,"_",i3.3)') time_array(1:7)
  else
     write(unittmp_,'(a)')'Southern Hemisphere'
  end if
  do i=1,nTheta
     do j=1,nPsi
        write(unittmp_, '(5e13.5)')Theta(i,j),Psi(i,j),Phi_C(i,j), PhiOld_CB(i,j,iBlock), HighLatBoundary(j)
     end do
  end do

contains
  !============================================================================
  subroutine check_solution(x_I, y_I, b_I, n)

    use IE_ModMain, ONLY: Time_Simulation, nSolve
    use ModIoUnit,  ONLY: UnitTmp_

    ! Calculate y = A.x where A is the (pentadiagonal) matrix

    integer, intent(in) :: n          ! number of unknowns
    real(real8_), intent(in) :: x_I(n)        ! vector of unknowns
    real(real8_), intent(out):: y_I(n)        ! y = A.x
    real(real8_), intent(in) :: b_I(n)        ! rhs

    integer :: iTheta, iPsi, i, iError(1)
    !-------------------------------------------------------------------------
    call matvec_ionosphere(x_I, y_I, n)

    if(north)then
       open(UnitTmp_, FILE='iono_north.out')
    else
       open(UnitTmp_, FILE='iono_south.out')
    end if
    write(UnitTmp_,'(a79)')'Iono Solver_var22'
    write(UnitTmp_,'(i7,es13.5,3i3)') nSolve,Time_Simulation,2,1,13
    write(UnitTmp_,'(2i4)') nThetaSolver, nPsiUsed
    write(UnitTmp_,'(es13.5)') 0.0
    write(UnitTmp_,'(a79)') &
         'Theta Psi A B C D E d e f e1 f1 Solution Rhs Error Param'

    i = 0
    iError = maxloc(abs(y_I-b_I))
    do iPsi = 1, nPsiUsed; do iTheta = iMin, iMax; i = i+1
       if(i==iError(1))write(*,*) &
            '!!! Check: max(abs(y-b)), iPsi, iTheta, nPtsTheta, north=',&
            iPsi, iTheta, nThetaSolver, abs(y_I(iError)-b_I(iError)), north
       write(UnitTmp_,'(100es18.10)') &
            cRadToDeg*Theta(iTheta,1), cRadToDeg*Psi(1,iPsi), &
            C_A(iTheta, iPsi),  C_B(iTheta, iPsi),  C_C(iTheta, iPsi), &
            C_D(iTheta, iPsi),  C_E(iTheta, iPsi),  &
            d_I(i), e_I(i), f_I(i), e1_I(i), f1_I(i), &
            y_I(i), b_I(i), y_I(i)-b_I(i)
    end do; end do
    close(UnitTmp_)

  end subroutine check_solution

end subroutine ionosphere_solver

!============================================================================
subroutine matvec_ionosphere(x_I, y_I, n)

  use ModIonosphere,   ONLY: IONO_nPsi, IONO_nTheta, nThetaUsed,nThetaSolver, &
                             HighLatBoundary, north, DoPrecond, d_I, e_I, f_I, &
                             e1_I, f1_I, C_A, C_B, C_C, C_D, C_E, PhiIono_Weimer, &
                             UseWeimer, nTheta_sol, nPsi_sol, iHighBnd
  use ModLinearsolver, ONLY: Uhepta, Lhepta
  use ModKind
  implicit none

  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1

  ! Calculate y = A.x where A is the (pentadiagonal) matrix

  integer, intent(in) :: n          ! number of unknowns
  real(real8_), intent(in) :: x_I(n)        ! vector of unknowns
  real(real8_), intent(out):: y_I(n)        ! y = A.x

  integer :: iTheta, iTheta2, iPsi, i, iMin, iMax,j
  real(real8_) :: x_G(0:nTheta+1, 0:nPsi) ! 2D array with ghost cells
  !-------------------------------------------------------------------------

  if(north)then
     iMin = nTheta - floor(maxval(HighLatBoundary)); iMax = nThetaUsed+1  
  else
     iMin = nTheta - nThetaUsed; iMax = ceiling(maxval(HighLatBoundary))
  end if

  nThetaSolver = iMax - iMin + 1

  x_G = 0.0

  ! Put 1D vector into 2D solution
  i = 0;
  do iPsi = 1, nPsi-1; do iTheta = iMin, iMax; i = i+1
     x_G(iTheta, iPsi) = x_I(i)
  enddo; enddo

  if(north)then
     if(UseWeimer)then ! apply the boundary condition at high latitude
        do j=1, nPsi-1
           x_G(iMin-1:iHighBnd(j), j) = PhiIono_Weimer(iMin-1:iHighBnd(j),j)
        end do
     else
        do j=1, nPsi-1
           x_G(iMin-1:iHighBnd(j), j) = 0.0
        end do
     end if
     x_G(iMax+1,1:nPsi-1) = 0.0
  else
     if(UseWeimer)then
        do j=1, nPsi-1
           x_G(iHighBnd(j):iMax+1,j) = PhiIono_Weimer(iHighBnd(j):iMax+1,j)
        end do
     else
        do j=1, nPsi-1
           x_G(iHighBnd(j):iMax+1,j) = 0.0
        end do
     end if
     x_G(iMin-1,1:nPsi-1) = 0.0 
  end if

!  write(*,*)'PhiIono_weimer(iMin-1:iHighBnd(npsi/2),npsi/2):',PhiIono_weimer(iMin-1:iHighBnd(nPsi/2),nPsi/2)
!  write(*,*)'x_G(iMin-1:iMax+1,nPsi/2)',x_G(iMin-1:iMax+1,nPsi/2)
  
  ! Apply periodic boundary conditions in Psi direction
  x_G(:,nPsi) = x_G(:,1)
  x_G(:,0)    = x_G(:,nPsi-1)

  i = 0
  do iPsi = 1, nPsi-1; do iTheta = iMin, iMax; i = i+1     
     y_I(i) = &
          C_A(iTheta, iPsi)*x_G(iTheta,   iPsi)   + &
          C_B(iTheta, iPsi)*x_G(iTheta-1, iPsi)   + &
          C_C(iTheta, iPsi)*x_G(iTheta+1, iPsi)   + &
          C_D(iTheta, iPsi)*x_G(iTheta,   iPsi-1) + &
          C_E(iTheta, iPsi)*x_G(iTheta,   iPsi+1)
  end do; end do



  if(UseWeimer)then ! apply the boundary condition at high latitude
     i = 0
     do iPsi=1, nPsi-1; do iTheta = iMin, iMax; i = i + 1
        if(north)then
           if (iTheta <= iHighBnd(iPsi)) y_I(i) = 0.0
        else
           if (iTheta >= iHighBnd(iPsi)) y_I(i) = 0.0
        end if
!        if (iPsi == 1 .and. iTheta == iMax)then 
!           write(*,*)'y_I(i) = ', y_I(i)
!           write(*,*)'C_A, x_G, C_B...=', &
!                C_A(iTheta, iPsi), x_G(iTheta,   iPsi)   , &
!                C_B(iTheta, iPsi), x_G(iTheta-1, iPsi)   , &
!                C_C(iTheta, iPsi), x_G(iTheta+1, iPsi)   , &
!                C_D(iTheta, iPsi), x_G(iTheta,   iPsi-1) , &
!                C_E(iTheta, iPsi), x_G(iTheta,   iPsi+1)
!        end if
     end do; end do

  end if

  ! Preconditioning: y'= U^{-1}.L^{-1}.y
  if(DoPrecond)then
     call Lhepta(       n,1,nThetaSolver,n,y_I,d_I,e_I,e1_I)
     call Uhepta(.true.,n,1,nThetaSolver,n,y_I,    f_I,f1_I)
  end if

end subroutine matvec_ionosphere

!========================================================================
subroutine iono_potential_weimer(lat, mlt, PhiWeimer)
  
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi, IONO_NORTH_Theta, IONO_NORTH_Psi, &
                           IONO_SOUTH_Theta, IONO_SOUTH_Psi, cRadToDeg, cPi
  USE ModIoUnit,     ONLY: UNITTMP_
  use ModRamIndices, ONLY: NameOmniFile
  use ModKind
  use w05
!  use Module1, ONLY: ist3, tilt

  use ModRamTiming,    ONLY: TimeRamElapsed
  use ModScbVariables, ONLY: tilt

  implicit none

  integer :: ist3
  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi
  real(real8_), intent(in) :: lat(nTheta,nPsi), mlt(nTheta,nPsi)
  real(real8_), intent(out) :: PhiWeimer(nTheta, nPsi)
  character(LEN = 100) :: header
  integer :: i,j,ierr
  integer :: iYear_l, iDoy_l, iHour_l, iMin_l, iLines
  real(real8_) :: bzimfArr_l, byimf_l, pdyn_l, Nk_l, Vk_l, bTot_l
  real(real8_) :: angle, bzimf_l
  integer :: symh_l,al_l
  real(real8_), external :: EpotVal
  logical :: UseAL
  !---------------------------------------------------------------------
     ist3 = TimeRamElapsed/300

     open(Unittmp_, file=NameOmniFile, status = 'unknown',action='read')
     i = 0
     ierr = 0
     do while (ierr == 0)
        read(unittmp_, *, iostat=ierr) header
        i = i + 1
     end do
     ilines = i-1

     rewind(unittmp_)
!     print*, 'IP: year, day, hour, min, bt, by, bz, v, n, pdyn, al. symh'
     read_data_from_file1: do i=1,ilines
        read(unittmp_,101) iyear_l,idoy_l,ihour_l, imin_l, btot_l, byimf_l, bzimf_l, vk_l, nk_l, pdyn_l, al_l, symh_l
        if(i==ist3+1)then ! the omni data is saved every 5 minutes, so the same as the iteration ist3
!           write(*,101)iyear_l,idoy_l,ihour_l, imin_l, btot_l, byimf_l, bzimf_l, vk_l, nk_l, pdyn_l, al_l, symh_l
           exit read_data_from_file1
        end if
     end do read_data_from_file1
     close(unittmp_)
101  FORMAT(2I4, 2I3, 3F8.2, F8.1, F7.2, F6.2, 2I6)
     
     UseAL = .true.
     angle = cRadToDeg*ACOS(bzimf_l/bTot_l)

     ! use the weimer 2001 model
     CALL SetModel(angle, bTot_l, tilt, Vk_l, Nk_l, real(Al_l,kind=8), UseAL); call flush(6)
     ! use the weimer 2005 model
!     CALL SetModel05(byimf_l, bzimf_l, tilt, Vk_l, Nk_l)

     ! calculate weimer potential from w05 for the high-latitude potential used in the solver
     DO j = 1, nPsi
        DO i = 1, nTheta
           IF (ABS(pdyn_l-99.99) < 1E-3) THEN
              PRINT*, 'IP: not calling Weimer model, bad SW data'
              EXIT ! Bad SW data, do not call Weimer model
           END IF
           ! use weimer 2001 model
           PhiWeimer(i,j) =  EpotVal(lat(i,j), mlt(i,j))*1.0e3 ! in Volts
           ! use the weimer 2005 model
!           call EpotVal05(lat(i,j), mlt(i,j), real(0.0,kind=8), PhiWeimer(i,j))
!           PhiWeimer(i,j) = PhiWeimer(i,j) * 1.0e3 ! in Volts

        END DO
     END DO

   end subroutine iono_potential_weimer
