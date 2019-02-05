!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!

!============================================================================
module ModLinearSolverWrapper

  use ModLinearSolver, only: LinearSolverParamType
  implicit none
  type(LinearSolverParamType), public :: Param
  integer, public :: iBlock, nVarIjk
  integer, public :: nDim       ! Number of spatial dimensions
  integer, public :: nVar       ! Number of impl. variables/cell
  integer, public :: nI, nJ, nK ! Number of cells in a block
  integer, public :: nBlock     ! Number of impl. blocks on current proc
  real, pointer,   public :: PrecondMatrix_II(:,:)

  contains
    subroutine linear_solver_matvec(x_I, y_I, n)

      ! Fortran subroutine calling C matvec routine
      use iso_c_binding
      use ModLinearSolver, ONLY: multiply_left_precond


      implicit none

      interface
         subroutine linear_wrapper_matvec_c(x_I, y_I, n)  bind(C) 
           integer, VALUE:: n
           real          :: x_I(n)
           real          :: y_I(n)
         end subroutine linear_wrapper_matvec_c
      end interface

      integer, intent(in) :: n
      real,    intent(in) :: x_I(n) 
      real,    intent(out):: y_I(n)

      integer :: iStart
      !-----------------------------------------------------------------------

      call linear_wrapper_matvec_c(x_I, y_I, n)

      if(Param%DoPrecond) then
         do iBlock = 1, nBlock

            ! Starting index in the linear arrays
            iStart = nVarIjk*(iBlock-1)+1
            call multiply_left_precond(&
                 Param%TypePrecond, Param%TypePrecondSide, &
                 nVar, nDim, nI, nJ, nK, PrecondMatrix_II, &
                 y_I(iStart))
         end do
      end if

    end subroutine linear_solver_matvec
    !=========================================================================

  end module ModLinearSolverWrapper
!============================================================================
subroutine linear_solver_wrapper(TypeSolver, Tolerance, nIteration, &
     nVarIn, nDimIn, nIIn, nJIn, nKIn, nBlockIn, iComm, Rhs_I, x_I, &
     PrecondParam, PrecondMatrixIn_II, lTest) bind(C)

  ! Wrapper to call ModLinearSolver methods from C/C++ code

  use iso_c_binding
  use ModLinearSolverWrapper
  use ModLinearSolver, ONLY: get_precond_matrix, multiply_left_precond,&
       multiply_initial_guess, bicgstab,gmres,cg,LinearSolverParamType
  use ModUtilities, ONLY: char_array_to_string, CON_stop
 
  implicit none

  character(c_char), intent(in) :: TypeSolver(*) 
  real,    intent(in):: Tolerance  ! tolerance for the solver
  integer, intent(in):: nIteration ! max iteration number
  integer, intent(in):: nVarIn       ! Number of impl. variables/cell
  integer, intent(in):: nDimIn       ! Number of spatial dimensions
  integer, intent(in):: nIIn, nJIn, nKIn ! Number of cells in a block
  integer, intent(in):: nBlockIn     ! Number of impl. blocks on current proc
  integer, intent(in):: iComm      ! MPI communicator for processors
  real,    intent(in):: PrecondParam ! Parameter for the preconditioner
  integer, intent(in):: lTest      ! if lTest==1 DoTest=.true.
  real, target, intent(inout):: PrecondMatrixIn_II(nVarIn*nVarIn*nIIn*nJIn*nKIn,2*nDimIn+1)
  !precond_matrix contains diagonal and super/sub diagonal elements from
  !the matrix A, which is in the equation Ax = b 

  real, intent(inout):: Rhs_I(nVarIn*nIIn*nJIn*nKIn*nBlockIn) ! RHS vector
  real, intent(inout):: x_I(nVarIn*nIIn*nJIn*nKIn*nBlockIn)   ! Initial guess/solution

  logical :: DoTest     ! show Krylov iterations and convergence

  ! Local variables
  integer:: n, i, j, k, iVar
  integer:: nImpl

  logical:: DoDebug = .false.
  character(len=*), parameter:: NameSub = 'linear_solver_wrapper'

  !---------------------------------------------------------------------------
  
  nBlock = nBlockIn
  nVar   = nVarIn
  nDim   = nDimIn
  nI     = nIIn
  nJ     = nJIn
  nK     = nKIn
  PrecondMatrix_II => PrecondMatrixIn_II

  ! PrecondParam:      The parameter for Gustafsson modification:
  !
  !           +2 DILU  prec: LU for diagonal, keep off-diagonal blocks
  !           +1 BILU  prec: LU for diagonal, premultiply U with D^-1
  !           <0 MBILU prec: Gustaffson modification of diagonal blocks
  !                          using -1 <= PrecondParam < 0 parameter

  call char_array_to_string(TypeSolver, Param%TypeKrylov)

  Param%PrecondParam = PrecondParam

  Param%DoPrecond = .true.

  if(PrecondParam == 0.0)then
     Param%DoPrecond = .false.
  else if(abs(PrecondParam - 2.0) < 1e-6)then
     Param%TypePrecond = 'DILU'
  else if(abs(PrecondParam - 1.0) < 1e-6)then
     Param%TypePrecond = 'BILU'
  else if(PrecondParam < 0.0 .and. PrecondParam >= -1.0)then
     Param%TypePrecond = 'MBILU'
  else
     call CON_stop(NameSub//': invalid value for PrecondParam=', PrecondParam)
  end if

  if (Param%DoPrecond .and. Param%TypeKrylov == 'CG') &
       call CON_stop(NameSub//' CG solver does not have preconditioner')

  Param%TypeStop      = 'rel'
  Param%ErrorMax      = Tolerance
  Param%MaxMatvec     = nIteration
  Param%nKrylovVector = nIteration
  Param%UseInitialGuess = .false. 

  ! Number of variables per block
  nVarIjk = nVar*nI*nJ*nK

  ! Number of variables per processor
  nImpl   = nVarIjk*nBlock

  DoTest = lTest==1 

  ! Make sure that left preconditioning is used when necessary
  Param%TypePrecondSide = 'left'

  ! Initialize solution vector to zero
  x_I = 0.0

  ! Get preconditioning matrix if required. 
  ! Precondition RHS and initial guess (for symmetric prec only)
  if(Param%DoPrecond)then
     ! Preconditioning  matrix
     call get_precond_matrix(                         &
          Param%PrecondParam, nVar, nDim, nI, nJ, nK, PrecondMatrix_II)

     do iBlock = 1, nBlock
        ! Starting index in the linear arrays
        n = nVarIjk*(iBlock-1)+1

        ! rhs --> P_L.rhs, where P_L=U^{-1}.L^{-1}, L^{-1}, or I
        ! for left, symmetric, and right preconditioning, respectively
        call multiply_left_precond(&
             Param%TypePrecond, Param%TypePrecondSide, &
             nVar, nDim, nI, nJ, nK, PrecondMatrix_II, &
             Rhs_I(n))
     end do
  endif

  ! Initialize stopping conditions. Solver will return actual values.
  Param%nMatVec = Param%MaxMatvec
  Param%Error   = Param%ErrorMax

  if(DoTest)write(*,*)NameSub,': Before ', Param%TypeKrylov, &
       ' nMatVec, Error:', Param%nMatVec, Param%Error

  ! Solve linear problem
  !call timing_start('krylov solver')
  select case(Param%TypeKrylov)
  case('BICGSTAB')
     call bicgstab(linear_solver_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl, &
          Param%Error, Param%TypeStop, Param%nMatvec, &
          Param%iError, DoTest, iComm)
  case('GMRES')
     call  gmres(linear_solver_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl, &
          Param%nKrylovVector, &
          Param%Error, Param%TypeStop, Param%nMatvec, &
          Param%iError, DoTest, iComm)
  case('CG')
     if(.not. Param%DoPrecond)then
        call cg(linear_solver_matvec, Rhs_I, x_I, Param%UseInitialGuess, nImpl, &
             Param%Error, Param%TypeStop, Param%nMatvec, &
             Param%iError, DoTest, iComm)      
     end if
  case default
     call CON_stop(NameSub//': Unknown TypeKrylov='//Param%TypeKrylov)
  end select
  !call timing_stop('krylov solver')

  if(DoTest)write(*,*)NameSub,&
       ': After nMatVec, Error, iError=',&
       Param%nMatvec, Param%Error, Param%iError

  ! Converging without any iteration is not a real error, so set iError=0
  if(Param%iError==3) Param%iError=0

end subroutine linear_solver_wrapper
!=========================================================================
