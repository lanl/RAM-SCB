!============================================================================
module ModRamMpi
  ! Hold number processors, RAM-SCB Mpi communicator, and processor number.
  ! These values are set one of two ways:
  ! 1- by Main.f90 (standalone only)
  ! 2- by the SWMF (component mode only)
  ! In either case, accessing the values here will give you the correct
  ! values for running RAM-SCB as a parallelized code.
  !    Copyright (c) 2016, Los Alamos National Security, LLC
  !    All rights reserved.

  implicit none; save
  integer::iProc,nProc,iComm,iGroup,iError


end module ModRamMpi

!============================================================================
