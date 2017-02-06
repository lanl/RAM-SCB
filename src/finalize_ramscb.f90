!******************************************************************************
subroutine finalize_ramscb()
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!******************************************************************************
  
  use ModRamTiming, ONLY: finalize_timing
  use ModRamMain,   ONLY: TimeRamElapsed, iCal, IsComponent, DoSaveFinalRestart
  use ModRamMpi,    ONLY: iProc
  use ModRamIO,     ONLY: write_restart

  implicit none
  !---------------------------------------------------------------------------
  if((.not.IsComponent) .and. (iProc==0)) then
     write(*,*) &
          '==============================================================='
     write(*,*) '                     RAM_SCB finished!'
     
     write(*,'(f12.2,a,i10.10,a)') &
          TimeRamElapsed, 's simulated in ', iCal-1, ' iterations.'

     if(DoSaveFinalRestart)call write_restart

     ! Stop timing Ram; write final timing report.
     call finalize_timing

     write(*,*) &
          '==============================================================='
  end if

end subroutine finalize_ramscb
