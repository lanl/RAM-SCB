!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

module ModRamNCDF

  implicit none
  save
  
  contains
!==============================================================================
  subroutine ncdf_check(iStatusIn, NameSubIn, NameFileIn)
    ! Check the status flag of a netcdf operation, CON_stop as necessary.

    use netcdf
    
    integer, intent(in) :: iStatusIn
    character(len=*), intent(in) :: NameSubIn
    character(len=*), intent(in), optional :: NameFileIn

    character(len=*), parameter :: NameSub = 'ncdf_check'
    !------------------------------------------------------------------------
    if(iStatusIn /= nf90_noerr) then
       write(*,*) 'ERROR WITH NETCDF:'
       write(*,*) '      Subroutine: ', trim(NameSubIn)
       if(present(NameFileIn)) &
            write(*,*) '      File: ', trim(NameFileIn)
       write(*,*) '      Error code: ', iStatusIn
       write(*,*) '      ', trim(nf90_strerror(iStatusIn))
       call CON_stop(NameSub // ': NetCDF error.')
    end if

  end subroutine ncdf_check

!==============================================================================
  subroutine write_ncdf_globatts(iFileID)
    ! For a netcdf file that is still in define mode, write a plethora of 
    ! global attributes that all RAM-SCB output NetCDFs should have.

    use netcdf
    use ModRamTiming, ONLY: TimeRamStart
    use ModRamParams, ONLY: StrRamDescription

    integer, intent(in) :: iFileID  ! NetCDF file ID number for open file.

    integer :: iStatus
    character(len=100) :: NameSub = 'write_ncdf_globatts'
    !------------------------------------------------------------------------
    ! Description of run:
    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'description', &
         StrRamDescription)
    call ncdf_check(iStatus, NameSub)

    ! Time of run:
!    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'run_date', &
!         StringRunDate)
!    call ncdf_check(iStatus, NameSub)

    ! Start time of run:
    iStatus = nf90_put_att(iFileID, NF90_GLOBAL, 'start_time', &
         TimeRamStart % String)
    call ncdf_check(iStatus, NameSub)

  end subroutine write_ncdf_globatts

!==============================================================================
!  subroutine write_ncdf_var(iFileID, varChar, varDim, var)
!      
!    use netcdf
!
!    implicit none
!
!    integer, intent(in) :: iFileID
!    character(len=*), intent(in) LL varChar
!
!    integer :: iVar
!    iStatus = nf90_def_var(iFileID, varChar, nf90_double, varDim, iVar)
!    call ncdf_check(iStatus, 'Define variable: '//varChar)
!    iStatus = nf90_put_var(iFileID, iVar, var)
!    call ncdf_check(iStatus, 'Put variable: '//varChar)
!
!  end subroutine write_ncdf_var
!==============================================================================

END MODULE ModRamNCDF
