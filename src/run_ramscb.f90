subroutine run_ramscb
!==============================================================================
!    Call RAM for two timesteps, run SCB model and couple as necessary.
!    This file replaces the later half of Main.f.
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!==============================================================================
  
  use ModRamMain, ONLY: iCal, DTsMax,DTsNext, Kp, UTs, DTs, Dt_hI,&
       PPerE,PParE, PPerO, PParO, PPerH, PParH, PPerHe, PParHe, PPerT, PParT, &
       event, electric, iPressure, iDomain, TimeRamElapsed, TimeMax, &
       TimeRamNow, TimeRamStart, Real8_, DoUseScb, NR, NT, DtsFramework, &
       DoVarDt, DtsMin, f107, electrons, nIter
  use ModRamTiming,  ONLY: do_timing
  use ModRamIndices, ONLY: get_indices
  use ModRamIO,      ONLY: max_output_timestep, write_prefix
  use ModRamFunctions,ONLY: ram_sum_pressure
  use ModTimeConvert
  use ModRamMpi
  use ModMpi

  implicit none

  real(kind=Real8_) :: DtOutputMax, DtEndMax

  character(len=*), parameter :: NameSub='run_ramscb'
  character(len=4)            :: StrScbIter = '0001'
  !---------------------------------------------------------------------------

  ! Set UTs to TimeRamElapsed + starting hour.
  UTs = TimeRamElapsed + TimeRamStart%iHour*3600.0

  ! Set the half-timestep based on CFL, SWMF, and Maximum allowed.
  ! Use CFL Number if #VARIABLEDT is set in PARAM file.
  if(DoVarDt) then
     DtEndMax   =(TimeMax-TimeRamElapsed)/2.0
     DtOutputMax=max_output_timestep(TimeRamElapsed)
     DTs=min(DTsNext, DTsmax, DtsFramework, DtOutputMax, DtEndMax)
     if(Kp.gt.6.0 .AND. DTs.gt.1.) DTs = 1.   !1. or 0.25 
 else if(mod(UTs, Dt_hI) .eq. 0) then
     DTs = 5.0
     if(Kp .ge. 5.0) DTs = min(DTsMin,DTs)
     if(Kp .gt. 6.0) DTs = 1.0   !1. or 0.25
  endif
     
  DTsNext = DTsmax

  ! Set indices for this timestep.
  call get_indices(TimeRamNow%Time, Kp, f107)
  if (ical.EQ.1) Kp=0.0

  ! Broadcast current call to ram_all
  if(iProc==0) then
     call write_prefix
     write(*,'(1x, a, 1x, F8.1, 3x, F5.1, 3x, F5.3)') &
          'Calling ram_all for UTs, DTs,Kp = ', UTs, DTs, Kp
  end if

  ! Call RAM for each species.
  ! Try to parallelize.
  if(nProc .lt. -1)then ! Disabled for the time being.
     if(iProc == 0) then
        call ram_all(4)  ! Oxygen
        PPerO = PPerT(4,:,:); PParO = PParT(4,:,:)
     end if
     if(iProc == 1) then
        call ram_all(3)  ! Helium
        PPerHe = PPerT(3,:,:); PParHe = PParT(3,:,:)
     end if
     if(iProc == 2) then   
      if (electrons) then
        call ram_all(1)  ! Electron
        PPerE = PPerT(1,:,:); PParE = PParT(1,:,:)
      end if
     end if
     if(iProc == 3) then   
        call ram_all(2)  ! Hydrogen
        PPerH = PPerT(2,:,:); PParH = PParT(2,:,:)
     end if
     call MPI_BCAST(PPerO, NR*NT, MPI_DOUBLE_PRECISION,0,iComm,iError)
     call MPI_BCAST(PParO, NR*NT, MPI_DOUBLE_PRECISION,0,iComm,iError)
     call MPI_BCAST(PPerHe,NR*NT, MPI_DOUBLE_PRECISION,1,iComm,iError)
     call MPI_BCAST(PParHe,NR*NT, MPI_DOUBLE_PRECISION,1,iComm,iError)
     call MPI_BCAST(PPerE, NR*NT, MPI_DOUBLE_PRECISION,2,iComm,iError)
     call MPI_BCAST(PParE, NR*NT, MPI_DOUBLE_PRECISION,2,iComm,iError)
     call MPI_BCAST(PPerH, NR*NT, MPI_DOUBLE_PRECISION,3,iComm,iError)
     call MPI_BCAST(PParH, NR*NT, MPI_DOUBLE_PRECISION,3,iComm,iError)
     ! What variables to share???
     ! Do that here!
  else if(iProc == 0) then
     call ram_all(4)  ! Oxygen
     PPerO = PPerT(4,:,:); PParO = PParT(4,:,:)
     call ram_all(3)  ! Helium
     PPerHe = PPerT(3,:,:); PParHe = PParT(3,:,:)
     if (electrons) then
       call ram_all(1)  ! Electron
       PPerE = PPerT(1,:,:); PParE = PParT(1,:,:)
     endif
     call ram_all(2)  ! Hydrogen
     PPerH = PPerT(2,:,:); PParH = PParT(2,:,:)
  end if

  ! Ensure all calls to ram_all have completed.
  call MPI_BARRIER(iComm, iError)

  ! Increment time
  TimeRamElapsed = TimeRamElapsed + 2.0 * DTs

  ! Update time
  TimeRamNow % Time = TimeRamStart % Time + TimeRamElapsed
  call time_real_to_int(TimeRamNow)

  ! Update B every five minutes using SCB model.
  if(mod(TimeRamElapsed, Dt_hI) .eq. 0)then
     if(iProc==0)then
        call write_prefix
        write(*,*) 'Running SCB model to update B-field...'
     end if
     
     ! Write Pressure file from Ram for SCB.
     write(StrScbIter,'(I4.4)') int(TimeRamElapsed/Dt_hI)
     call ram_sum_pressure
     
     if(iProc==0) call ram_write_pressure(StrScbIter)
     call MPI_BARRIER(iComm, iError)

     if(DoUseScb)then
        call parallel_3d_eq(UTs/3600.0D0, StrScbIter, iPressure, &
             iDomain, electric, event)
        if(iProc==0) then
           call write_prefix
           write(*,*) 'Finished 3D Equilibrium code.'
        end if
     end if
  end if

  ! Ensure that re-initialization doesn't occur.
  iCal = iCal + 1
  nIter=nIter + 1

  ! Do these on head node only:
  if(iProc == 0)then
     ! Do timing.
     call do_timing
     ! Check and write output, as necessary.
     call handle_output(TimeRamElapsed)
  end if
  
contains
  !===========================================================================
  subroutine ram_write_pressure(StringIter)
    
    use ModIoUnit,  ONLY: UNITTMP_
    use ModRamIO
    use ModRamMain, ONLY: NR, NT, Kp, PathRamOut, lz, phi, Pi, UTs, PPerO, &
         PParO, PPerH, PParH, PPerHe, PParHe, PPerE, PParE, PAllSum, &
         TimeRamNow, TimeRamStart, TimeRamElapsed, PparSum, DoAnisoPressureGMCoupling
    
    implicit none
    
    character(len=23) :: StringTime
    character(len=*), intent(in) :: StringIter
    character(len=*), parameter  :: NameSub = 'ram_write_pressure'
    character(len=200)           :: FileName
    integer                      :: iError, i, j
    !------------------------------------------------------------------------
    ! Create Ram pressure output file.
    if(UseNewFmt)then
       FileName = trim(PathRamOut)// &
            trim(RamFileName('/pressure','in',TimeRamNow))
    else
       FileName = trim(PathRamOut)//'/pressure_'//StringIter//'.in'
    end if
    open(unit=UNITTMP_, status='REPLACE', iostat=iError, file=FileName)
    if(iError /= 0) call CON_stop &
         (NameSub//' Error opening file '//FileName)
    
    ! Prepare date/time for header.
    write(StringTime,"(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,2(':',i2.2)'.',i3.3)") &
         TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
         TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond, &
         floor(TimeRamNow%FracSecond*1000.0)

    ! Write pressure to file.
    write(UNITTMP_,'(a, a, a3,f8.3,2x,a4,f3.1)') &
         'Date=', StringTime, ' T=', &
         TimeRamElapsed/3600.0 + TimeRamStart%iHour, ' Kp=', Kp
    if(.not.DoAnisoPressureGMCoupling)then
       write(UNITTMP_,'(2a)') '  Lsh    MLT     PPER_H      PPAR_H     ', &
            'PPER_O    PPAR_O    PPER_He    PPAR_He    PPER_E    PPAR_E   PTotal   [keV/cm3]'
    else
       write(UNITTMP_,'(2a)') '  Lsh    MLT     PPER_H      PPAR_H     ', &
            'PPER_O    PPAR_O    PPER_He    PPAR_He    PPER_E    PPAR_E   PTotal   [keV/cm3] Ppar [keV/cm3]'
       end if
    
    do i=2, NR; do j=1, NT
       if(.not.DoAnisoPressureGMCoupling)then
          write(UNITTMP_,'(f5.2,f8.2,9e12.3)') &
               LZ(I),PHI(J)*12/PI,PPERH(I,J),PPARH(I,J),PPERO(I,J),PPARO(I,J), &
               PPERHE(I,J),PPARHE(I,J),PPERE(I,J),PPARE(I,J),PAllSum(i,j)
       else
          write(UNITTMP_,'(f5.2,f8.2,10e12.3)') &
               LZ(I),PHI(J)*12/PI,PPERH(I,J),PPARH(I,J),PPERO(I,J),PPARO(I,J), &
               PPERHE(I,J),PPARHE(I,J),PPERE(I,J),PPARE(I,J),PAllSum(i,j),PparSum(i,j)
       end if
    enddo; enddo
    close(UNITTMP_)
    
  end subroutine ram_write_pressure
  
  !===========================================================================
  subroutine handle_output(TimeIn)
    ! For every type of output, check the following:
    ! 1) if that output is activated,
    ! 2) if it's time to write that output
    ! If so, call the proper routines to write the output.

    use ModRamMain, ONLY: TimeRamNow, nR, nT, DtRestart, Dt_hI, &
         PPerO, PParO, PPerH, PParH, PParHe, PPerHe, PPerE, PParE, VT
    use ModRamSats, ONLY: fly_sats
    use ModRamIO,   ONLY: iUnitLog, write_restart, DtLogfile, &
         DoSaveRamSats, DtWriteSat, NameFileLog, NameFileDst
    use ModRamIndices, ONLY: get_ramdst
    use Module1, ONLY : DstBiot, iUnitDst, DstDPS, DstDPSInsideGeo, &
         DstBiot, DstBiotInsideGeo

    implicit none

    real(kind=Real8_), intent(in) :: TimeIn

    real(kind=Real8_) :: dst
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSubSub = NameSub // '::handle_output'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSubSub, DoTest, DoTestMe)

    ! Get current Dst.
    call get_ramdst(dst)

    ! Write Logfile
    if(mod(TimeIn, DtLogfile)==0.0)then
       open(iUnitLog, FILE=NameFileLog, POSITION='APPEND')
       write(iUnitLog, '(es13.5,i5,5i3,1x,i3.3,10es13.5)') &
            TimeIn, &
            TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
            TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond, &
            floor(TimeRamNow%FracSecond*1000.0), dst, DstBiot, &
            sum(PParH)/(nR*nT),  sum(PPerH)/(nR*nT), &
            sum(PParO)/(nR*nT),  sum(PPerO)/(nR*nT), &
            sum(PParHe)/(nR*nT), sum(PPerHe)/(nR*nT), &
            sum(PParE)/(nR*nT), sum(PPerE)/(nR*nT)
       call flush(iUnitLog)
    end if

    ! Write SCB Dst file
    if(mod(TimeIn, Dt_hI)==0.0)then
       open(iUnitDst, FILE=NameFileDst, POSITION='APPEND')
       write(iUnitDst, '(es13.5,i5,5i3,1x,i3.3,6es13.5)') &
            TimeIn, &
            TimeRamNow%iYear, TimeRamNow%iMonth, TimeRamNow%iDay, &
            TimeRamNow%iHour, TimeRamNow%iMinute, TimeRamNow%iSecond, &
            floor(TimeRamNow%FracSecond*1000.0), dstDPS, dstDPSInsideGeo,&
            DstBiot, DstBiotInsideGeo, maxval(VT), minval(VT)
       call flush(iUnitDst)
    end if

    ! Check satellites.
    if(DoSaveRamSats .and. (mod(TimeIn,DtWriteSat) .eq. 0)) call fly_sats
    
    ! Write restarts.
    if(mod(TimeIn,DtRestart).eq.0) call write_restart()

  end subroutine handle_output

!=============================================================================
end subroutine run_ramscb

!=============================================================================
! Notes on conversion from Main.f to these subroutines:
!
! Variable iCal: Used to track iteration number in Main.f, now only needs
!                to indicate if this is the first iteration or not so that
!                RAM can initialize variables.  First set in IM_init_session.
!
! Variables UTs and TimeRamElapsed: UTs is the seconds from the starting day;
!                TimeRamElapsed is seconds from the start of the simulation.
!                Though essentially unnecessary, TimeRamElapsed is a 
!                convenience variable while UTs is the original time format
!                used in ram_all.f.  
!=============================================================================
