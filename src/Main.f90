program ram_scb

!============================================================================
!    RAM_SCB, Version 2.X
!
!    This is the main program unit for RAM_SCB; for compilation & execution,
!    please see accompanying documentation.
!
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

use ModRamMain,      ONLY: Real8_, S, iCal, nIter
use ModRamParams,    ONLY: DoSaveFinalRestart, DoVarDt, IsComponent, NameBoundMag
use ModRamGrids,     ONLY: NR, NT, NPA, NE
use ModRamTiming,    ONLY: DtsFramework, DtsMax, DtsMin, DtsNext, Dts, Dt_hI, &
                           TimeRamStart, TimeRamNow, TimeMax, TimeRamElapsed, &
                           TimeRamStop, T, UTs, Dt_bc, DtEfi
use ModRamVariables, ONLY: Kp, F107, PPARH, PPERH, PPARO, PPERO, PPARE, PPERE, &
                           PPARHE, PPERHE, PPART, PPERT, BNES, FNHS, FNIS, HDNS, &
                           BOUNHS, BOUNIS, dIdt, dIbndt, dBdt, EIR, EIP, VT, F2, &
                           FLUX, LZ, MU, Pabn, FFACTOR

use ModScbGrids,     ONLY: npsi, nzeta
use ModScbVariables, ONLY: h_Cart, I_Cart, bZEq_Cart, hdens_Cart, h_Cart_interp, &
                           I_Cart_interp


use ModRamFunctions, ONLY: FUNT, FUNI, COSD, ram_sum_pressure
use ModRamIndices,   ONLY: get_indices
use ModRamTiming,    ONLY: max_output_timestep, init_timing, finalize_timing, do_timing
use ModRamIO,        ONLY: init_output, handle_output, ram_write_pressure
use ModRamRestart,   ONLY: write_restart
use ModRamInit,      ONLY: ram_init, init_input
use ModRamRun,       ONLY: ram_run
use ModRamBoundary,  ONLY: get_boundary_flux
use ModRamEField,    ONLY: get_electric_field

use ModScbInit, ONLY: scb_init
use ModScbRun,  ONLY: scb_run

use ModRamScb,  ONLY: computehI

use ModReadParam
use CON_planet,      ONLY: set_planet_defaults
use CON_axes,        ONLY: init_axes, test_axes
use ModPlanetConst,  ONLY: init_planet_const
use ModTimeConvert,  ONLY: time_real_to_int

use ModMpi
use ModRamMpi

implicit none

integer :: I, J, L, K
real(kind=Real8_) :: BNESPrev(NR+1,NT), FNISPrev(NR+1,NT,NPA), BOUNISPrev(NR+1,NT,NPA)
real(kind=Real8_) :: DtOutputMax, DtEndMax
real(kind=Real8_) :: fluxVolume(npsi,nzeta)

character(len=4) :: StrScbIter

!----------------------------------------------------------------------------
! Ensure code is set to StandAlone mode.
IsComponent = .false.

!!!!!!!!!! START RAM-SCB INITIALIZATION !!!!!!!!!!!\
! Initialize MPI.
call MPI_INIT(ierror)
iComm = MPI_COMM_WORLD
call MPI_COMM_RANK(iComm,iProc,ierror)
call MPI_COMM_SIZE(iComm,nProc,ierror)

! Begin tracking code efficiency.
call init_timing()

! Read and set parameters
call read_file('PARAM.in', iComm)
call read_init('  ', iSessionIn=1, iLineIn=0)
call IM_set_parameters

! Initialize RAM_SCB
call ram_init
call scb_init


! Initialize planet and axes if in stand-alone mode.
! In component mode, the framework takes care of this.
! This is needed for correct coordinate transformations.
if (.not.IsComponent) then
   call init_planet_const
   call set_planet_defaults
   call init_axes(TimeRamStart % Time)
end if

! Get initial Data (checks for restart)
call init_input

! Initialize output files as necessary.
call init_output

iCal = 1
!!!!!!!!!!! END RAM-SCB INITIALIZATION !!!!!!!!!!!


!!!!!!!!!!! START RAM-SCB RUN !!!!!!!!!!!!
if (TimeRamElapsed .lt. TimeMax) then ! No wasted cycles, please.
   MAIN: do
!!!!!!!!!! SET CURRENT TIME AND TIMESTEP
      UTs = TimeRamElapsed + TimeRamStart%iHour*3600.0
      ! Set the half-timestep based on CFL, SWMF, and Maximum allowed.
      ! Use CFL Number if #VARIABLEDT is set in PARAM file.
      if (DoVarDt) then
         DtEndMax   = (TimeMax-TimeRamElapsed)/2.0
         DtOutputMax = max_output_timestep(TimeRamElapsed)
         DTs = min(DTsNext, DTsmax, DtsFramework, DtOutputMax, DtEndMax)
         if (Kp.gt.6.0 .AND. DTs.gt.1.) DTs = 1.   !1. or 0.25 
      else if(mod(UTs, Dt_hI) .eq. 0) then
         DTs = 5.0
         if(Kp .ge. 5.0) DTs = min(DTsMin,DTs)
         if(Kp .gt. 6.0) DTs = 1.0   !1. or 0.25
      endif
      DTsNext = DTsmax
!!!!!!!!!!

!!!!!!!!!! UPDATES AS NEEDED
      call get_indices(TimeRamNow%Time, Kp, f107)

      ! Update Boundary Flux if Dt_bc has passed
      if (mod(TimeRamElapsed, Dt_bc).eq.0) then
         call get_boundary_flux
      end if

      ! Update Electric Fields if DtEfi has passed
      if (mod(TimeRamElapsed, DtEfi).eq.0) then
         call get_electric_field
      end if
!!!!!!!!!

!!!!!!!!!! RUN RAM
      ! Broadcast current call to ram_all
      call write_prefix
      write(*,'(1x, a, 1x, F8.1, 3x, F5.1, 3x, F5.3)') &
              'Calling ram_all for UTs, DTs,Kp = ', UTs, DTs, Kp

      ! Call RAM for each species.
      call ram_run

      ! Update species pressures
      PPerO  = PPerT(4,:,:); PParO  = PParT(4,:,:)
      PPerHe = PPerT(3,:,:); PParHe = PParT(3,:,:)
      PPerE  = PPerT(1,:,:); PParE  = PParT(1,:,:)
      PPerH  = PPerT(2,:,:); PParH  = PParT(2,:,:)

      ! Update flux
      DO I = 2, NR
         DO K = 2, NE
            DO L = 2, NPA
               DO J = 1, NT-1
                  FLUX(S,I,J,K,L) = F2(S,I,J,K,L)/FFACTOR(S,I,K,L)/FNHS(I,J,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ! Increment and update time
      TimeRamElapsed = TimeRamElapsed + 2.0 * DTs
      TimeRamNow % Time = TimeRamStart % Time + TimeRamElapsed
      call time_real_to_int(TimeRamNow)
!      T = TimeRamElapsed + TimeRamStart%iHour*3600.0
!!!!!!!!!

!!!!!!!!! RUN SCB
      ! Call SCB if Dt_hI has passed
      if (mod(TimeRamElapsed, Dt_hI) .eq. 0) then
         call write_prefix
         write(*,*) 'Running SCB model to update B-field...'

         write(StrScbIter,'(I4.4)') int(TimeRamElapsed/Dt_hI)
         call ram_sum_pressure
!         call ram_write_pressure(StrScbIter)

         call scb_run(fluxVolume)

         ! Couple SCB -> RAM
         call computehI(fluxVolume)

         call write_prefix
         write(*,*) 'Finished 3D Equilibrium code.'

         ! Update h and I values if Dt_hI has passed
         IF (NameBoundMag .NE. 'DIPL') THEN
          DO I=2,NR+1
            DO J=1,NT
              BNESPrev(I,J)=BNES(I,J)
              DO L=1,NPA
                FNISPrev(I,J,L)=FNIS(I,J,L)
                BOUNISPrev(I,J,L)=BOUNIS(I,J,L)
                ! SCB Variables -> RAM Variables
                FNHS(I,J,L)   = h_Cart(I-1,J,L)
                FNIS(I,J,L)   = i_Cart(I-1,J,L)
                BNES(I,J)     = bZEq_Cart(I-1,J)
                HDNS(I,J,L)   = hdens_Cart(I-1,J,L)
                BOUNHS(I,J,L) = h_Cart_interp(I-1,J,L)
                BOUNIS(I,J,L) = i_Cart_interp(I-1,J,L)
                !
                dIdt(I,J,L)=(FNIS(I,J,L)-FNISPrev(I,J,L))/DT_hI
                dIbndt(I,J,L)=(BOUNIS(I,J,L)-BOUNISPrev(I,J,L))/DT_hI
              ENDDO
              IF (J.LT.8.or.J.GT.18) THEN
                DO L=15,2,-1
                  if (FNHS(I,J,L-1).gt.FNHS(I,J,L)) then
                    FNHS(I,J,L-1)=0.99*FNHS(I,J,L)
                  endif
                ENDDO
              ENDIF
              BNES(I,J)=BNES(I,J)/1e9 ! to convert in [T]
              dBdt(I,J)=(BNES(I,J)-BNESPrev(I,J))/DT_hI
            ENDDO
          ENDDO
          DO J=1,NT ! use dipole B at I=1
            BNES(1,J)=0.32/LZ(1)**3/1.e4
            dBdt(1,J) = 0.
            EIR(1,J) = 0.
            EIP(1,J) = 0.
            DO L=1,NPA
              FNHS(1,J,L) = FUNT(MU(L))
              FNIS(1,J,L) = FUNI(MU(L))
              BOUNHS(1,J,L)=FUNT(COSD(PAbn(L)))
              BOUNIS(1,J,L)=FUNI(COSD(PAbn(L)))
              HDNS(1,J,L)=HDNS(2,J,L)
              dIdt(1,J,L)=0.
              dIbndt(1,J,L)=0.
            ENDDO
          ENDDO
         ENDIF
      end if
!!!!!!!!!

!!!!!!!!! CREATE OUTPUTS
      ! Do timing.
      call do_timing
      ! Check and write output, as necessary.
      call handle_output(TimeRamElapsed)
!!!!!!!!!

!!!!!!!!! FINISH TIME STEP
      iCal  = iCal + 1
      nIter = nIter + 1
      if(TimeRamElapsed .ge. TimeMax) exit MAIN
!!!!!!!!!
   end do MAIN
end if

if ((.not.IsComponent)) then
     write(*,*) &
          '==============================================================='
     write(*,*) '                     RAM_SCB finished!'

     write(*,'(f12.2,a,i10.10,a)') &
          TimeRamElapsed, 's simulated in ', iCal-1, ' iterations.'

     if(DoSaveFinalRestart) call write_restart

     ! Stop timing Ram; write final timing report.
     call finalize_timing

     write(*,*) &
          '==============================================================='
end if


call MPI_Finalize(iError)

end program ram_scb
!==============================================================================
  subroutine write_prefix

    use ModRamParams, ONLY: IsComponent

    implicit none

    character(len=7) :: StringPrefix = 'IM:'

    if(.not. IsComponent) RETURN
    write(*,'(a)',ADVANCE='NO')trim(StringPrefix)
   
  end subroutine write_prefix

!============================================================================
subroutine CON_stop(String)
  ! "Safely" stop RAM-SCB on all nodes.
  use ModRamTiming, ONLY: TimeRamElapsed
  implicit none

  character(len=*), intent(in) :: String
  write(*,*)'Stopping execution! at time=',TimeRamElapsed,&
       ' with msg:'
  write(*,*)String
  stop
end subroutine CON_stop

!============================================================================

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  ! Replaces the SWMF testing routine.
  
  use ModRamParams, ONLY: StringTest

  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest   = i_sub_string(' '//StringTest,' '//String//' ')>0
  DoTestMe = DoTest !.and. i_proc()==iProcTest
contains
  !===========================================================================
  integer function i_sub_string(StringA,StringB)

    ! This is needed to avoid some SGI f90 compiler bug 
    ! (which results in a memory leak) if we use
    !
    ! index(' '//StringTest,' '//str//' ')
    !
    ! directly.

    implicit none

    character (len=*), intent(in) :: StringA, StringB

    i_sub_string=index(StringA, StringB)

  end function i_sub_string

end subroutine CON_set_do_test

!============================================================================
