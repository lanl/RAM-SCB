!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

Module ModRamInjection

  implicit none

  contains
!==============================================================================
  subroutine injection(iCal)

!    use ModRamMain,      ONLY: Real8_, S
!    use ModRamVariables, ONLY: F2
!    use ModRamGrids,     ONLY: nR, nT, nE, nPa
!    use ModRamRun,       ONLY: ANISCH, SUMRC
!    use ModRamFunctions, ONLY: ram_sum_pressure
!    use ModScbRun,       ONLY: scb_run
!    use ModRamScb,       ONLY: computehI

!    implicit none
!
    integer, intent(in) :: iCal
!    integer :: iS, I, J, K, L, R
!    real(kind=Real8_) :: fluxAdd, M
!
!    R = 1
!    M = 2.0
    if (iCal.eq.100) then
return
!       R = 1
!       M = 2.0
    elseif (iCal.eq.500) then
return
!       R = 3
!       M = 2.0
    elseif (iCal.eq.1200) then
return
!       R = 1
!       M = 1.0
    elseif (iCal.eq.2000) then
return
!       R = 2
!       M = 2.0
    endif
!
!    write(*,*) '!!!! INJECTION !!!!'
!    do iS=1,4
!       S = iS
!       do I=nR-6,nR
!          do K=1,nE
!             do L=1,nPa
!                fluxAdd = maxval(F2(S,nR-6:nR,nT-R:nT,:,L))
!                do J=nT-R,nT
!                   F2(S,I,J,K,L) = M*fluxAdd
!                enddo
!                do J=1,1+R
!                   F2(S,I,J,K,L) = M*fluxAdd
!                enddo
!             enddo
!          enddo
!       enddo
!       call SUMRC(S)
!       call ANISCH(S)
!    enddo
!    call ram_sum_pressure
!    call scb_run
!    call computehI

  end subroutine injection

!==============================================================================
End Module ModRamInjection
