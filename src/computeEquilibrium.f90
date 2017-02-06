SUBROUTINE computeEquilibrium(ST3)

  !****************************************************************************
  !  This is the main iteration routine for equilibrium solution
  !  finds flux coordinates inside plasma for a given flux boundary
  !  input: definition of boundary and pressure distribution

  !  Author: S. Zaharia
  !  Copyright (c) 2016, Los Alamos National Security, LLC
  !  All rights reserved.
  !****************************************************************************


  USE nrtype
  USE Module1
  USE mpi
  USE ModIoUnit, ONLY: UNITTMP_, io_unit_new
  USE ModRamMpi, ONLY: iComm
  IMPLICIT NONE

  CHARACTER*4  :: ST3  ! Vania

  INTERFACE RHS_alpha
     SUBROUTINE newk(vec_x, vec_dk)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       INTEGER         :: i, j, k 
       REAL(DP), DIMENSION(:, :, :), INTENT(OUT) :: vec_x
       REAL(DP), DIMENSION(:, :, :), INTENT(IN OUT) ::  vec_dk
     END SUBROUTINE newk
  END INTERFACE

  INTERFACE RHS_psi
     SUBROUTINE newj(vec_r, vec_dj)
       USE Module1
       USE nrtype
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: vec_r
       REAL(DP), DIMENSION(:,:,:), INTENT(IN OUT) :: vec_dj
     END SUBROUTINE newj
  END INTERFACE

  INTERFACE LHS_alpha
     SUBROUTINE metrica(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       REAL(DP), DIMENSION(:, :, :), INTENT(OUT) :: vecd, vec1, vec2, vec3, &
            vec4, vec6, vec7, vec8, vec9
     END SUBROUTINE metrica
  END INTERFACE

  INTERFACE LHS_psi
     SUBROUTINE metric(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: vecd, vec1, vec2, vec3, &
            vec4, vec6, vec7, vec8, vec9
     END SUBROUTINE metric
  END INTERFACE

  INTERFACE solutionEquilibrium_alpha
     SUBROUTINE iterateAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vec1, vec2, vec3, vec4, vec6, vec7, &
            vec8, vec9, vecx
       REAL(DP), DIMENSION(:,:,:), INTENT(IN OUT) :: vecd
     END SUBROUTINE iterateAlpha
  END INTERFACE
  INTERFACE solutionEquilibrium_dalpha
     SUBROUTINE directAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)
       USE nrtype
       USE Module1 
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd, vec1, vec2, vec3, vec4, vec6, vec7, &
            vec8, vec9, vecx
     END SUBROUTINE directAlpha
  END INTERFACE
  INTERFACE solutionEquilibrium_psi
     SUBROUTINE iteratePsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd,vec1,vec2,vec3,vec4,vec6, &
            vec7,vec8,vec9,vecr
     END SUBROUTINE iteratePsi
  END INTERFACE
  INTERFACE solutionEquilibrium_dpsi  
     SUBROUTINE directPsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)
       USE nrtype
       USE Module1
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: vecd,vec1,vec2,vec3,vec4, &
            & vec6,vec7,vec8,vec9,vecr
     END SUBROUTINE directPsi
  END INTERFACE

  INTERFACE pressure
     SUBROUTINE pressure(entropy_local, vol_local, icount_local)
       USE nrtype
       USE Module1
       REAL(DP), INTENT(IN) :: entropy_local(:,:), vol_local(:,:)
       INTEGER, INTENT(IN) :: icount_local
     END SUBROUTINE pressure
  END INTERFACE

  INTERFACE entropy
     SUBROUTINE entropy(ent_local, vol_local, iteration_local)
       USE nrtype
       USE Module1
       REAL(DP) :: ent_local(:,:), vol_local(:,:)
       INTEGER, INTENT (IN) :: iteration_local
     END  SUBROUTINE entropy
  END INTERFACE

  INTERFACE hRAM
     SUBROUTINE hRAM(ST3, flux_volume)
       USE nrtype, ONLY : DP, pi_d
       IMPLICIT NONE
       CHARACTER*4, INTENT(IN) :: ST3
       REAL(DP), INTENT(IN) :: flux_volume(:,:)
     END SUBROUTINE hRAM
     END INTERFACE

  INTEGER, PARAMETER :: NPA = 72
  INTEGER :: iconv, nisave1, ierr, iCountEntropy
  INTEGER :: i
  REAL(DP) :: sumdbconv, errorfirstalpha, diffmxfirstalpha, &
       errorfirstpsi, diffmxfirstpsi
  REAL(DP) :: sumb1, sumdb1, diffmx1
  REAL(DP) :: entropyFixed(npsi,nzeta), fluxVolume(npsi,nzeta) 
  REAL(DP)  :: hValue(npsi, nzeta, NPA)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: vecd, vec1, vec2, vec3, vec4, vec6, &
       vec7, vec8, vec9, vecx, vecr
  REAL(DP), ALLOCATABLE, SAVE :: xPrev(:,:,:), yPrev(:,:,:), zPrev(:,:,:)

  sumb1 = 0._dp
  sumdb1 = 0._dp
  diffmx1 = 0._dp
  entropyFixed = 0._dp
  fluxVolume = 0._dp

  ALLOCATE(vecd(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec1(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec2(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec3(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec4(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec6(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec7(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec8(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vec9(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vecx(nthe,npsi,nzeta), stat = ierr)
  ALLOCATE(vecr(nthe,npsi,nzeta), stat = ierr)

  IF (ierr /= 0) STOP 'Allocation problem in computeEquilibrium.'

  IF (rank == 0) THEN ! Developed only for isotropic pressure so far
     iUnitLog = io_unit_new()
     OPEN(iUnitLog, file = 'convergence_Details', status = 'replace')
     WRITE(iUnitLog, 9876); CALL FLUSH(iUnitLog)
     WRITE(*, 9876); CALL FLUSH(iUnitLog)
     IF (isFBDetailNeeded == 1) THEN
        WRITE(iUnitLog, *) 'Norms of |jxB-grad P|, |jxB|, |gradP| '; CALL FLUSH(iUnitLog)
        WRITE(*,*) 'Norms of |jxB-grad P|, |jxB|, |gradP| '; CALL FLUSH(6)  
     END IF
  END IF

  iteration = 1
  iCountEntropy = 1
  iconv = 0
  sumdbconv = 0.0_dp
  iConvGlobal = 0

  IF (iAMR == 1) THEN 
     CALL findR
     CALL InterpolatePsiR
     CALL mappsi(0) ! Full mapping needed, psis changed
     CALL psifunctions
     CALL maptheta
  ENDIF

  Outeriters: DO 

     CALL computeBandJacob_initial
     CALL metrica(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)


     !   define the right-hand side of the betaEuler equation
     CALL pressure(entropyFixed, fluxVolume, iCountEntropy)   ! Entropy etc. for some isotropic press. calculations, not used here
     
     IF (method /= 3) CALL newk(vecx, vecd)

     IF (isotropy /= 1 .AND. rank == 0 .AND. iteration == 1 .AND. isFBDetailNeeded == 1) THEN
        CALL test_Convergence_anisotropic 
        ! to see how far from equilibrium we are before computing
     END IF

     IF (isEnergDetailNeeded == 1) THEN
        SELECT CASE(isotropy)
        CASE(1) ! Isotropic
           CALL energy
        CASE default ! Anisotropic
           !C   CALL dsp_general ! Computes energies and Dst from DPS relation
        END SELECT
     END IF

     CALL MPI_BARRIER(iComm, ierr)


     NotTurPressure: IF (iPressureChoice /= 10) THEN ! Tur comparison, magnetic flux computation only
        errorAlphaPrev = errorAlpha
        
        IF (iteration==1) THEN
           alfaSav1(:,:,:) = alfa(:,:,:)
           alfaSav2(:,:,:) = alfa(:,:,:) ! Before the calculation
        END IF

        blendAlpha = MAX(blendAlpha,blendMin)
        blendAlpha = MIN(blendAlpha,blendMax)
        SELECT CASE (method)
        CASE(1)
           CALL directAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)
        CASE(2)
           CALL iterateAlpha(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecx)
        CASE(3)
           PRINT*, 'No calculation.'
           GO TO 1000
        END SELECT

        sumb1 = sumb
        sumdb1 = sumdb
        diffmx1 = diffmx
        nisave1 = nisave

        errorAlpha = diffmx

        IF (iteration == 1) THEN
           errorfirstalpha = sumdb1
           diffmxfirstalpha = diffmx1
           errorAlphaPrev = errorAlpha
        END IF
        IF (sumdb1 < sumdbconv) sumdbconv = sumdb1 

        iAlphaMove = 1

        IF (errorAlpha/errorAlphaPrev > thresh) THEN
           PRINT*, 'CE: decrease blendAlpha at iteration ', iteration
           blendAlpha = damp * blendAlpha
           IF (MOD(iteration,2)==1) THEN
              alfa(:,:,:) = alfaSav1(:,:,:)
           ELSE
              alfa(:,:,:) = alfaSav2(:,:,:)
           END IF
        END IF

        IF (method /= 3) THEN
           IF (.NOT. ALLOCATED(xPrev)) ALLOCATE(xPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
           IF (.NOT. ALLOCATED(yPrev)) ALLOCATE(yPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
           IF (.NOT. ALLOCATED(zPrev)) ALLOCATE(zPrev(SIZE(x,1), SIZE(x,2), SIZE(x,3)), STAT = ierr)
           xPrev = x
           yPrev = y
           zPrev = z

           Move_points_in_alpha_theta: DO
              !cc  move zeta grid points along constant alphaEuler and theta lines
              CALL mapalpha(iSm)
              !cc  move theta grid points along constant alphaEuler and zeta lines
              CALL maptheta
              CALL metric(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)
              IF (MINVAL(jacobian) < 0._dp) THEN
                 !C             iAlphaMove = iAlphaMove+1
                 blendAlpha = damp * blendAlpha
                 PRINT*, 'CE: Cycling alpha_theta pts, blendAlpha = ', blendAlpha
                 IF (MOD(iteration,2)==1) THEN
                    alfa(:,:,:) = alfa(:,:,:)*blendAlpha + (1.-blendAlpha)*alfaSav1(:,:,:)
                 ELSE
                    alfa(:,:,:) = alfa(:,:,:)*blendAlpha + (1.-blendAlpha)*alfaSav2(:,:,:)
                 END IF
                 ! Revert to previous point configuration
                 x = xPrev
                 y = yPrev
                 z = zPrev
                 CYCLE Move_points_in_alpha_theta
              END IF
              EXIT Move_points_in_alpha_theta
           END DO Move_points_in_alpha_theta
        END IF

        IF (MOD(iteration,2)==1) THEN
           alfaSav1(:,:,:) = alfa(:,:,:)
        ELSE
           alfaSav2(:,:,:) = alfa(:,:,:)
        END IF

        IF (iAMR == 1) THEN
           CALL findR
           CALL InterpolatePsiR
           CALL mappsi(0)  ! Full mapping needed, changed psis
           CALL psifunctions
           CALL maptheta
        ENDIF

     END IF NotTurPressure
     
     IF (MOD(iteration,nrelax) == 0) blendAlpha = relax * blendAlpha
     ! PRINT*, 'CE: blendAlpha = ', blendAlpha

     !     print*, 'CE: ST3 = ', ST3
     !     if (ST3.EQ.'0280') STOP

1001 CONTINUE

     IF (isotropy == 1) CALL entropy(entropyFixed, fluxVolume, iCountEntropy)

     CALL computeBandJacob_initial
     CALL metric(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)
     IF (MINVAL(jacobian) < 0._dp) STOP 'CE: metric problem.'
     CALL pressure(entropyFixed, fluxVolume, iCountEntropy)  

     !c  define the right-hand side of the alphaEuler equation
     IF (method /= 3) CALL newj(vecr, vecd)

     errorPsiPrev = errorPsi
     IF (iteration==1) THEN
        psiSav1(:,:,:) = psi(:,:,:)
        psiSav2(:,:,:) = psi(:,:,:) ! Before the calculation
     END IF

     blendPsi = MAX(blendPsi,blendMin)
     blendPsi = MIN(blendPsi,blendMax)
     SELECT CASE (method)
     CASE(1)
        CALL directPsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)
     CASE(2)
        CALL iteratePsi(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vecr)
     CASE(3)
        PRINT*, 'No iteration'
        GO TO 2000
     CASE default
        PRINT*, 'Wrong choice.'
        STOP
     END SELECT

     errorPsi = diffmx
     IF (iteration==1) errorPsiPrev = errorPsi

     IF (errorPsi/errorPsiPrev > thresh) THEN
        PRINT*, 'CE: decrease blendPsi at iteration ', iteration
        blendPsi = damp * blendPsi
        IF (MOD(iteration,2)==1) THEN
           psi(:,:,:) = psiSav1(:,:,:)
        ELSE
           psi(:,:,:) = psiSav2(:,:,:)
        END IF
     END IF

     ! if (rank == 0) print*, 'eralpha/er1stalpha, erpsi/er1stpsi = ', sumdb1/errorfirstalpha, sumdb/errorfirstpsi

     IF (rank == 0 .AND. method /= 3 .AND. errorAlpha/twopi_d < decreaseConvAlpha &
          .AND. errorPsi/MAXVAL(ABS(psival)) < decreaseConvPsi .AND. iconv == 0) iConvGlobal=1

     ! Broadcast iConvGlobal to all processors
     CALL MPI_BARRIER(iComm, ierr)
     CALL MPI_BCAST(iConvGlobal, 1, MPI_INTEGER, 0, iComm, ierr)
     CALL MPI_BARRIER(iComm, ierr)


     iPsiMove = 1

     IF (method /= 3) THEN
        xPrev = x
        yPrev = y
        zPrev = z
        Move_points_in_psi_theta: DO
           ! Call metrica to find out if the jacobian is well behaved
           !          PRINT*, 'iPsiMove, Jacobian after moving points (after psi equation):', iPsiMove
           CALL mappsi(iSm)
           CALL maptheta
           CALL metric(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9)
           IF (MINVAL(jacobian) < 0._dp) THEN
             !C iPsiMove = iPsiMove+1
              blendPsi = damp * blendPsi
              PRINT*, 'CE: Cycling psi_theta pts, blendPsi = ', blendPsi
              IF (MOD(iteration,2)==1) THEN
                 psi(:,:,:) = psi(:,:,:)*blendPsi + (1.-blendPsi)*psiSav1(:,:,:)
              ELSE
                 psi(:,:,:) = psi(:,:,:)*blendPsi + (1.-blendPsi)*psiSav2(:,:,:)
              END IF
              ! Revert to previous point configuration
              x = xPrev
              y = yPrev
              z = zPrev
              CYCLE Move_points_in_psi_theta
           END IF
           EXIT Move_points_in_psi_theta
        END DO Move_points_in_psi_theta
     END IF

  IF (MOD(iteration,2)==1) THEN
        psiSav1(:,:,:) = psi(:,:,:)
     ELSE
        psiSav2(:,:,:) = psi(:,:,:)
     END IF

     IF (MOD(iteration,nrelax) == 0) blendPsi = relax * blendPsi
     ! PRINT*, 'CE: blendPsi = ', blendPsi


     IF (rank == 0) THEN
        IF (iteration == 1) THEN
           errorfirstpsi = sumdb
           ! print*, 'error1stpsi = ', errorfirstpsi
           diffmxfirstpsi = diffmx
        END IF
9876    FORMAT(1x,'itout',1x,'blendAlpha', 1x, 'blendPsi', 1x, 'itAlpha',3x,'diffAlpha',4x,'errorAlpha' &
             &  ,3x,'itPsi',4x,'diffPsi',4x,'errorPsi')
        WRITE(*,9875) iteration, blendAlpha, blendPsi, nisave1,sumdb1,errorAlpha/twopi_d,nisave,sumdb,errorPsi/MAXVAL(ABS(psival))  !C relative diffmx errors now
        CALL FLUSH(6)
9875    FORMAT(1x,i4,1x, F5.2, 1x, F5.2, 4x,i5, 3X, E11.3, 2X, E11.3, 4X, i4,1x,3X,2e11.3)
        WRITE(iUnitLog, 9875) iteration, blendAlpha, blendPsi, nisave1,sumdb1,diffmx1/twopi_d*REAL(nzeta,dp),nisave,sumdb,diffmx
        CALL FLUSH(iUnitLog)  !C relative diffmx errors now
     END IF
     
     IF (iAMR == 1) THEN
        CALL findR
        CALL InterpolatePsiR
        CALL mappsi(0)  ! Full mapping needed, changed psis
        CALL psiFunctions
        CALL maptheta
     ENDIF

     IF (iConvGlobal == 1) THEN
        IF (rank == 0) PRINT*, 'Approaching convergence.'; CALL FLUSH(6)
        sumdbconv = sumdb1
        CALL MPI_BARRIER(iComm, ierr)
        EXIT Outeriters
     END IF

5012 FORMAT(1E12.4)

2000 CONTINUE

     CALL MPI_BARRIER(iComm, ierr)

     IF (isotropy /= 1 .AND. rank == 0 .AND. isFBDetailNeeded == 1) CALL test_Convergence_anisotropic
     CALL MPI_BARRIER(iComm, ierr)

     IF (iteration < numit .AND. iConvGlobal == 0) THEN
        iteration = iteration + 1
        CYCLE Outeriters
     END IF

     EXIT Outeriters
  END DO Outeriters

  iConv = 1
  iConvGlobal = 1

  !   The end of the iterative calculation
  IF (rank == 0) THEN
     PRINT*, iteration, " outer iterations performed."
     PRINT*, ' '
     CALL FLUSH(6)
  END IF

1000 IF (iPressureChoice /= 5) PRINT*, "End of calculation."

  IF (iteration > numit) lconv = 1
  nitry = nisave

  IF (numit == 0) THEN  ! In case no iteration was performed; still need pressure for metrics
     CALL pressure(entropyFixed, fluxVolume, iCountEntropy) 
     IF (isotropy == 1) CALL entropy(entropyFixed, fluxVolume, iteration)
     IF ((isotropy == 1) .AND. (isEnergDetailNeeded == 1)) THEN
        CALL energy
     END IF
  END IF

  ! The following block should be uncommented for applications where equal-arc-length is needed 
  ! or desirable
  !C  constTheta = 0.0_dp ! For equal-arc length
  !C  chiVal = (thetaVal + constTheta * SIN(2.*thetaVal)) 
  !C  CALL maptheta
  !C  CALL metrica(vecd,vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9) 
  ! This will compute the new Bfield on the new grid, to get the right pressure mapping

  CALL pressure(entropyFixed, fluxVolume, iCountEntropy)  ! To write to disk
  CALL MPI_BARRIER(iComm, ierr)

  CALL entropy(entropyFixed, fluxVolume, iCountEntropy)  

  IF (iElectric == 1 .OR. iElectric == 2 .or. iElectric == 3) THEN
     CALL ionospheric_potential
     PRINT*, '3DEQ: mapping iono. potentials along SCB-field lines'
  END IF
  
 
  ! Compute physical quantities: currents, field components etc..
  CALL metrics
  ! extrapolate to the fixed boundary
  IF (rank == 0) CALL bounextp 

   IF(isotropy == 0 .AND. isEnergDetailNeeded == 1 .AND. rank == 0) &
       CALL dps_general ! Computes energies and Dst from DPS relation, write to disk (+ Biot-Savart values) ! Remove for speed 

  CALL MPI_BARRIER(iComm, ierr)

  !C  if (rank == 0) print*, 'Calling plot_physical...'; call flush(6)
  IF (rank == 0) CALL plot_physical   ! Writes physical quantities to disk

  IF (isotropy == 0 .AND. iOutput /= 0 .AND. iPressureChoice /= 8) CALL hRAM(ST3, fluxVolume) ! Coupling with RAM
  CALL MPI_BARRIER(iComm, ierr)


  ! Deallocate 
  DEALLOCATE(vecd, stat = ierr)
  DEALLOCATE(vec1, stat = ierr)
  DEALLOCATE(vec2, stat = ierr)
  DEALLOCATE(vec3, stat = ierr)
  DEALLOCATE(vec4, stat = ierr)
  DEALLOCATE(vec6, stat = ierr)
  DEALLOCATE(vec7, stat = ierr)
  DEALLOCATE(vec8, stat = ierr)
  DEALLOCATE(vec9, stat = ierr)
  DEALLOCATE(vecx, stat = ierr)
  DEALLOCATE(vecr, stat = ierr) 

  end_time = MPI_WTIME( )

  CLOSE(iUnitLog)

  RETURN

END SUBROUTINE computeEquilibrium




