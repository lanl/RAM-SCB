SUBROUTINE Interpolation_natgrid_2D(r_local, azim_local, hFlux_l, IFlux_l, bZEqFlux_l, &
     flux_vol_l, hdensFlux_l, hCart_l, ICart_l, bZEqCart_l, flux_vol_Cart_l, hdensCart_l, &
     Epot_l, EpotCart_l)

  ! gives values on polar grid, using interpolation

  ! .. Use Statements ..
  USE nrtype, ONLY : DP, pi_d
  USE Module1, ONLY : x, y, z, nThetaEquator, nZetaMidnight, rank, numProc,prefixOut
  USE mpi

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r_local(:), azim_local(:)
  REAL(DP), INTENT(IN) :: hFlux_l(:,:,:), IFlux_l(:,:,:), bZEqFlux_l(:,:), &
       flux_vol_l(:,:), hdensFlux_l(:,:,:)
  REAL(DP), OPTIONAL, INTENT(IN) :: Epot_l(:,:)
  REAL(DP), INTENT(OUT) :: hCart_l(:,:,:), ICart_l(:,:,:), bZEqCart_l(:,:), &
       hdensCart_l(:,:,:), flux_vol_Cart_l(:,:)
  REAL(DP), OPTIONAL, INTENT(OUT) :: EpotCart_l(:,:)
  INTEGER :: nxl, nyl, npsil, nzetal, npal

  INTEGER :: iTriangle(SIZE(r_local), SIZE(azim_local))
  INTEGER :: nTotal 

  REAL(DP) :: xMax, yMax, zMax, xMin, yMin, zMin, dx, dy, dz, bet, del, end_time
  INTEGER :: ix, ier, i, j, k, iTotal, iTotal2, ierr, idealerr
  INTEGER :: nq, nr, nw

  REAL(DP) :: eps, eq, eqx, eqy, eqz, p2y, p2z, q1, rmax, xtemp, pMin
  REAL(DP) :: triangle(2,3)
  INTEGER :: l

  ! .. Intrinsic Functions ..
  INTRINSIC KIND, REAL, SUM
  ! .. Parameters ..
  INTEGER, PARAMETER :: wp = KIND(1.0D0)
  ! .. Local Scalars ..
  INTEGER :: m, numberTriangles
  REAL (DP) :: q, qx, qy, qz, u, v, w, rad, radSqLoc
  REAL(DP) :: center(2)

  REAL(DP), ALLOCATABLE :: xScatter(:), yScatter(:), radSqScatter(:), bZScatter(:), EpotScatter(:)
  REAL(DP), ALLOCATABLE :: hScatter(:,:), IScatter(:,:), hdensScatter(:,:), fluxVolScatter(:)

  REAL(DP) :: coordPoint(2)

  INTEGER :: jy, it, isInTriangle,iDim, itMin, inBigSphere

  npsil = SIZE(hFlux_l,1)
  nzetal = SIZE(hFlux_l,2) + 1 ! Since only 2:nzeta is being passed
  npal = SIZE(hFlux_l,3)
  nxl = SIZE(hCart_l,1)
  nyl = SIZE(hCart_l,2)

  nTotal = npsil*(nzetal-1) + 1

  IF(.NOT.ALLOCATED(xScatter))   ALLOCATE(xScatter(nTotal), stat = ierr)
  IF(.NOT.ALLOCATED(yScatter))   ALLOCATE(yScatter(nTotal), stat = ierr)
  IF(.NOT.ALLOCATED(radSqScatter))   ALLOCATE(radSqScatter(nTotal), stat = ierr)

  IF(.NOT.ALLOCATED(hScatter))   ALLOCATE(hScatter(nTotal, NPAL), stat = ierr)
  IF(.NOT.ALLOCATED(hDensScatter))   ALLOCATE(hDensScatter(nTotal, NPAL), stat = ierr)
  IF(.NOT.ALLOCATED(IScatter))   ALLOCATE(IScatter(nTotal, NPAL), stat = ierr)
  IF(.NOT.ALLOCATED(bZScatter))   ALLOCATE(bZScatter(nTotal), stat = ierr)  
  IF (.NOT. ALLOCATED(fluxVolScatter)) ALLOCATE(fluxVolScatter(nTotal), stat = ierr)
  IF (PRESENT(Epot_l)) ALLOCATE(EpotScatter(nTotal), stat=ierr)

  iTotal = 0

  j = 1
  j_loop: DO WHILE (j <= npsil)     
     k = 2
     k_Loop:  DO WHILE(k < nzetal+1) 
        iTotal = iTotal+1
        xScatter(iTotal) = x(nThetaEquator,j,k)
        yScatter(iTotal) = y(nThetaEquator,j,k)
        radSqScatter(iTotal) = xScatter(iTotal)**2 + yScatter(iTotal)**2
        bZScatter(iTotal) = bZEqFlux_l(j,k-1)
        fluxVolScatter(iTotal) = flux_vol_l(j,k-1)
        IF (ALLOCATED(EpotScatter)) EPotScatter(iTotal) = Epot_l(j,k-1)
        Pitch_angle_loop: DO L = 1, NPAL
           hScatter(iTotal, L) = hFlux_l(j,k-1, L)
           IScatter(iTotal, L) = IFlux_l(j,k-1, L)
           hdensScatter(iTotal, L) = hdensFlux_l(j,k-1, L)
        END DO Pitch_angle_loop
        k = k+1
     END DO k_Loop
     j = j+1
  END DO j_loop

  ! Add a point in the center for proper triangulation
  ! iTotal = iTotal + 1
  ! xScatter(iTotal) = 0._dp
  ! yScatter(iTotal) = 0._dp
  ! radSqScatter(iTotal) = 0._dp
  ! bZScatter(iTotal) = bZEqFlux_l(1,nZetaMidnight) 
  ! DO L = 1, NPAL
  !   hScatter(iTotal, L) = hFlux_l(1,nZetaMidnight,L)
  !   IScatter(iTotal, L) = IFlux_l(1,nZetaMidnight,L)
  ! end DO

  CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), bZScatter(1:iTotal))
  DO ix = 1, nxl
     DO jy = 1, nyl
        coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
        coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
        CALL NNPNTD(coordPoint(1),coordPoint(2),bZEqCart_l(ix,jy))
     END DO
  END DO
  CALL NNPNTENDD()

  CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), fluxVolScatter(1:iTotal))
  DO ix = 1, nxl
     DO jy = 1, nyl
        coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
        coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
        CALL NNPNTD(coordPoint(1),coordPoint(2),flux_vol_Cart_l(ix,jy))
     END DO
  END DO
  CALL NNPNTENDD()

  IF (PRESENT(Epot_l)) THEN
     CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), EpotScatter(1:iTotal))
     DO ix = 1, nxl
        DO jy = 1, nyl
           coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           CALL NNPNTD(coordPoint(1),coordPoint(2),EpotCart_l(ix,jy))
        END DO
     END DO
     CALL NNPNTENDD()
  END IF

  DO L = 1, NPAL
     CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), hScatter(1:iTotal,L))
     DO ix = 1, nxl
        DO jy = 1, nyl
           coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           CALL NNPNTD(coordPoint(1),coordPoint(2),hCart_l(ix,jy,L))
        END DO
     END DO
     CALL NNPNTENDD()

     CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), hdensScatter(1:iTotal,L))
     DO ix = 1, nxl
        DO jy = 1, nyl
           coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           CALL NNPNTD(coordPoint(1),coordPoint(2),hdensCart_l(ix,jy,L))
        END DO
     END DO
     CALL NNPNTENDD()

     CALL NNPNTINITD(iTotal, xScatter(1:iTotal), yScatter(1:iTotal), IScatter(1:iTotal,L))
     DO ix = 1, nxl
        DO jy = 1, nyl
           coordPoint(1) = r_local(ix) * COS(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           coordPoint(2) = r_local(ix) * SIN(azim_local(jy)*2._dp*pi_d/24._dp - pi_d)
           CALL NNPNTD(coordPoint(1),coordPoint(2),ICart_l(ix,jy,L))
        END DO
     END DO
     CALL NNPNTENDD()

  END DO

  IF(ALLOCATED(xScatter))  DEALLOCATE(xScatter, stat = ierr)
  IF (ALLOCATED(yScatter))  DEALLOCATE(yScatter, stat = ierr)
  IF(ALLOCATED(radSqScatter))  DEALLOCATE(radSqScatter, stat = ierr)

  IF(ALLOCATED(hScatter))  DEALLOCATE(hScatter, stat = ierr)
  IF(ALLOCATED(IScatter))  DEALLOCATE(IScatter, stat = ierr)
  IF(ALLOCATED(bZScatter))  DEALLOCATE(bZScatter, stat = ierr)
  IF (ALLOCATED(EpotScatter)) DEALLOCATE(EpotScatter, stat = ierr)
  IF(ALLOCATED(hdensScatter))  DEALLOCATE(hdensScatter, stat = ierr)

  IF (ALLOCATED(fluxVolScatter)) DEALLOCATE(fluxVolScatter, stat = ierr)

  RETURN

END SUBROUTINE Interpolation_natgrid_2D


