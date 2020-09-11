!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbIO
  
  use nrtype, ONLY: DP

  implicit none
 
  REAL(DP), SAVE :: PARMOD(10)
  INTEGER, SAVE  :: IOPT

  contains

!=============================================================================!
!============================= INPUT ROUTINES ================================!
!=============================================================================!

  subroutine computational_domain
    !!!! Module Variables  
    USE ModRamVariables, ONLY: Kp, MLT
    use ModRamParams,    ONLY: IsComponent, NameBoundMag, verbose
    use ModRamTiming,    ONLY: TimeRamNow
    use ModRamGrids,     ONLY: nRExtend, nT, nR, RadiusMax, RadiusMin
    USE ModScbParams,    ONLY: Symmetric, constZ, constTheta
    USE ModScbGrids,     ONLY: npsi, nthe, nzeta
    USE ModScbVariables, ONLY: xpsiin, xpsiout, r0Start, bf, nThetaEquator, &
                               x, y, z, thetaVal, zetaVal, left, right, chiVal, &
                               kmax, nZetaMidnight, xzero3, f, fzet, fp, &
                               fzetp, psiVal, alphaVal, psiin, psiout, psitot, SORFail
    !!!! SMWF Variables
    use ModRamCouple,    ONLY: generate_field
    !!!! Module Subroutines/Functions
    use ModRamGSL,       ONLY: GSL_Interpolation_1D
    use ModRamFunctions, ONLY: RamFileName
    use ModScbEuler,     ONLY: psiges, alfges, InterpolatePsiR, mappsi, maptheta, &
                               psifunctions
    use ModScbFunctions, ONLY: get_dipole_lines
    !!!! Share Modules
    use ModTimeConvert, ONLY: n_day_of_year, time_int_to_real
    USE ModIoUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, pi_d, twopi_d

    implicit none

    character(len=200) :: FileName

    logical :: check, c(24)
    INTEGER :: i, ii, j, k, scanLeft, scanRight, GSLerr
    integer :: Year, Month, Day, Hour, Mins, sec
    integer :: iMinute, iSecond

    REAL(DP) :: dphi, phi, psis, xpsitot, xpl, tempVal
    REAL(DP) :: r0, t0, t1, tt, zt, b, rr, rt, psitemp, rc(1)
    REAL(DP) :: Pdyn, Dst, ByIMF, BzIMF, G(3), W(6)
    REAL(DP), DIMENSION(1000) :: distance, xx, yy, zz, bx, by, bz
    INTEGER :: LMAX = 1000
    INTEGER :: LOUT, iYear, iMonth, iDay, iHour, iMin, iSec, ifail
    REAL(DP) :: x0, y0, z0, xf, yf, zf, rf, tf, b0

    integer :: time1, clock_rate, clock_max
    real(dp) :: starttime,stoptime
    clock_rate = 1000
    clock_max = 100000

    left = 1
    right = npsi
    r0Start = 1.0
    if ((NameBoundMag.eq.'DIPL').or.(NameBoundMag.eq.'DIPS').or.(NameBoundMag.eq.'DIPC')) then
       ! For generating x, y, and z arrays using analytic dipole and analytic compressed dipole
       ! the variable b controls the compression with 0 being no compression
       Symmetric = .true.
       xpsiin = 1.75
       xpsiout = 7.50
       call get_dipole_lines(xpsiin,xpsiout,constTheta,nthe,npsi,nzeta,x,y,z,bf,.false.)
       chival = (thetaVal + constTheta*sin(2.*thetaVal))
       kmax = nZetaMidnight
       x(:,:,nzeta+1) = x(:,:,2)
       y(:,:,nzeta+1) = y(:,:,2)
       z(:,:,nzeta+1) = z(:,:,2)

    elseif (NameBoundMag.eq.'SWMF') then
       ! For generating x, y, and z arrays using the Space Weather Modelling Framework
       call generate_field(xpsiin,xpsiout)
       if (SORFail) return
    else
       ! For generating x, y, and z arrays using field line tracing
       !! Define the inputs needed for the magnetic field models for the tracing
       Symmetric = .false. 

       ! Get correct model inputs and place them in cooresponding variables
       call get_model_inputs(Pdyn,Dst,ByIMF,BzIMF,G,W)

       ! Start tracing timing
       write(*,'(1x,a)',ADVANCE='NO') NameBoundMag//' tracing starting'
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)
   
       ! Find the correct starting point for the outer edge.
       ! We have to scan the night sector for the most stretched
       ! location since it isn't always at midnight
       xpsiout = 8.0
       scanLeft  = nZetaMidnight/2
       scanRight = nZetaMidnight + scanLeft
       r0 = 0.0
       Find_Outer: Do
          do k = scanLeft,scanRight
             x0 = (8._dp-r0)*cos(zetaVal(k))
             y0 = (8._dp-r0)*sin(zetaVal(k))
             z0 = 0._dp
             call trace(x0,y0,z0,1.0_dp,xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
             psitemp = 1./(1.-zf**2/(xf**2+yf**2+zf**2))
             if (psitemp.lt.xpsiout) then
                xpsiout = psitemp
                kmax = k
             endif
          enddo
          if (xpsiout.gt.3.0) then
             exit Find_Outer
          else
             r0 = r0 + 0.25
             xpsiout = 8.0-r0
          endif
          if (r0.gt.3.0) exit Find_Outer
       enddo Find_Outer
       xpsiin = 1.75 ! Don't need to do tracing at inner boundary

       ! Calculate dipole starting points for given xpsiin and xpsiout
       ! in chosen field and perform nzeta*npsi traces to create grid
       do k=2,nzeta
          do j=1,npsi
             r0 = xpsiin + REAL(j-1,DP)/REAL(npsi-1,DP)*(xpsiout-xpsiin)
             tt = pi_d-asin(sqrt(1.0/r0))
             rt = r0*sin(tt)**2
             zt = zetaVal(k)+constZ*sin(zetaVal(k))
             x0 = rt*cos(zt)*sin(tt)
             y0 = rt*sin(zt)*sin(tt)
             z0 = rt*cos(tt)
             CALL trace(x0,y0,z0,-1.0_dp,xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
             distance(1) = 0._dp
             do i = 2,LOUT
                distance(i) = distance(i-1) + SQRT((xx(i)-xx(i-1))**2 &
                              +(yy(i)-yy(i-1))**2 +(zz(i)-zz(i-1))**2)
             enddo
             chiVal = (thetaVal + constTheta * SIN(2.*thetaVal)) * distance(LOUT)/pi_d

             CALL GSL_Interpolation_1D(distance(1:LOUT),xx(1:LOUT),chiVal(2:nthe),x(2:nthe,j,k),GSLerr)
             CALL GSL_Interpolation_1D(distance(1:LOUT),yy(1:LOUT),chiVal(2:nthe),y(2:nthe,j,k),GSLerr)
             CALL GSL_Interpolation_1D(distance(1:LOUT),zz(1:LOUT),chiVal(2:nthe),z(2:nthe,j,k),GSLerr)
             x(1,j,k) = x0
             y(1,j,k) = y0
             z(1,j,k) = z0
             x(nthe,j,k) = x0
             y(nthe,j,k) = y0
             z(nthe,j,k) = -z0
          enddo
       enddo

       ! Finish tracing timing
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       write(*,'(a,1x,f6.2,1x,a)') ': Completed in',stoptime-starttime,'seconds'

       ! Periodic in zeta
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzeta+1) = x(:,:,2)
       y(:,:,nzeta+1) = y(:,:,2)
       z(:,:,nzeta+1) = z(:,:,2)
       chival = (thetaVal + constTheta*sin(2.*thetaVal))
    endif

    ! Output magnetic field for testing
    !call Write_MAGxyz('ComputeDomain')

    ! Get the Psi (Alpha) Euler Potential
    ! This is done by assuming a dipole on the field line foot points
    ! and then assigning the value of where they would cross the equator
    ! to the actual equatorial cross point
    psiin   = -xzero3/xpsiin
    psiout  = -xzero3/xpsiout
    psitot  = psiout-psiin
    xpsitot = xpsiout - xpsiin
    DO j = 1, npsi
       psis = REAL(j-1, DP) / REAL(npsi-1, DP)
       xpl = xpsiin + xpsitot * psis
       psival(j) = -xzero3 / xpl
       f(j) = (xzero3 / xpl**2) * xpsitot !dPsi/dR -- For converting between euler potential
       fp(j) = 0._dp                      !           and the curvilinear coordinate
    END DO
    call psiges

    dphi  = twopi_d/REAL(nzeta-1, DP)
    DO k = 1, nzeta+1
       phi         = REAL(k-2, DP) * dphi
       alphaVal(k) = phi + constZ*sin(phi)
       fzet(k)     = 1._dp ! dAlpha/dPhi -- For converting between the euler potential
       fzetp(k)    = 0._dp !                and the curvilinar coordinate
    END DO
    call alfges

    return

  end subroutine computational_domain

!=============================================================================!
  subroutine update_domain(updated)

    !!! Module Variables
    use ModRamParams,    ONLY: NameBoundMag, verbose
    use ModRamTiming,    ONLY: TimeRamNow
    use ModRamVariables, ONLY: Kp
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbParams,    ONLY: constTheta, constZ
    use ModScbVariables, ONLY: x, y, z, psiVal, alphaVal, psi, psiin, psiout, &
                               psitot, xpsiin, xpsiout, f, fp, nThetaEquator, &
                               fzet, fzetp, thetaVal, xzero3, kmax, zetaVal, &
                               nZetaMidnight, chiVal, SORFail
    !!! Module Subroutines/Functions
    use ModRamGSL,       ONLY: GSL_Interpolation_1D
    use ModRamFunctions, ONLY: RamFileName
    use ModScbEuler,     ONLY: mapAlpha, mapPsi, InterpolatePsiR, mapTheta, &
                               psiges, alfges, psifunctions
    use ModScbFunctions, ONLY: extap
    use ModScbCompute,   ONLY: ComputeBandJacob
    !!!! Share Modules
    use ModTimeConvert, ONLY: n_day_of_year
    USE ModIOUnit,      ONLY: UNITTMP_
    !!! NR Modules
    use nrtype, ONLY: DP, pi_d, twopi_d

    implicit none

    LOGICAL, INTENT(OUT) :: updated
    LOGICAL :: outside
    INTEGER :: i, j, k, GSLerr, i1, i2, jout, ktemp
    INTEGER, ALLOCATABLE :: outer(:,:)
    REAL(DP) :: xpsitot, xpl, psis, psitemp, xpsitemp, rtest, dout, psiRatio
    REAL(DP), ALLOCATABLE :: xOldTheta(:), yOldTheta(:), zOldTheta(:), chiValOld(:), &
                             radius(:), xOldPsi(:), yOldPsi(:), zOldPsi(:), psiOld(:), &
                             xtemp(:), ytemp(:), ztemp(:), psiValTemp(:), rtemp(:), dj(:), &
                             xatemp(:), yatemp(:), zatemp(:), phi(:), xPhi(:), yPhi(:), &
                             zPhi(:), xout(:,:), yout(:,:), zout(:,:), rout(:,:), cVal(:), &
                             xmid(:,:), ymid(:,:), zmid(:,:), rmid(:,:), rold(:,:,:)

    ! Variables for tracing
    REAL(DP) :: Pdyn, Dst, ByIMF, BzIMF, G(3), W(6)
    REAL(DP) :: x0, y0, z0, xf, yf, zf, r0, rt, tt, zt, dut
    REAL(DP), DIMENSION(1000) :: xx, yy, zz, bx, by, bz, distance
    INTEGER :: LMAX = 1000, LOUT, scanLeft, scanRight
    INTEGER :: iYear, iMonth, iDay, iHour, iMin, iSec, ifail
    REAL(DP) :: AA, SPS, CPS, PS, AB
    COMMON /GEOPACK1/ AA(10),SPS,CPS,AB(3),PS

    INTEGER :: ID
    REAL(DP) :: Dist, XMGNP, YMGNP, ZMGNP, tempy

    Updated = .false.

    ! No need to do anything when using a dipole
    if ((NameBoundMag.eq.'DIPL').or.(NameBoundMag.eq.'DIPS')) return

    if (NameBoundMag.eq.'SWMF') then
       call computational_domain
       return
    endif

    ! Get correct model inputs and place them in cooresponding variables
    call get_model_inputs(Pdyn,Dst,ByIMF,BzIMF,G,W)

    xpsitemp = xpsiout
    xpsiout = 8.0
    scanLeft  = nZetaMidnight/2
    scanRight = nZetaMidnight + scanLeft
    Find_PsiOut: do k = scanLeft,scanRight
       x0 = 8._dp*cos(zetaVal(k))
       y0 = 8._dp*sin(zetaVal(k))
       z0 = 0._dp
       call trace(x0,y0,z0,1.0_dp,xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
       r0 = xf**2+yf**2+zf**2
       if ((r0.gt.1.1).or.(LOUT.ge.LMAX)) then
          xpsiout = xpsitemp
          exit Find_PsiOut
       endif
       psitemp = 1./(1.-zf**2/r0)
       if (psitemp.lt.xpsiout) then
          xpsiout = psitemp
          ktemp = k
       endif
    enddo Find_PsiOut
    
    !if ((abs(xpsiout-xpsitemp).le.1e-9).or.(xpsitemp.eq.-1._dp)) return

    Updated = .true.

    ! Since T89 tracing is so quick, just retrace everything
    ! this is helpful because T89 jumps are very large
    IF ((NameBoundMag.eq.'T89I').or.(NameBoundMag.eq.'T89D')) THEN
       call Computational_Domain
       return
    ENDIF

    ALLOCATE(xOldTheta(nthe),yOldTheta(nthe),zOldTheta(nthe),chiValOld(nthe),&
             radius(npsi),xOldPsi(npsi),yOldPsi(npsi),zOldPsi(npsi),psiOld(npsi),&
             xtemp(npsi+1),ytemp(npsi+1),ztemp(npsi+1),psiValTemp(npsi+1),rtemp(npsi+1),&
             dj(npsi+1),xatemp(nzeta-1),yatemp(nzeta-1),zatemp(nzeta-1),phi(nzeta+1),&
             xPhi(nzeta+1),yPhi(nzeta+1),zPhi(nzeta+1),xout(nthe,nzeta+1),yout(nthe,nzeta+1),&
             zout(nthe,nzeta+1),rout(nthe,nzeta+1),xmid(nthe,npsi),ymid(nthe,npsi),cVal(nthe),&
             zmid(nthe,npsi),rmid(nthe,npsi),rold(nthe,npsi,nzeta),outer(nthe,nzeta))
    xOldTheta = 0.0; yOldTheta = 0.0; zOldTheta = 0.0; chiValOld = 0.0; radius = 0.0
    xOldPsi = 0.0; yOldPsi = 0.0; zOldPsi = 0.0; psiOld = 0.0; xtemp = 0.0; ytemp = 0.0
    ztemp = 0.0; psiValTemp = 0.0; rtemp = 0.0; dj = 0.0; xatemp = 0.0; yatemp = 0.0
    zatemp = 0.0; phi = 0.0; xPhi = 0.0; yPhi = 0.0; zPhi = 0.0; xout = 0.0; yout = 0.0
    zout = 0.0; rout = 0.0; xmid = 0.0; ymid = 0.0; zmid = 0.0; rmid = 0.0; rold = 0.0
    outer = 0; cVal = 0.0

    ! Trace the outer shell and identify the points that lay outside the new outer boundary
    rtest = 0._dp
    do k=1,nzeta
       r0 = xpsiout
       tt = pi_d-asin(sqrt(1.0/r0))
       rt = r0*sin(tt)**2
       zt = zetaVal(k)
       x0 = rt*cos(zt)*sin(tt)
       y0 = rt*sin(zt)*sin(tt)
       z0 = rt*cos(tt)
       CALL trace(x0,y0,z0,-1.0_dp,xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
       r0 = xf**2+yf**2+zf**2
       if ((r0.gt.1.1).or.(LOUT.ge.LMAX)) then ! If somethings wrong then abort
          Updated = .false.
          DEALLOCATE(xOldTheta,yOldTheta,zOldTheta,chiValOld,radius,xOldPsi,yOldPsi, &
                     zOldPsi,psiOld,xtemp,ytemp,ztemp,psiValTemp,rtemp,dj,xatemp,yatemp, &
                     zatemp,phi,xPhi,yPhi,zPhi,xout,yout,zout,rout,xmid,ymid,zmid,rmid, &
                     rold,outer,cVal)
          return
       endif

       distance(1) = 0._dp
       do i = 2,LOUT
          distance(i) = distance(i-1) + SQRT((xx(i)-xx(i-1))**2 &
                                            +(yy(i)-yy(i-1))**2 &
                                            +(zz(i)-zz(i-1))**2)
       enddo
       cVal = (thetaVal + constTheta * SIN(2.*thetaVal)) * distance(LOUT)/pi_d

       CALL GSL_Interpolation_1D(distance(1:LOUT),xx(1:LOUT),cVal(2:nthe),xout(2:nthe,k),GSLerr)
       CALL GSL_Interpolation_1D(distance(1:LOUT),yy(1:LOUT),cVal(2:nthe),yout(2:nthe,k),GSLerr)
       CALL GSL_Interpolation_1D(distance(1:LOUT),zz(1:LOUT),cVal(2:nthe),zout(2:nthe,k),GSLerr)
       xout(1,k) = x0
       yout(1,k) = y0
       zout(1,k) = z0
       xout(nthe,k) = xf
       yout(nthe,k) = yf
       zout(nthe,k) = zf
       rout(:,k) = xout(:,k)**2+yout(:,k)**2+zout(:,k)**2
       ! Find the points outside of the new magnetic field domain
       outer(:,k) = 1
       do i = 1,nthe
          do j = 1,npsi
             outside = .false.
             rold(i,j,k) = x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2
             if ((rold(i,j,k) > rout(i,k)).and.(j.gt.15)) outside = .true.
             if ((outside).and.(outer(i,k).eq.1)) then
                outer(i,k) = j-2
             endif
          enddo
          if (outer(i,k).le.1) outer(i,k) = npsi-1
       enddo
    enddo

   ! Now move points radially to match the new xpsiout value
    ! Get the new psi values
    kmax = ktemp
    psiin   = -xzero3/xpsiin
    psiout  = -xzero3/xpsiout
    psitot  = psiout-psiin
    xpsitot = xpsiout - xpsiin
    DO j = 1, npsi
       psis = REAL(j-1, DP) / REAL(npsi-1, DP)
       xpl = xpsiin + xpsitot * psis
       psival(j) = -xzero3 / xpl
       f(j) = (xzero3 / xpl**2) * xpsitot !dPsi/dR -- For converting between euler potential
       fp(j) = 0._dp                      !           and the curvilinear coordinate
    END DO

    DO k = 1, nzeta+1
       alphaVal(k) = twopi_d*(REAL(k-2, DP)/REAL(nzeta-1, DP))
       fzet(k)     = 1._dp ! dAlpha/dPhi -- For converting between the euler potential
       fzetp(k)    = 0._dp !                and the curvilinar coordinate
    END DO
   
    ! Map all interior points onto the new psi values
    do k = 1,nzeta
       do i = 2,nthe-1
          jout = outer(i,k)
          xtemp(1:jout) = x(i,1:jout,k)
          ytemp(1:jout) = y(i,1:jout,k)
          ztemp(1:jout) = z(i,1:jout,k)
          psiValTemp(1:jout) = psi(i,1:jout,k)

          xtemp(jout+1) = x(i,jout,k) + (rout(i,k) - rold(i,jout,k)) &
                                       /(rold(i,jout-1,k)-rold(i,jout,k)) &
                                       *(x(i,jout-1,k) - x(i,jout,k))
          if (abs(1-xtemp(jout+1)/x(i,jout,k)).gt.abs(1-rout(i,k)/rold(i,jout,k))) then
             xtemp(jout+1) = x(i,jout,k)*rout(i,k)/rold(i,jout,k)
          endif
          ytemp(jout+1) = y(i,jout,k) + (rout(i,k) - rold(i,jout,k)) &
                                       /(rold(i,jout-1,k)-rold(i,jout,k)) &
                                       *(y(i,jout-1,k) - y(i,jout,k))
          if (abs(1-ytemp(jout+1)/y(i,jout,k)).gt.abs(1-rout(i,k)/rold(i,jout,k))) then
             ytemp(jout+1) = y(i,jout,k)*rout(i,k)/rold(i,jout,k)
          endif
          ztemp(jout+1) = z(i,jout,k) + (rout(i,k) - rold(i,jout,k)) &
                                       /(rold(i,jout-1,k)-rold(i,jout,k)) &
                                       *(z(i,jout-1,k) - z(i,jout,k))
          if (abs(1-ztemp(jout+1)/z(i,jout,k)).gt.abs(1-rout(i,k)/rold(i,jout,k))) then
             ztemp(jout+1) = z(i,jout,k)*rout(i,k)/rold(i,jout,k)
          endif

          dj(1) = 0._dp
          do j = 2,jout
             dj(j) = dj(j-1) + SQRT((x(i,j,k)-x(i,j-1,k))**2 &
                                   +(y(i,j,k)-y(i,j-1,k))**2 &
                                   +(z(i,j,k)-z(i,j-1,k))**2)
          enddo
          dout = dj(jout) + SQRT((xtemp(jout+1)-x(i,jout,k))**2 &
                                +(ytemp(jout+1)-y(i,jout,k))**2 &
                                +(ztemp(jout+1)-z(i,jout,k))**2)
          psiValTemp(jout+1) = psiValTemp(jout) + (dout-dj(jout)) &
                                                 /(dj(jout-1)-dj(jout)) &
                                                 *(psiValTemp(jout-1)-psiValTemp(jout))
          psiRatio = (psiValTemp(jout+1)-psiValTemp(1))/(psiVal(npsi)-psiVal(1))
          do j = 1,jout+1
             psiValTemp(j) = (psiValTemp(j)-psiValTemp(1))/psiRatio + psiVal(1)
          enddo
          call GSL_Interpolation_1D(psiValTemp(1:jout+1),xtemp(1:jout+1),psiVal(1:npsi),x(i,1:npsi,k),GSLerr)
          call GSL_Interpolation_1D(psiValTemp(1:jout+1),ytemp(1:jout+1),psiVal(1:npsi),y(i,1:npsi,k),GSLerr)
          call GSL_Interpolation_1D(psiValTemp(1:jout+1),ztemp(1:jout+1),psiVal(1:npsi),z(i,1:npsi,k),GSLerr)
          x(i,npsi,k) = xtemp(jout+1)
          y(i,npsi,k) = ytemp(jout+1)
          z(i,npsi,k) = ztemp(jout+1)
       enddo
       do j = 1,npsi
          r0 = xpsiin + REAL(j-1,DP)/REAL(npsi-1,DP)*(xpsiout-xpsiin)
          tt = pi_d-asin(sqrt(1.0/r0))
          rt = r0*sin(tt)**2
          zt = atan2(y(2,j,k),x(2,j,k))
          x0 = rt*cos(zt)*sin(tt)
          y0 = rt*sin(zt)*sin(tt)
          z0 = rt*cos(tt)
          x(1,j,k) = x0 
          y(1,j,k) = y0
          z(1,j,k) = z0
          x(nthe,j,k) = x0
          y(nthe,j,k) = y0
          z(nthe,j,k) = -z0
       enddo
    enddo

    DO k = 2, nzeta
       DO j = 1, npsi
          distance(1) = 0._dp
          xOldTheta(:) = x(1:nthe,j,k)
          yOldTheta(:) = y(1:nthe,j,k)
          zOldTheta(:) = z(1:nthe,j,k)
          chiValOld(1) = 0._dp

          DO i = 2, nthe
             distance(i) = distance(i-1) + SQRT((x(i,j,k)-x(i-1,j,k))**2 &
                  & +(y(i,j,k)-y(i-1,j,k))**2 +(z(i,j,k)-z(i-1,j,k))**2)
          END DO

          chiValOld = distance(1:nthe) / distance(nthe) * pi_d

          i1 = 2
          i2 = nthe-1
          CALL GSL_Interpolation_1D(chiValOld,xOldTheta,chiVal(i1:i2),x(i1:i2,j,k),GSLerr)
          CALL GSL_Interpolation_1D(chiValOld,yOldTheta,chiVal(i1:i2),y(i1:i2,j,k),GSLerr)
          CALL GSL_Interpolation_1D(chiValOld,zOldTheta,chiVal(i1:i2),z(i1:i2,j,k),GSLerr)
       END DO
    END DO
    call psiges
    call psifunctions

    !  periodic boundary conditions
    x(:,:,1) = x(:,:,nzeta)
    y(:,:,1) = y(:,:,nzeta)
    z(:,:,1) = z(:,:,nzeta)
    x(:,:,nzeta+1) = x(:,:,2)
    y(:,:,nzeta+1) = y(:,:,2)
    z(:,:,nzeta+1) = z(:,:,2)

    ! For testing the magnetic field update
    !call Write_MAGxyz('UpdateDomain')

    ! Check if the newly calculated magnetic field has any problems
    SORFail = .false.
    call ComputeBandJacob
    do k = 2,nzeta
       do j = 1,npsi
          i = nThetaEquator
          if ((sqrt(x(i,j,k)**2+y(i,j,k)**2+z(i,j,k)**2) < 1.7).or.(SORFail)) then
             write(*,*) 'Issue with calculating new magnetic boundary, using old configuration'
             updated = .false.
          endif
       enddo
    enddo

    DEALLOCATE(xOldTheta,yOldTheta,zOldTheta,chiValOld,radius,xOldPsi,yOldPsi, &
               zOldPsi,psiOld,xtemp,ytemp,ztemp,psiValTemp,rtemp,dj,xatemp,yatemp, &
               zatemp,phi,xPhi,yPhi,zPhi,xout,yout,zout,rout,xmid,ymid,zmid,rmid, &
               rold,outer,cVal)
    return

  end subroutine update_domain
!=============================================================================!
  subroutine trace(x0,y0,z0,RIN,xf,yf,zf,xx,yy,zz,LOUT,LMD,bx,by,bz)
    ! Use Geopack to trace magnetic field lines
    use ModRamParams, ONLY: NameBoundMag

    use nrtype, ONLY: DP

    implicit none

    EXTERNAL :: DIP_08, IGRF_GSW_08, SMGSW_08, T89C, T96_01, T01_01, T04_s
    EXTERNAL :: TS07D_JULY_2017
    integer, intent(in)   :: LMD
    REAL(DP), intent(in)  :: x0, y0, z0, RIN
    integer, intent(out)  :: LOUT
    REAL(DP), intent(inout) :: xf, yf, zf, xx(:), yy(:), zz(:), bx(:), by(:), bz(:)

    integer  :: i, LMAX
    REAL(DP) :: xGSW, yGSW, zGSW, xfGSW, yfGSW, zfGSW, R0, DIR, DSMAX, ER, RLIM
    REAL(DP), ALLOCATABLE :: xxGSW(:), yyGSW(:), zzGSW(:), BxGSW(:), ByGSW(:), BzGSW(:)

    DIR = -1.0
    DSMAX = 0.1
    ER = 0.001
    RLIM = 30.0
    LMAX = LMD

    if (LMAX.lt.0) then
       LMAX = -LMAX
       DIR = -DIR
    endif

    ALLOCATE(xxGSW(LMAX), yyGSW(LMAX), zzGSW(LMAX), BxGSW(LMAX), ByGSW(LMAX), BzGSW(LMAX))
    xxGSW = 0.0; yyGSW = 0.0; zzGSW = 0.0; BxGSW = 0.0; ByGSW = 0.0; BzGSW = 0.0

    CALL SMGSW_08(x0,y0,z0,XGSW,YGSW,ZGSW,1)
    if (RIN.lt.0) then
       R0 = SQRT(XGSW**2+YGSW**2+ZGSW**2)
    else
       R0 = RIN
    endif

    IF (NameBoundMag.eq.'T89D') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T89C,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T89I') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T89C,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T96D') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T96_01,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T96I') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T96_01,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T02D') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T01_01,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T02I') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T01_01,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T04D') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T04_s,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T04I') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     T04_s,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T07D') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     TS07D_JULY_2017,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'T07I') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     TS07D_JULY_2017,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF (NameBoundMag.eq.'IGRF') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     DUMMY,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ELSEIF ((NameBoundMag.eq.'DIPL').or.(NameBoundMag.eq.'DIPS')) then
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     DUMMY,DIP_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ENDIF

    CALL SMGSW_08(xf,yf,zf,xfGSW,yfGSW,zfGSW,-1)
    do i=1,LOUT
       CALL SMGSW_08(xx(i),yy(i),zz(i),xxGSW(i),yyGSW(i),zzGSW(i),-1)
       CALL SMGSW_08(Bx(i),By(i),Bz(i),BxGSW(i),ByGSW(i),BzGSW(i),-1)
    enddo

    DEALLOCATE(xxGSW,yyGSW,zzGSW,BxGSW,ByGSW,BzGSW)
    return

  end subroutine trace

!=============================================================================!
  subroutine DUMMY(IOPT,PARMOD,PSI,X,Y,Z,BXGSW,BYGSW,BZGSW)
    use nrtype, ONLY: DP
    
    implicit none

    integer :: iopt
    real(DP) :: parmod(10), x, y, z, bxgsw, bygsw, bzgsw, psi

    IOPT = IOPT
    PARMOD = PARMOD
    PSI = PSI
    X = X; Y = Y; Z = Z
    BXGSW = 0.0
    BYGSW = 0.0
    BZGSW = 0.0

  end subroutine DUMMY

!=============================================================================!
  subroutine get_model_inputs(Pdyn,Dst,ByIMF,BzIMF,G,W)
    use ModRamTiming,    ONLY: TimeRamNow
    use ModRamVariables, ONLY: Kp
    use ModRamParams,    ONLY: NameBoundMag
    use ModSCBParams,    ONLY: QinDentonPath

    use ModTimeConvert, ONLY: time_int_to_real, n_day_of_year
    USE ModIoUnit, ONLY: UNITTMP_

    use nrtype, ONLY: DP

    implicit none

    real(DP), intent(out) :: Pdyn, Dst, ByIMF, BzIMF, G(3), W(6)

    logical :: lExist
    character(len=4)   :: StringFileFolder
    character(len=8)   :: StringFileDate
    character(len=25)  :: TimeBuffer, StringHeader
    character(len=500) :: QDFile
    integer :: Year, Month, Day, FileIndexStart, FileIndexEnd, nIndex
    integer :: iYear, iMonth, iDay, iHour, iMinute, iSecond
    integer :: i, iError, Hour, Mins, Sec, iFail
    real(DP) :: dsA, dsI
    real(DP), allocatable :: Buffer(:,:), BufferA(:), nSeconds(:)
    REAL(DP) :: AA, SPS, CPS, PS, AB, dut
    COMMON /GEOPACK1/ AA(10),SPS,CPS,AB(3),PS

    ! Usual debug variables.
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='get_model_inputs'
    !==============================================================================================
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    Year   = TimeRamNow%iYear
    Month  = TimeRamNow%iMonth
    Day    = TimeRamNow%iDay
    Hour   = TimeRamNow%iHour
    Mins   = TimeRamNow%iMinute
    Sec    = TimeRamNow%iSecond
    call time_int_to_real(TimeRamNow)
    call RECALC_08(Year,n_day_of_year(Year,Month,Day),Hour,Mins,Sec,-400._dp,0._dp,0._dp)
    ! For now set dipole tile angle to 0
    SPS = 0.0
    CPS = 1.0
    PS = 0.0

    Pdyn  = 0.0 
    Dst   = 0.0 
    ByIMF = 0.0 
    BzIMF = 0.0 
    IF ((NameBoundMag.ne.'T89I').and.(NameBoundMag.ne.'T89D')) THEN
       write(StringFileDate,'(i4.4,i2.2,i2.2)') Year, Month, Day
       write(StringFileFolder,'(i4.4)') Year
       QDFile = trim(QinDentonPath)//'/QinDenton_'//StringFileDate//'_1min.txt'
       INQUIRE(File= trim(QDFile), EXIST=LExist)
       IF (.not.LExist) then
          QDFile = trim(QinDentonPath)//StringFileFolder//'/QinDenton_'//StringFileDate//'_1min.txt'
       ENDIF
       open(unit=UNITTMP_, file=QDFile, status='OLD', iostat=iError)
       if(iError/=0) call CON_stop('get_model_inputs: Error opening file '//trim(QDFile))
       FileIndexStart = 0
       FileIndexEnd = 0
       nIndex = 0
       Read_QDFile_Dates: DO
          read(UNITTMP_,*,IOSTAT=iError) TimeBuffer
          if ((trim(TimeBuffer).ne.'#').and.(FileIndexStart.eq.0)) FileIndexStart = nIndex
          if (iError.lt.0) then
             FileIndexEnd = nIndex
             exit Read_QDFile_Dates
          else
             nIndex = nIndex + 1
             cycle Read_QDFile_Dates
          endif
       ENDDO Read_QDFile_Dates
       nIndex = FileIndexEnd-FileIndexStart-1
       close(UNITTMP_)

       open(unit=UNITTMP_, file=QDFile, status='OLD', iostat=iError)
       do i=1,FileIndexStart
          read(UNITTMP_,*) StringHeader
       enddo

       allocate(nSeconds(nIndex),Buffer(nIndex,36), BufferA(36))
       nSeconds = 0.0; Buffer = 0.0; BufferA = 0.0

       i = 1
       Cycle_QDFile: do
          read(UNITTMP_,*) TimeBuffer, iYear, iMonth, iDay, iHour, iMinute, iSecond, Buffer(i,:)
          if (iSecond.eq.60) then
             iMinute = iMinute + 1
             iSecond = 0
          endif
          if (iMinute.eq.60) then
             iHour = iHour + 1
             iMinute = 0
          endif
          call time_int_to_real((/iYear,iMonth,iDay,iHour,iMinute,iSecond,0/),nSeconds(i))
          if (nSeconds(i).ge.TimeRamNow%Time) then  ! Check that we are on or past the time we want
             dsA = nSeconds(i) - TimeRamNow%Time
             if (abs(dsA).le.1e-9) then                     ! Check if we are exactly on the time or past
                BufferA(:) = Buffer(i,:)
             else
                if (i.eq.1) then                    ! Check if we are on the first time step
                   BufferA(:) = Buffer(1,:)
                elseif (i.eq.nIndex) then
                   BufferA(:) = Buffer(nIndex,:)
                else
                   dsA = TimeRamNow%Time - nSeconds(i-1)
                   dsI = nSeconds(i) - nSeconds(i-1)
                   BufferA(:) = Buffer(i-1,:) + (dsA/dsI)*(Buffer(i,:)-Buffer(i-1,:))
                endif
             endif
             !DenP = BufferA(4)
             Pdyn = BufferA(5)
             Dst = BufferA(19)
             ByIMF = BufferA(1)
             BzIMF = BufferA(2)
             G(:) = BufferA(6:8)
             W(:) = BufferA(26:31)
             if (DoTestMe) then
                write(*,*) 'Testing get_model_inputs'
                write(*,*) 'Buffer = ', BufferA(:)
                write(*,*) 'Values = ', Pdyn, Dst, ByIMF, BzIMF, G, W
             endif
             exit Cycle_QDFile
          endif
          i = i + 1
       enddo Cycle_QDFile
       close(UNITTMP_)
    endif

    IOPT = 0
    PARMOD(1) = Pdyn
    PARMOD(2) = Dst
    PARMOD(3) = ByIMF
    PARMOD(4) = BzIMF
    IF ((NameBoundMag.eq.'T89I').or.(NameBoundMag.eq.'T89D')) THEN
       IOPT = min(floor(Kp+0.5),6)+1
    ELSEIF ((NameBoundMag.eq.'T96I').or.(NameBoundMag.eq.'T96D')) THEN
       ! No extra parameters
    ELSEIF ((NameBoundMag.eq.'T02I').or.(NameBoundMag.eq.'T02D')) THEN
       PARMOD(5) = G(1)
       PARMOD(6) = G(2)
    ELSEIF ((NameBoundMag.eq.'T04I').or.(NameBoundMag.eq.'T04D')) THEN
       PARMOD(5) = W(1)
       PARMOD(6) = W(2)
       PARMOD(7) = W(3)
       PARMOD(8) = W(4)
       PARMOD(9) = W(5)
       PARMOD(10) = W(6)
    ELSEIF ((NameBoundMag.eq.'T07I').or.(NameBoundMag.eq.'T07D')) THEN
       dut = Sec+Mins*60+Hour*3600
       call INIT_TS07D_COEFFS(Year,n_day_of_year(Year,Month,Day),dut,ifail)
       call INIT_TS07D_TLPR
    ENDIF

    if (allocated(nSeconds)) deallocate(nSeconds,Buffer,BufferA)
    return
  end subroutine get_model_inputs
!=============================================================================!
!============================= OUTPUT ROUTINES ===============================!
!=============================================================================!
  SUBROUTINE Write_MAGxyz(FileName)
    use ModRamMain,      ONLY: PathScbOut
    use ModRamTiming,    ONLY: TimeRamNow
    use ModScbGrids,     ONLY: nthe, npsi, nzeta
    use ModScbVariables, ONLY: x, y, z, nThetaEquator

    use ModRamFunctions, ONLY: RamFileName

    use ModIoUnit, ONLY: UNITTMP_

    implicit none

    INTEGER :: i, j, k
    character(len=*), intent(in) :: FileName

    open(UNITTMP_,FILE=trim(PathScbOut)//RamFileName(trim(FileName),'dat',TimeRamNow))
    write(UNITTMP_,*) nthe, npsi, nzeta
    do i = 1, nthe
     do j = 1, npsi
      do k = 1, nzeta
       write(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k)
      enddo
     enddo
    enddo
    close(UNITTMP_)

  END SUBROUTINE Write_MAGxyz

!=============================================================================!
SUBROUTINE Write_ionospheric_potential

  use ModRamTiming, ONLY: TimeRamElapsed

  use ModScbMain,      ONLY: prefixOut
  use ModScbGrids,     ONLY: npsi, nzeta
  use ModScbVariables, ONLY: x, y, PhiIono, dPhiIonodAlpha, dPhiIonodBeta, &
                             alphaVal, psiVal, nThetaEquator, bnormal

  USE nrtype
  USE netcdf

  implicit none

  CHARACTER*500 :: filename

  INTEGER :: alphaid, betaid, timeid, alphavarid, betavarid, &
             xeqid, yeqid, xionoid, yionoid, phiionoid, dphiionodalphaid, &
             dphiionodbetaid, timevarid, ncid

  integer :: START(1), COUNT(1)
  integer :: START1D(2)
  INTEGER :: START2D(3) ! For 2-D arrays (+ time)

  INTEGER, SAVE :: iCALLIP = 0

  REAL :: time(1)

  time(1) = TimeRamElapsed

  START = (/iCALLIP+1/)
  COUNT = (/1/)

  START1D = (/1,iCALLIP+1/)
  START2D = (/1,1,iCALLIP+1/)


  fileName = TRIM(ADJUSTL(prefixOut))//'ionospheric_potential.nc'

  First_time_call : IF(iCALLIP == 0) THEN
     CALL check (nf90_create(filename, nf90_clobber, ncid))

     ! Define dimensions
     CALL check(nf90_def_dim(ncid, 'alpha', npsi, alphaid))
     CALL check(nf90_def_dim(ncid, 'beta', nzeta, betaid))
     CALL check(nf90_def_dim(ncid, 'time', nf90_unlimited, timeid))

     ! Define variables
     CALL check(nf90_def_var(ncid, 'alpha', nf90_float,(/alphaid,timeid/),alphavarid))
     CALL check(nf90_put_att(ncid,alphavarid,'title','Magnetic flux-like Euler potential'))

     CALL check(nf90_def_var(ncid, 'beta', nf90_float,(/betaid,timeid/),betavarid))
     CALL check(nf90_put_att(ncid,betavarid,'title','Azimuthal angle-like Euler potential'))

     CALL check(nf90_def_var(ncid, 'time', nf90_float,timeid,timevarid))
     CALL check(nf90_put_att(ncid,timevarid,'title','Time'))

     CALL check(nf90_def_var(ncid, 'xEq', nf90_float, (/alphaid,betaid,timeid/),xeqid))
     CALL check(nf90_put_att(ncid,xeqid,'title','2D array of xEq locations '))

     CALL check(nf90_def_var(ncid, 'yEq', nf90_float, (/alphaid,betaid,timeid/),yeqid))
     CALL check(nf90_put_att(ncid,yeqid,'title','2D array of yEq locations '))

     CALL check(nf90_def_var(ncid, 'xIono', nf90_float, (/alphaid,betaid,timeid/),xionoid))
     CALL check(nf90_put_att(ncid,xionoid,'title','2D array of xIono locations '))

     CALL check(nf90_def_var(ncid, 'yIono', nf90_float, (/alphaid,betaid,timeid/),yionoid))
     CALL check(nf90_put_att(ncid,yionoid,'title','2D array of yIono locations '))

     CALL check(nf90_def_var(ncid, 'PhiIono', nf90_float, (/alphaid,betaid,timeid/),phiionoid))
     CALL check(nf90_put_att(ncid,phiionoid,'title','2D array of phiIono values'))

     CALL check(nf90_def_var(ncid, 'dPhiIonodAlpha', nf90_float, (/alphaid,betaid,timeid/),dphiionodalphaid))
     CALL check(nf90_put_att(ncid,dphiionodalphaid,'title','2D array of dPhi/dAlpha values'))

     CALL check(nf90_def_var(ncid, 'dPhiIonodBeta', nf90_float, (/alphaid,betaid,timeid/),dphiionodbetaid))
     CALL check(nf90_put_att(ncid,dphiionodbetaid,'title','2D array of dPhi/dBeta values'))

    ! End define mode
     CALL check(nf90_enddef(ncid))

  ELSE ! Open existing NetCDF file
     CALL check(nf90_open(filename, nf90_write, ncid))

     CALL check( nf90_inq_dimid(ncid, 'alpha', alphaid))
     CALL check( nf90_inq_dimid(ncid, 'beta', betaid))
     CALL check( nf90_inq_dimid(ncid, 'time', timeid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'alpha',alphavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'beta', betavarid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'time',timevarid))

     CALL CHECK ( NF90_INQ_VARID (NCID, 'xEq', xeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yEq', yeqid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'xIono', xionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'yIono', yionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'PhiIono', phiionoid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'dPhiIonodAlpha', dphiionodalphaid))
     CALL CHECK ( NF90_INQ_VARID (NCID, 'dPhiIonodBeta', dphiionodbetaid))

  END IF First_time_call

 ! Write mode - write at all times
  CALL check(nf90_put_var(ncid, alphavarid, REAL(psiVal(1:npsi)*bnormal), START1D))
  CALL check(nf90_put_var(ncid, betavarid, REAL(alphaVal(1:nzeta)), START1D))
  CALL check(nf90_put_var(ncid, xeqid, REAL(x(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yeqid, REAL(y(nThetaEquator,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, xionoid, REAL(x(1,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, yionoid, REAL(y(1,1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, phiionoid, REAL(PhiIono(1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, dphiionodalphaid, REAL(dPhiIonodAlpha(1:npsi,1:nzeta)),START2D))
  CALL check(nf90_put_var(ncid, dphiionodbetaid, REAL(dPhiIonodBeta(1:npsi,1:nzeta)),START2D))

  CALL check(nf90_put_var(ncid, timevarid, time, START, COUNT))

  CALL check(nf90_close(ncid))

  iCALLIP = iCALLIP+1

  RETURN

CONTAINS
  SUBROUTINE check(status)
    INTEGER, INTENT ( in) :: status

    IF(status /= nf90_noerr) THEN
       PRINT*, 'STATUS = ', status
       PRINT *, TRIM(nf90_strerror(status))
       STOP 2
    END IF
  END SUBROUTINE check

END SUBROUTINE Write_ionospheric_potential

!==============================================================================
  ! Previously test_Convergence_anisotropic
  SUBROUTINE Write_convergence_anisotropic(iter)
  !!!! Module Variables
  use ModRamTiming,    ONLY: TimeRamNow
  use ModScbMain,      ONLY: prefixOut
  USE ModScbGrids,     ONLY: nthe, npsi, nzeta
  USE ModScbVariables, ONLY: x, y, z, bf, bnormal, Jx, Jy, Jz, Bx, By, Bz, &
                             GradPx, GradPy, GradPz, JCrossB, GradP, &
                             vecd, vec1, vec2, vec3, vec4, vec6, vec7, vec8, &
                             vec9, vecr, vecx, alfa, psi, nZetaMidnight, &
                             pnormal, bnormal, pjconst, nThetaEquator
  !!!! Module Subroutine/Function
  use ModRamFunctions, ONLY: RamFileName
  use ModRamGSL,       ONLY: GSL_Derivs
  use ModScbEquation,  ONLY: metric, metrica, newk, newj
  !!!! Share Modules
  USE ModIoUnit, ONLY: UNITTMP_
  !!!! NR Modules
  use nrtype,    ONLY: DP

  implicit none

  INTEGER :: i, j, k
  CHARACTER(len=200) :: FileName

  REAL(DP), ALLOCATABLE :: xRHS(:,:,:), xLHS(:,:,:), rRHS(:,:,:), rLHS(:,:,:)

  character(len=2), intent(in) :: iter
  !**********************************************************************************************************!
  ALLOCATE(xRHS(nthe,npsi,nzeta), xLHS(nthe,npsi,nzeta), &
           rRHS(nthe,npsi,nzeta), rLHS(nthe,npsi,nzeta))
  xRHS = 0.0; xLHS = 0.0; rRHS = 0.0; rLHS = 0.0

  call metric
  call newj
  DO i = 2,nthe-1
     DO j = 2, npsi-1
        DO k = 2, nzeta
           rLHS(i,j,k) = - vecd(i,j,k)*psi(i,j,k) &
                         + vec1(i,j,k)*psi(i-1,j-1,k) &
                         + vec2(i,j,k)*psi(i,j-1,k) &
                         + vec3(i,j,k)*psi(i+1,j-1,k) &
                         + vec4(i,j,k)*psi(i-1,j,k) &
                         + vec6(i,j,k)*psi(i+1,j,k) &
                         + vec7(i,j,k)*psi(i-1,j+1,k) &
                         + vec8(i,j,k)*psi(i,j+1,k) &
                         + vec9(i,j,k)*psi(i+1,j+1,k)
           rRHS(i,j,k) = vecr(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  rLHS(1,:,:) = 0
  rLHS(nthe,:,:) = 0
  rRHS(1,:,:) = 0
  rRHS(nthe,:,:) = 0
  rLHS(:,:,1) = rLHS(:,:,nzeta)
  rRHS(:,:,1) = rRHS(:,:,nzeta)

  call metrica
  call newk
  DO i = 2, nthe-1
     DO j = 2, npsi-1
        DO k = 2, nzeta
           xLHS(i,j,k) = - vecd(i,j,k)*alfa(i,j,k)  &
                         + vec1(i,j,k)*alfa(i-1,j,k-1) &
                         + vec2(i,j,k)*alfa(i,j,k-1) &
                         + vec3(i,j,k)*alfa(i+1,j,k-1) &
                         + vec4(i,j,k)*alfa(i-1,j,k)  &
                         + vec6(i,j,k)*alfa(i+1,j,k)  &
                         + vec7(i,j,k)*alfa(i-1,j,k+1) &
                         + vec8(i,j,k)*alfa(i,j,k+1) &
                         + vec9(i,j,k)*alfa(i+1,j,k+1)
           xRHS(i,j,k) = vecx(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  xLHS(1,:,:) = 0
  xLHS(nthe,:,:) = 0
  xRHS(1,:,:) = 0
  xRHS(nthe,:,:) = 0
  xLHS(:,:,1) = xLHS(:,:,nzeta)
  xRHS(:,:,1) = xRHS(:,:,nzeta)
 
  ! Force balance quantities
  FileName = trim(prefixOut)//'FB_equatorial_'//iter
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) npsi, nzeta
  i = nThetaEquator
  DO j = 1, npsi
     DO k = 1, nzeta
        WRITE(UNITTMP_, *) x(i,j,k), y(i,j,k), z(i,j,k),    &
                           bf(i,j,k)*bnormal,               &
                           jCrossB(i,j,k), gradP(i,j,k),    &
                           rLHS(i,j,k), rRHS(i,j,k),        &
                           xLHS(i,j,k), xRHS(i,j,k),        &
                           Jx(i,j,k)*pjconst, &
                           Jy(i,j,k)*pjconst, &
                           Jz(i,j,k)*pjconst, &
                           Bx(i,j,k)*bnormal, &
                           By(i,j,k)*bnormal, &
                           Bz(i,j,k)*bnormal, &
                           GradPx(i,j,k)*pnormal/6.4, &
                           GradPy(i,j,k)*pnormal/6.4, &
                           GradPz(i,j,k)*pnormal/6.4
     END DO
  END DO
  CLOSE(UNITTMP_)

  FileName = trim(prefixOut)//'FB_midnight_'//iter
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) nthe, npsi
  DO i = 1, nthe
     DO j = 1, npsi
        k = nZetaMidnight
        WRITE(UNITTMP_, *) x(i,j,k), y(i,j,k), z(i,j,k),    &
                           bf(i,j,k)*bnormal,               &
                           jCrossB(i,j,k), gradP(i,j,k),    &
                           rLHS(i,j,k), rRHS(i,j,k),        &
                           xLHS(i,j,k), xRHS(i,j,k),        &
                           Jx(i,j,k)*pjconst, &
                           Jy(i,j,k)*pjconst, &
                           Jz(i,j,k)*pjconst, &
                           Bx(i,j,k)*bnormal, &
                           By(i,j,k)*bnormal, &
                           Bz(i,j,k)*bnormal, &
                           GradPx(i,j,k)*pnormal/6.4, &
                           GradPy(i,j,k)*pnormal/6.4, &
                           GradPz(i,j,k)*pnormal/6.4
        k = 2
        WRITE(UNITTMP_, *) x(i,j,k), y(i,j,k), z(i,j,k),    &  
                           bf(i,j,k)*bnormal,               &
                           jCrossB(i,j,k), gradP(i,j,k),    &  
                           rLHS(i,j,k), rRHS(i,j,k),        &
                           xLHS(i,j,k), xRHS(i,j,k),        &
                           Jx(i,j,k)*pjconst, &
                           Jy(i,j,k)*pjconst, &
                           Jz(i,j,k)*pjconst, &
                           Bx(i,j,k)*bnormal, &
                           By(i,j,k)*bnormal, &
                           Bz(i,j,k)*bnormal, &
                           GradPx(i,j,k)*pnormal/6.4, &
                           GradPy(i,j,k)*pnormal/6.4, &
                           GradPz(i,j,k)*pnormal/6.4
     END DO
  END DO
  CLOSE(UNITTMP_)

  DEALLOCATE(xRHS, xLHS, rRHS, rLHS)
  RETURN

  END SUBROUTINE Write_convergence_anisotropic

!==================================================================================================
  subroutine write_scb_pressure
    !!!! Module Variables
    USE ModRamTiming,    ONLY: TimeRamNow
    USE ModScbMain,      ONLY: prefixOut
    USE ModScbParams,    ONLY: Isotropy
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: x, y, z, pressure3D, pnormal, dPdPsi, dPdAlpha, &
                               dPPerdPsi, dPPerdAlpha
    
    !!!! Module Subroutine/Function
    use ModRamFunctions, ONLY: RamFileName
  
    !!!! Share Modules
    USE ModIoUnit, ONLY: UNITTMP_

    implicit none

    integer :: i, j, k
    character(len=200) :: FileName

    FileName = trim(prefixOut)//'Pressure3D'
    OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
    write(UNITTMP_, *) nthe, npsi, nzeta
    WRITE(UNITTMP_, *) "X (Re)    Y (Re)    Z (Re)    P (nPa)"
    if (isotropy == 1) then
       DO i = 1,nthe
          DO j = 1,npsi
             DO k = 1,nzeta
                WRITE(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k), pressure3D(i,j,k)*pnormal, &
                                  dPdPsi(i,j,k), dPdAlpha(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    else
       DO i = 1,nthe
          DO j = 1,npsi
             DO k = 1,nzeta
                WRITE(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k), pressure3D(i,j,k)*pnormal, &
                                  dPPerdPsi(i,j,k), dPPerdAlpha(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    endif
    CLOSE(UNITTMP_)

    return

  end subroutine write_scb_pressure
END MODULE ModScbIO
