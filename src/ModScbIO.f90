MODULE ModScbIO
  
  use ezcdf
  
  implicit none
  save
  
  contains

!=============================================================================!
!============================= INPUT ROUTINES ================================!
!=============================================================================!

  subroutine computational_domain
    !!!! Module Variables  
    USE ModRamVariables, ONLY: Kp
    use ModRamConst,     ONLY: Re
    use ModRamParams,    ONLY: IsComponent, NameBoundMag, boundary
    use ModRamTiming,    ONLY: TimeRamNow
    USE ModScbMain,      ONLY: PathScbIn, blendInitial, tsygcorrect
    USE ModScbParams,    ONLY: decreaseConvAlphaMin, decreaseConvPsiMin, &
                               decreaseConvAlphaMax, decreaseConvPsiMax
    USE ModScbGrids,     ONLY: npsi, nthe, nzeta
    USE ModScbVariables, ONLY: by_imf, bz_imf, dst_global, p_dyn, wTsyg, tilt, constZ, &
                               constTheta, xpsiin, xpsiout, r0Start, byimfglobal, &
                               bzimfglobal, pdynglobal, blendGlobal, blendGlobalInitial, &
                               x, y, z, rhoVal, thetaVal, zetaVal, left, right, chiVal, &
                               kmax, nThetaEquator
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: RamFileName
    use ModScbSpline, ONLY: spline, splint
    use ModScbCouple, ONLY: build_scb_init
    !!!! Share Modules
    use ModTimeConvert, ONLY: n_day_of_year
    USE ModIoUnit, ONLY: UNITTMP_
    !!!! NR Modules
    use nrtype, ONLY: DP, SP, pi_d

    IMPLICIT NONE

    INTEGER :: i, j, k, scans

    REAL(DP) :: ratioFl=1, r0, t0, t1, tt, zt, b, rr, rt
    REAL(DP) :: Pdyn, Dst, ByIMF, BzIMF, G(3), W(6)
    REAL(DP), DIMENSION(1000) :: distance, xx, yy, zz, distance2derivsX, &
                                   distance2derivsY, distance2derivsZ, xxGSW, &
                                   yyGSW, zzGSW, bx, by, bz
    INTEGER :: LMAX = 1000
    INTEGER :: IOPT, LOUT, iYear, iMonth, iDay, iHour, iMin, iSec
    REAL(DP) :: ER, DSMAX, RLIM, PARMOD(10), xf, yf, zf, xf2, yf2, zf2, DIR
    REAL(DP) :: x0, y0, z0, XGSW, YGSW, ZGSW, xfGSW, yfGSW, zfGSW, RIN
    REAL(DP) :: distConsecFluxSqOld, distConsecFluxSq
    REAL(DP) :: AA, SPS, CPS, PS, AB, tVal(nthe), cVal(nthe)
    COMMON /GEOPACK1/ AA(10),SPS,CPS,AB(3),PS

    integer :: time1, clock_rate = 1000, clock_max = 100000
    real(dp) :: starttime,stoptime


    left = 1
    right = npsi
    if ((NameBoundMag.eq.'DIPS').or.(NameBoundMag.eq.'DIPC').or.(NameBoundMag.eq.'DIPL')) then
       ! For generating x, y, and z arrays using analytic dipole and analytic compressed dipole
       ! the variable b controls the compression with 0 being no compression
       constZ = 0.0
       constTheta = 0.0
       xpsiin = 1.75
       xpsiout = 7.00
       b = 0
       if (NameBoundMag.eq.'DIPC') b = 0.4
       DO i = 1, nthe
          tVal(i) = pi_d * REAL(i-1, DP)/REAL(nthe-1, DP)
       END DO
       chival = (tVal + constTheta*sin(2.*tVal))

       do k=1,nzeta
          do j=1,npsi
             r0 = xpsiin + REAL(j-1, DP)/REAL(npsi-1, DP)*(xpsiout-xpsiin)
             rr = (2-b*cos(zetaVal(k)))/(1+b*cos(zetaVal(k)))
             t0 = pi_d-dasin((1.0/r0)**(1./rr))
             t1 = dasin((1.0/r0)**(1./rr))
             do i=1,nthe
                tt = t0 + REAL(i-1,DP)/REAL(nthe-1,DP)*(t1-t0)
                tt = tt + constTheta * SIN(2._dp*tt)
                zt = zetaVal(k)!+constZ*SIN(zetaVal(k))
                rt = r0*dsin(tt)**rr
                x(i,j,k) = (rt)*dcos(zt)*dsin(tt)
                y(i,j,k) = (rt)*dsin(zt)*dsin(tt)
                z(i,j,k) = (rt)*dcos(tt)
             enddo
          enddo
       enddo
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzeta+1) = x(:,:,2)
       y(:,:,nzeta+1) = y(:,:,2)
       z(:,:,nzeta+1) = z(:,:,2)

    else
       ! For generating x, y, and z arrays using field line tracing
       !! Define the inputs needed for the magnetic field models for the tracing
       DIR = -1.0
       DSMAX = 0.1
       ER = 0.001
       RLIM = 20.0
       IOPT = 1
       PARMOD = (/1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
       iYear  = TimeRamNow%iYear
       iMonth = TimeRamNow%iMonth
       iDay   = TimeRamNow%iDay
       iHour  = TimeRamNow%iHour
       iMin   = TimeRamNow%iMinute
       iSec   = TimeRamNow%iSecond
       call RECALC_08(iYear,n_day_of_year(iYear,iMonth,iDay),iHour,iMin,iSec,-400._dp,0._dp,0._dp)

       ! For now set dipole tile angle to 0
       SPS = 0.0
       CPS = 1.0
       PS = 0.0
       !

       ! Get correct model inputs and place them in cooresponding variables
       call get_model_inputs(Pdyn,Dst,ByIMF,BzIMF,G,W)
       IF ((NameBoundMag.eq.'T89I').or.(NameBoundMag.eq.'T89D')) THEN
          IOPT = floor(Kp)
       ELSEIF ((NameBoundMag.eq.'T96I').or.(NameBoundMag.eq.'T96D')) THEN
          PARMOD(1) = Pdyn
          PARMOD(2) = Dst
          PARMOD(3) = ByIMF
          PARMOD(4) = BzIMF
       ELSEIF ((NameBoundMag.eq.'T02I').or.(NameBoundMag.eq.'T02D')) THEN
          PARMOD(1) = Pdyn
          PARMOD(2) = Dst
          PARMOD(3) = ByIMF
          PARMOD(4) = BzIMF
          PARMOD(5) = G(1)
          PARMOD(6) = G(2)
       ELSEIF ((NameBoundMag.eq.'T04I').or.(NameBoundMag.eq.'T04D')) THEN
          PARMOD(1) = Pdyn
          PARMOD(2) = Dst
          PARMOD(3) = ByIMF
          PARMOD(4) = BzIMF
          PARMOD(5) = W(1)
          PARMOD(6) = W(2)
          PARMOD(7) = W(3)
          PARMOD(8) = W(4)
          PARMOD(9) = W(5)
          PARMOD(10) = W(6)
       ELSEIF (NameBoundMag.eq.'IGRF') THEN
          ! Don't need to do anything, just want it to not fail
       ELSE
          CALL CON_STOP('Unrecognized magnetic boundary')
       ENDIF

       ! Start tracing timing
       call system_clock(time1,clock_rate,clock_max)
       starttime=time1/real(clock_rate,dp)
   
       ! Find the correct starting point for the outer edge.
       ! We have to scan the night sector for the most stretched
       ! location since it isn't always at midnight
       scans = 100
       xpsiout = 8.0
       do k = 1,scans
          x0 = 8._dp*dcos(pi_d/2._dp+(k-1)*pi_d/(scans-1))
          y0 = 8._dp*dsin(pi_d/2._dp+(k-1)*pi_d/(scans-1))
          z0 = 0._dp
          call trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD, &
                     xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
          xpsiout = min(xpsiout,1./(1.-zf**2/(xf**2+yf**2+zf**2)))
       enddo

       ! Find the correct starting point for the inner edge.
       ! No need to scan through the night sector since the field
       ! lines won't vary much longitudinally
       !x0 = -1.75_dp
       !y0 = 0._dp
       !z0 = 0._dp
       !call trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,1._dp,IOPT,PARMOD, &
       !           xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
       xpsiin = 1.75  !1./(1.-zf**2/(xf**2+yf**2+zf**2))

       ! Calculate dipole starting points for given xpsiin and xpsiout
       ! in chosen field and perform nzeta*npsi traces to create grid
       constZ = 0.0
       constTheta = 0.0
       DO i = 1, nthe
          tVal(i) = pi_d * REAL(i-1, DP)/REAL(nthe-1, DP)
       END DO
       do k=2,nzeta
          do j=1,npsi
             r0 = xpsiin + REAL(j-1,DP)/REAL(npsi-1,DP)*(xpsiout-xpsiin)
             tt = pi_d-dasin(dsqrt(1.0/r0))
             rt = r0*dsin(tt)**2
             zt = zetaVal(k)!+constZ*sin(zetaVal(k))
             x0 = rt*dcos(zt)*dsin(tt)
             y0 = rt*dsin(zt)*dsin(tt)
             z0 = rt*dcos(tt)
             CALL trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,-1._dp,IOPT,PARMOD, &
                        xf,yf,zf,xx(:),yy(:),zz(:),LOUT,LMAX,bx,by,bz)
             distance(1) = 0._dp
             do i = 2,LOUT
                distance(i) = distance(i-1) + SQRT((xx(i)-xx(i-1))**2 &
                              +(yy(i)-yy(i-1))**2 +(zz(i)-zz(i-1))**2)
             enddo
             cVal = (tVal + constTheta * SIN(2.*tVal)) * distance(LOUT)/pi_d

             call spline(distance(1:LOUT),xx(1:LOUT),1.E31_DP,1.E31_DP,distance2derivsX(1:LOUT))
             call spline(distance(1:LOUT),yy(1:LOUT),1.E31_DP,1.E31_DP,distance2derivsY(1:LOUT))
             call spline(distance(1:LOUT),zz(1:LOUT),1.E31_DP,1.E31_DP,distance2derivsZ(1:LOUT))
             do i = 2,nthe
                x(i,j,k) = splint(distance(1:LOUT),xx(1:LOUT),distance2derivsX(1:LOUT),cval(i))
                y(i,j,k) = splint(distance(1:LOUT),yy(1:LOUT),distance2derivsY(1:LOUT),cval(i))
                z(i,j,k) = splint(distance(1:LOUT),zz(1:LOUT),distance2derivsZ(1:LOUT),cval(i))
             enddo
             x(1,j,k) = x0
             y(1,j,k) = y0
             z(1,j,k) = z0
          enddo
       enddo

       ! Finish tracing timing
       call system_clock(time1,clock_rate,clock_max)
       stoptime=time1/real(clock_rate,dp)
       write(*,*) 'trace: ',stoptime-starttime

       ! Periodic in zeta
       x(:,:,1) = x(:,:,nzeta)
       y(:,:,1) = y(:,:,nzeta)
       z(:,:,1) = z(:,:,nzeta)
       x(:,:,nzeta+1) = x(:,:,2)
       y(:,:,nzeta+1) = y(:,:,2)
       z(:,:,nzeta+1) = z(:,:,2)
       chival = (tVal + constTheta*sin(2.*tVal))
    endif

    distConsecFluxSq = 0._dp
    distConsecFluxSqOld = 0._dp
    DO k = 2, nzeta
       do j = npsi-1, npsi
          distConsecFluxSq = (x(nThetaEquator,j,k) - x(nThetaEquator,j-1,k))**2 &
                            +(y(nThetaEquator,j,k) - y(nThetaEquator,j-1,k))**2
          IF (distConsecFluxSq > distConsecFluxSqOld) THEN
             distConsecFluxSqOld = distConsecFluxSq
             kMax = k
          END IF
       END DO
    END DO

    ! For outputing the outer shell of the magnetic field
    open(UNITTMP_,FILE=RamFileName('MAGxyz','dat',TimeRamNow))
    write(UNITTMP_,*) nthe, npsi, nzeta
    do i = 1,nthe
     do j = 1,npsi
      do k = 1,nzeta
       write(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k)
      enddo
     enddo
    enddo
    close(UNITTMP_)

    return

  end subroutine computational_domain

!=============================================================================!
  subroutine trace(x0,y0,z0,DIR,DSMAX,ER,RLIM,RIN,IOPT,PARMOD,xf,yf,zf, &
                   xx,yy,zz,LOUT,LMAX,bx,by,bz)
    use ModRamParams, ONLY: NameBoundMag

    use nrtype, ONLY: DP
    implicit none

    EXTERNAL :: DIP_08, IGRF_GSW_08, SMGSW_08, T89C, T96_01, T01_01, T04_s
    integer, intent(in)   :: IOPT,LMAX
    REAL(DP), intent(in)  :: x0, y0, z0, ER, DSMAX, RLIM, PARMOD(10), DIR, RIN
    integer, intent(out)  :: LOUT
    REAL(DP), intent(out) :: xf, yf, zf, xx(:), yy(:), zz(:), bx(:), by(:), bz(:)

    integer :: i
    REAL(DP) :: xGSW, yGSW, zGSW, xfGSW, yfGSW, zfGSW, R0
    REAL(DP), DIMENSION(LMAX) :: xxGSW, yyGSW, zzGSW, BxGSW, ByGSW, BzGSW 

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
    ELSEIF (NameBoundMag.eq.'IGRF') THEN
       call TRACE_08(XGSW,YGSW,ZGSW,DIR,DSMAX,ER,RLIM,R0,IOPT,PARMOD, &
                     DUMMY,IGRF_GSW_08,xfGSW,yfGSW,zfGSW,xxGSW(:),yyGSW(:),zzGSW(:), &
                     LOUT,LMAX,BXGSW,BYGSW,BZGSW)
    ENDIF

    CALL SMGSW_08(xf,yf,zf,xfGSW,yfGSW,zfGSW,-1)
    do i=1,LOUT
       CALL SMGSW_08(xx(i),yy(i),zz(i),xxGSW(i),yyGSW(i),zzGSW(i),-1)
       CALL SMGSW_08(Bx(i),By(i),Bz(i),BxGSW(i),ByGSW(i),BzGSW(i),-1)
    enddo

    return

  end subroutine trace

!=============================================================================!
  subroutine DUMMY(IOPT,PARMOD,PSI,X,Y,Z,BXGSW,BYGSW,BZGSW)
    use nrtype, ONLY: DP
    
    implicit none

    integer :: iopt
    real(DP) :: parmod(10), x, y, z, bxgsw, bygsw, bzgsw, psi

    BXGSW = 0.0
    BYGSW = 0.0
    BZGSW = 0.0

  end subroutine DUMMY

!=============================================================================!
  subroutine get_model_inputs(Pdyn,Dst,ByIMF,BzIMF,G,W)
    use ModRamTiming, ONLY: TimeRamNow
    use ModRamParams, ONLY: Optim

    use ModTimeConvert, ONLY: time_int_to_real
    USE ModIoUnit, ONLY: UNITTMP_

    use nrtype, ONLY: DP
    implicit none

    real(DP), intent(out) :: Pdyn, Dst, ByIMF, BzIMF, G(3), W(6)

    character(len=4)   :: StringFileFolder
    character(len=8)   :: StringFileDate
    character(len=25)  :: TimeBuffer, StringHeader
    character(len=500) :: QDFile
    integer :: Year, Month, Day, FileIndexStart, FileIndexEnd, nIndex
    integer :: iYear, iMonth, iDay, iHour, iMinute, iSecond
    integer :: i, iError
    real(DP) :: dsA, dsI
    real(DP), allocatable :: Buffer(:,:), BufferA(:,:), nSeconds(:)

    Year   = TimeRamNow%iYear
    Month  = TimeRamNow%iMonth
    Day    = TimeRamNow%iDay
    call time_int_to_real(TimeRamNow)

    write(StringFileDate,'(i4.4,i2.2,i2.2)') Year, Month, Day
    write(StringFileFolder,'(i4.4)') Year
    if (Optim) QDFile = '/projects/space_data/MagModelInputs/OptimQinDenton_TS04/'
    if (.not.Optim) QDFile = '/projects/space_data/MagModelInputs/QinDenton/'
    QDFile = trim(QDFile)//StringFileFolder//'/QinDenton_'//StringFileDate//'_1min.txt'
write(*,*) 'Reading File: ', QDFile
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

    allocate(nSeconds(nIndex),Buffer(nIndex,36), BufferA(nIndex,36))

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
          if (dsA.eq.0) then                     ! Check if we are exactly on the time or past
             BufferA(i,:) = Buffer(i,:)
          else
             if (i.eq.1) then                    ! Check if we are on the first time step
                BufferA(i,:) = Buffer(i,:)
             else
                dsA = TimeRamNow%Time - nSeconds(i-1)
                dsI = nSeconds(i) - nSeconds(i-1)
                BufferA(i,:) = Buffer(i-1,:) + (dsA/dsI)*(Buffer(i,:)-Buffer(i-1,:))
             endif
          endif
          Pdyn = BufferA(i,5)
          Dst = BufferA(i,19)
          ByIMF = BufferA(i,1)
          BzIMF = BufferA(i,2)
          G(:) = BufferA(i,6:8)
          W(:) = BufferA(i,26:31)
          exit Cycle_QDFile
       elseif (i.eq.nIndex) then
          Pdyn = Buffer(i,5)
          Dst = Buffer(i,19)
          ByIMF = Buffer(i,1)
          BzIMF = Buffer(i,2)
          G(:) = Buffer(i,6:8)
          W(:) = Buffer(i,26:31)
          exit Cycle_QDFile
       endif
       i = i + 1
    enddo Cycle_QDFile

    deallocate(nSeconds,Buffer,BufferA)

  end subroutine get_model_inputs
!=============================================================================!
!============================= OUTPUT ROUTINES ===============================!
!=============================================================================!
SUBROUTINE Write_ionospheric_potential

  use ModRamTiming, ONLY: TimeRamElapsed

  use ModScbMain,      ONLY: prefixOut
  use ModScbGrids,     ONLY: npsi, nzeta
  use ModScbVariables, ONLY: x, y, z, PhiIono, dPhiIonodAlpha, dPhiIonodBeta, &
                             alphaVal, psiVal, nThetaEquator, bnormal

  USE nrtype
  USE netcdf

  IMPLICIT NONE

  CHARACTER*500 :: filename

  INTEGER :: alphaid, betaid, timeid, alphavarid, betavarid, &
       xeqid, yeqid, xionoid, yionoid, phiionoid, dphiionodalphaid, &
       dphiionodbetaid, timevarid, ncid

  integer :: START(1), COUNT(1)
  integer :: START1D(2), COUNT1D(2)
  INTEGER :: START2D(3), COUNT2D(3) ! For 2-D arrays (+ time)

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
  use ModScbParams,    ONLY: isotropy
  USE ModScbGrids,     ONLY: nthe, npsi, nzeta, dt, dr, dpPrime
  USE ModScbVariables, ONLY: thetaVal, rhoVal, zetaVal, x, y, z, &
                             jacobian, normDiff, normGradP, GradZetaSq, &
                             GradThetaGradZeta, GradRhoGradTheta, GradRhoSq, &
                             GradRhoGradZeta, ppar, pper, nThetaEquator, &
                             normJxB, f, fzet, nZetaMidnight, pnormal, &
                             dPPerdRho, dPPerdZeta, dPPerdTheta, bnormal, &
                             pjconst, dPdAlpha, dPdPsi, vecd, vec1, vec2, &
                             vec3, vec4, vec6, vec7, vec8, vec9, vecr, vecx, &
                             alfa, psi, fp, alphaVal, psiVal
  !!!! Module Subroutine/Function
  use ModRamFunctions, ONLY: RamFileName
  use ModScbEquation,  ONLY: metric, metrica, newk, newj
  use ModScbSpline,    ONLY: Spline_coord_derivs
  !!!! Share Modules
  USE ModIoUnit, ONLY: UNITTMP_
  !!!! NR Modules
  use nrtype,    ONLY: DP

  IMPLICIT NONE

  INTEGER :: i, j, k, id, ierr, idealerr
  CHARACTER(len=200) :: FileName

  REAL(DP) :: normDiffRel, volume, bf(nthe,npsi,nzeta+1), bsq(nthe,npsi,nzeta+1), &
              distance(nthe,npsi,nzeta+1)
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: derivXTheta, derivXRho, derivXZeta, &
       derivYTheta, derivYRho, derivYZeta, derivZTheta, derivZRho, derivZZeta, &
       gradRhoX, gradRhoY, gradRhoZ, gradZetaX, gradZetaY, gradZetaZ, gradThetaX, &
       gradThetaY, gradThetaZ, gradThetaSq, derivBsqTheta, derivBsqRho, derivBsqZeta, &
       derivNU1, derivNU2
  ! gradRhoSq, gradRhoGradZeta are global

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jGradRhoPartialTheta, derivjGradRhoPartialTheta, &
       jGradRhoPartialZeta, derivjGradRhoPartialZeta, jGradRho, jGradRhoFactor, jGradZetaPartialRho, &
       derivjGradZetaPartialRho, jGradZetaPartialTheta, derivjGradZetaPartialTheta, jGradZeta, &
       jGradZetaFactor, jGradThetaPartialRho, derivjGradThetaPartialRho, jGradThetaPartialZeta, &
       derivjGradThetaPartialZeta, jGradTheta, jGradThetaFactor, phiToroid, derivPhiRho, derivPhiZeta, &
       derivPhiTheta, derivDiffPTheta

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jCrossBUpRho, jCrossBUpZeta, jCrossBUpTheta, &
       derivjCrossBUpRho, derivjCrossBUpZeta, derivjCrossBUpTheta, jCrossBMinusGradPSq, &
       jCrossBMinusGradPMod, jCrossBSq, jCrossB, gradPSq, gradP

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: jrrInt, jrr, jzzInt, jzz, jrtInt, jrt, jztInt, jzt, &
       rhoCompSq, zetaCompSq, thetaCompSq, curlJCrossBSq, curlJCrossB

  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: xRHS, xLHS, rRHS, rLHS
  REAL(DP), DIMENSION(npsi,nzeta) :: erRHS, erLHS, exRHS, exLHS
  REAL(DP), DIMENSION(nthe,npsi,nzeta) :: Jx, Jy, Jz, Bx, By, Bz, JxBx, JxBy, &
                                          JxBz, GradPx, GradPy, GradPz

  character(len=2), intent(in) :: iter
!  LOGICAL, EXTERNAL :: isnand ! Intrinsic for Portland Group Fortran

  !**********************************************************************************************************!

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, x(1:nthe, 1:npsi, 1:nzeta), &
                           derivXTheta, derivXRho, derivXZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, y(1:nthe, 1:npsi, 1:nzeta), &
                           derivYTheta, derivYRho, derivYZeta)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, z(1:nthe, 1:npsi, 1:nzeta), &
                           derivZTheta, derivZRho, derivZZeta)
  ! Now I have all the point derivatives

  ! Time to build the Jacobian
  jacobian = derivXRho * (derivYZeta * derivZTheta - derivYTheta * derivZZeta) &
           + derivXZeta * (derivYTheta * derivZRho - derivYRho * derivZTheta) &
           + derivXTheta * (derivYRho * derivZZeta - derivYZeta * derivZRho)

  gradRhoX = (derivYZeta * derivZTheta - derivYTheta * derivZZeta) / jacobian
  gradRhoY = (derivZZeta * derivXTheta - derivZTheta * derivXZeta) / jacobian
  gradRhoZ = (derivXZeta * derivYTheta - derivXTheta * derivYZeta) / jacobian

  gradZetaX = (derivYTheta * derivZRho - derivYRho * derivZTheta) / jacobian
  gradZetaY = (derivZTheta * derivXRho - derivZRho * derivXTheta) / jacobian
  gradZetaZ = (derivXTheta * derivYRho - derivXRho * derivYTheta) / jacobian

  gradThetaX = (derivYRho * derivZZeta - derivYZeta * derivZRho) / jacobian
  gradThetaY = (derivZRho * derivXZeta - derivZZeta * derivXRho) / jacobian
  gradThetaZ = (derivXRho * derivYZeta - derivXZeta * derivYRho) / jacobian

  gradRhoSq = gradRhoX**2 + gradRhoY**2 + gradRhoZ**2
  gradRhoGradZeta = gradRhoX * gradZetaX + gradRhoY * gradZetaY + gradRhoZ * gradZetaZ
  gradRhoGradTheta = gradRhoX * gradThetaX + gradRhoY * gradThetaY + gradRhoZ * gradThetaZ

  gradThetaSq = gradThetaX**2 + gradThetaY**2 + gradThetaZ**2
  gradThetaGradZeta = gradThetaX * gradZetaX + gradThetaY * gradZetaY + gradThetaZ * gradZetaZ

  gradZetaSq = gradZetaX**2 + gradZetaY**2 + gradZetaZ**2

  ! B field squared is obtained now, then the B field bf; in the following, i, j, k can be 
  ! taken over the full domain because we have all the required quantities everywhere;
  ! thus, extrapolation for Bfield is not necessary anymore

  DO  k = 1,nzeta
     DO  j = 1,npsi
        DO  i = 1,nthe
           bsq(i,j,k) = (gradRhoSq(i,j,k) * gradZetaSq(i,j,k) - gradRhoGradZeta(i,j,k) **2) &
                      * (f(j) * fzet(k)) **2
           bf(i,j,k) = SQRT(bsq(i,j,k))
        END DO
     END DO
  END DO

  ! j dot gradRho
  DO j = 1, npsi
     DO k = 1, nzeta
        jGradRhoPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * &
             gradRhoGradZeta(:,j,k))
        jGradRhoPartialZeta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoSq(:,j,k) * gradZetaSq(:,j,k) - gradRhoGradZeta(:,j,k) **2)
     END DO
  END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialTheta, &
       & derivjGradRhoPartialTheta, derivNU1, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradRhoPartialZeta, &
       & derivNU1, derivNU2, derivjGradRhoPartialZeta)

  jGradRho = (derivjGradRhoPartialTheta + derivjGradRhoPartialZeta)/jacobian

  ! j dot gradZeta
  DO j = 1, npsi
     DO k = 1, nzeta
        jGradZetaPartialRho(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k)**2 - gradRhoSq(:,j,k) * gradZetaSq(:,j,k))
        jGradZetaPartialTheta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhoGradZeta(:,j,k) * gradThetaGradZeta(:,j,k) - gradRhoGradTheta(:,j,k) * gradZetaSq(:,j,k))
     END DO
  END DO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialRho, &
       & derivNU1, derivjGradZetaPartialRho, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradZetaPartialTheta, &
       & derivjGradZetaPartialTheta, derivNU1, derivNU2)

  jGradZeta = (derivjGradZetaPartialRho + derivjGradZetaPartialTheta)/jacobian

  ! j dot gradTheta
  DO j = 1,npsi
     DO k = 1,nzeta
        jGradThetaPartialRho(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhogradTheta(:,j,k) * gradRhogradZeta(:,j,k) - gradThetagradZeta(:,j,k) * gradRhoSq(:,j,k))
        jGradThetaPartialZeta(:,j,k) = jacobian(:,j,k) * f(j) * fzet(k) * &
             (gradRhogradTheta(:,j,k) * gradZetaSq(:,j,k) - gradRhogradZeta(:,j,k) * gradThetagradZeta(:,j,k))
     ENDDO
  ENDDO

  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradThetaPartialRho, &
       & derivNU1, derivjGradThetaPartialRho, derivNU2)
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jGradThetaPartialZeta, &
       & derivNU1, derivNU2, derivjGradThetaPartialZeta)

  jGradTheta = (derivjGradThetaPartialRho + derivjGradThetaPartialZeta)/jacobian

  ! The one below is for calculating \partial Pper/\partial \theta and
  ! \partial(J(Pper-Ppar))/\partial \theta for the anisotropic case  
  CALL Spline_coord_derivs(thetaVal, rhoVal, zetaVal, jacobian(1:nthe,1:npsi,1:nzeta) * &
       (pper(1:nthe,1:npsi,1:nzeta)-ppar(1:nthe,1:npsi,1:nzeta)), derivDiffPTheta, derivNU1, derivNU2)

  Jx = (JGradRho*derivXRho + JGradZeta*derivXZeta + JGradTheta*derivXTheta)
  Jy = (JGradRho*derivYRho + JGradZeta*derivYZeta + JGradTheta*derivYTheta)
  Jz = (JGradRho*derivZRho + JGradZeta*derivZZeta + JGradTheta*derivZTheta)
  DO k = 1, nzeta
     DO j = 1, npsi
        Bx(:,j,k) = (f(j)*fzet(k)*derivXTheta(:,j,k)/jacobian(:,j,k))
        By(:,j,k) = (f(j)*fzet(k)*derivYTheta(:,j,k)/jacobian(:,j,k))
        Bz(:,j,k) = (f(j)*fzet(k)*derivZTheta(:,j,k)/jacobian(:,j,k))
        if (isotropy.eq.0) then
           GradPx(:,j,k) = (dPperdRho(:,j,k)*GradRhoSq(:,j,k) &
                          + dPperdZeta(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradRhoGradTheta(:,j,k))*derivXRho(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradZetaSq(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaGradZeta(:,j,k))*derivXZeta(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradThetaGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaSq(:,j,k))*derivXTheta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradRhoGradTheta(:,j,k)*derivXRho(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaGradZeta(:,j,k)*derivXZeta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaSq(:,j,k)*derivXTheta(:,j,k)
           GradPy(:,j,k) = (dPperdRho(:,j,k)*GradRhoSq(:,j,k) &
                          + dPperdZeta(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradRhoGradTheta(:,j,k))*derivYRho(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradZetaSq(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaGradZeta(:,j,k))*derivYZeta(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradThetaGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaSq(:,j,k))*derivYTheta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradRhoGradTheta(:,j,k)*derivYRho(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaGradZeta(:,j,k)*derivYZeta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaSq(:,j,k)*derivYTheta(:,j,k)
           GradPz(:,j,k) = (dPperdRho(:,j,k)*GradRhoSq(:,j,k) &
                          + dPperdZeta(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradRhoGradTheta(:,j,k))*derivZRho(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradZetaSq(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaGradZeta(:,j,k))*derivZZeta(:,j,k) &
                         + (dPperdRho(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + dPperdZeta(:,j,k)*GradThetaGradZeta(:,j,k) &
                          + dPperdTheta(:,j,k)*GradThetaSq(:,j,k))*derivZTheta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradRhoGradTheta(:,j,k)*derivZRho(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaGradZeta(:,j,k)*derivZZeta(:,j,k) &
                         + derivDiffPTheta(:,j,k)*GradThetaSq(:,j,k)*derivZTheta(:,j,k)
        else
           GradPx(:,j,k) = (f(j)*dPdPsi(:,j,k)*GradRhoSq(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradRhoGradZeta(:,j,k))*derivXRho(:,j,k) &
                         + (f(j)*dPdPsi(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradZetaSq(:,j,k))*derivXZeta(:,j,k) &
                         + (f(j)*dPdPsi(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradThetaGradZeta(:,j,k))*derivXTheta(:,j,k)
           GradPy(:,j,k) = (f(j)*dPdPsi(:,j,k)*GradRhoSq(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradRhoGradZeta(:,j,k))*derivYRho(:,j,k) &
                         + (f(j)*dPdPsi(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradZetaSq(:,j,k))*derivYZeta(:,j,k) &
                         + (f(j)*dPdPsi(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradThetaGradZeta(:,j,k))*derivYTheta(:,j,k)
           GradPz(:,j,k) = (f(j)*dPdPsi(:,j,k)*GradRhoSq(:,j,k) &
                           + fzet(k)*dPdAlpha(:,j,k)*GradRhoGradZeta(:,j,k))*derivZRho(:,j,k) &
                          + (f(j)*dPdPsi(:,j,k)*GradRhoGradZeta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradZetaSq(:,j,k))*derivZZeta(:,j,k) &
                         + (f(j)*dPdPsi(:,j,k)*GradRhoGradTheta(:,j,k) &
                          + fzet(k)*dPdAlpha(:,j,k)*GradThetaGradZeta(:,j,k))*derivZTheta(:,j,k)
        endif
     ENDDO
  ENDDO
  JxBx = Jy*Bz-Jz*By
  JxBy = Jz*Bx-Jx*Bz
  JxBz = Jx*By-Jy*Bx

  jCrossB = sqrt(JxBx**2+JxBy**2+JxBz**2)
  gradP = sqrt(GradPx**2+GradPy**2+GradPz**2)

  GradPx = GradPx * pnormal*2
  GradPy = GradPy * pnormal*2
  GradPz = GradPz * pnormal*2
  Bx = Bx * bnormal
  By = By * bnormal
  Bz = Bz * bnormal
  Jx = Jx * pjconst*6.4
  Jy = Jy * pjconst*6.4
  Jz = Jz * pjconst*6.4
  jCrossB = jCrossB*bnormal*pjconst*6.4
  GradP   = GradP*pnormal
  jCrossBMinusGradPMod = jCrossB-GradP

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
  FileName = trim(prefixOut)//'Force_balance_equatorial_'//iter
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) npsi-2, nzeta
  i = nThetaEquator
  DO j = 2, npsi-1
     DO k = 1, nzeta
        WRITE(UNITTMP_, *) x(i, j, k), y(i, j, k), bf(i, j, k)*bnormal, &
                           jCrossB(i,j,k), jCrossBMinusGradPMod(i,j,k), &
                           gradP(i,j,k), rLHS(i,j,k), rRHS(i,j,k), xLHS(i,j,k), &
                           xRHS(i,j,k)
     END DO
  END DO
  CLOSE(UNITTMP_)

  FileName = trim(prefixOut)//'Force_balance_midnight_'//iter
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
  WRITE(UNITTMP_, *) npsi-2, nthe
  k = nZetaMidnight
  DO j = 2, npsi-1
     DO i = 1, nthe
        k = nZetaMidnight
        WRITE(UNITTMP_, *) x(i,j,k), z(i,j,k), bf(i,j,k)*bnormal, &
                           jCrossB(i,j,k), jCrossBMinusGradPMod(i,j,k), &
                           gradP(i,j,k), rLHS(i,j,k), rRHS(i,j,k), xLHS(i,j,k), &
                           xRHS(i,j,k)
        k = 2
        WRITE(UNITTMP_, *) x(i,j,k), z(i,j,k), bf(i,j,k)*bnormal, &
                           jCrossB(i,j,k), jCrossBMinusGradPMod(i,j,k), &
                           gradP(i,j,k), rLHS(i,j,k), rRHS(i,j,k), xLHS(i,j,k), &
                           xRHS(i,j,k)
     END DO
  END DO
  CLOSE(UNITTMP_)

  i = nThetaEquator
  k = 8
  FileName = trim(prefixOut)//'Force_balance_line_'//iter
  OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace') 
  WRITE(UNITTMP_, *) i,k
  DO j = 2, npsi-1
     WRITE(UNITTMP_, *) x(i,j,k), y(i,j,k), jCrossB(i,j,k), gradP(i,j,k), rLHS(i,j,k), &
                        rRHS(i,j,k), xLHS(i,j,k), xRHS(i,j,k), Jx(i,j,k), Jy(i,j,k), &
                        Jz(i,j,k), Bx(i,j,k), By(i,j,k), Bz(i,j,k), JxBx(i,j,k), &
                        JxBy(i,j,k), JxBz(i,j,k), GradPx(i,j,k), GradPy(i,j,k), &
                        GradPz(i,j,k), jacobian(i,j,k), z(i,j,k), derivXTheta(i,j,k), &
                        derivXRho(i,j,k), derivXZeta(i,j,k), derivYTheta(i,j,k), &
                        derivYRho(i,j,k), derivYZeta(i,j,k), derivZTheta(i,j,k), &
                        derivZRho(i,j,k), derivZZeta(i,j,k), JGradTheta(i,j,k), &
                        JGradRho(i,j,k), JGradZeta(i,j,k)

  END DO
  CLOSE(UNITTMP_)

  normDiff = 0.0_dp
  normDiffRel = 0.0_dp
  normJxB  = 0.0_dp
  normGradP = 0.0_dp
  volume   = 0.0_dp
  DO i = 2, nthe-1
     DO j = 2, npsi
        DO k = 2, nzeta
           !IF (2.*pper(i,j,k) > 1.E-1_dp*bsq(i,j,k)) THEN
              ! Only in regions with beta > 1E-2 
              ! (in regions of low plasma beta, the pressure does not change the magnetic field)
              normDiff = normDiff + jacobian(i,j,k) * dr * dpPrime * dt * jCrossBMinusGradPMod(i,j,k)
              normDiffRel = normDiffRel + jacobian(i,j,k) * dr * dpPrime * dt * jCrossBMinusGradPMod(i,j,k) / jCrossB(i,j,k)
              normJxB = normJxB + jacobian(i,j,k) * dr * dpPrime * dt * jCrossB(i,j,k)
              normGradP = normGradP + jacobian(i,j,k) * dr * dpPrime * dt * gradP(i,j,k)
              volume = volume + jacobian(i,j,k) * dr * dpPrime * dt
           !END IF
        END DO
     END DO
  END DO

  ! Normalize to total computational volume
  normDiff = normDiff/volume
  normJxB = normJxB/volume
  normGradP = normGradP/volume

  !  Norms of |jxB-grad P|,      |jxB|,      |gradP| 
  WRITE(*, *) normDiff, normJxB, normGradP
  if ((normDiff.ne.normDiff).or.(normJxB.ne.normJxB).or.(normGradP.ne.normGradP)) then
   call con_stop('NaN encountered in ModScbIO.write_convergence subroutine')
  endif
  RETURN

  END SUBROUTINE Write_convergence_anisotropic

!==================================================================================================
  subroutine write_scb_pressure
    !!!! Module Variables
    USE ModRamTiming,    ONLY: TimeRamNow
    USE ModScbMain,      ONLY: prefixOut
    USE ModScbGrids,     ONLY: nthe, npsi, nzeta
    USE ModScbVariables, ONLY: x, y, z, pressure3D, pnormal
    
    !!!! Module Subroutine/Function
    use ModRamFunctions, ONLY: RamFileName
  
    !!!! Share Modules
    USE ModIoUnit, ONLY: UNITTMP_

    implicit none

    integer :: i, j, k
    character(len=200) :: FileName

    FileName = trim(prefixOut)//'Pressure3D'
    OPEN(UNITTMP_, file = RamFileName(FileName,'dat',TimeRamNow), status='replace')
    WRITE(UNITTMP_, *) "X (Re)    Y (Re)    Z (Re)    P (nPa)"
    DO i = 1,nthe
       DO j = 1,npsi
          DO k = 1,nzeta
             WRITE(UNITTMP_,*) x(i,j,k), y(i,j,k), z(i,j,k), pressure3D(i,j,k)*pnormal
          ENDDO
       ENDDO
    ENDDO
    CLOSE(UNITTMP_)

    return

  end subroutine write_scb_pressure
END MODULE ModScbIO
