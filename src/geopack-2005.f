c

c
c          ##########################################################################
c          #                                                                        #
c          #                             GEOPACK-2005                               #
c          #                     (MAIN SET OF FORTRAN CODES)                        #
c          #                                                                        #
c          ##########################################################################
C
c
c  This collection of subroutines is a result of several upgrades of the original package
c  written by N. A. Tsyganenko in 1978-1979. This version is dated May 04, 2005. On that
c  date, the IGRF coefficients were updated according to the recently published table of
c  IGRF-10 coefficients, so that the main field model now extends through 2010 (a linear
c  extrapolation is used for 2005 - 2010, based on the table of secular velocities). For
c  more details, see  
c  http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html (revision of 03/22/2005).
c
c
c  Prefatory notes to the version of April 22, 2003:
c
c  This package represents an in-depth revision of the previous version, with significant
c  changes in the format of calling statements. Users should familiarize themselves with
c  the new formats and rules, and accordingly adjust their source codes, as specified
c  below. Please consult the documentation file
c  
c  http://modelweb.gsfc.nasa.gov/magnetos/data-based/Geopack-2005.doc for detailed
c  descriptions of individual subroutines.
c
c  The following changes were made to the previous release of GEOPACK (of Jan 5, 2001).
c
c (1) Subroutine IGRF, calculating the Earth's main field:

c   (a) Two versions of this subroutine are provided here. In the first one (IGRF_GSM)
c     both input (position) and output (field components) are in the Geocentric Solar-
c     Magnetospheric Cartesian coordinates, while the second one (IGRF_GEO) uses sphe-
c     rical geographical (geocentric) coordinates, as in the older releases.

c   (b) updating of all expansion coefficients is now made separately in the s/r RECALC,
c     which also takes into account the secular change of the coefficients within
c     a given year (at the Earth's surface, the rate of the change can reach 7 nT/month).

c   (c) the optimal length of spherical harmonic expansions is now automatically set
c     inside the code, based on the radial distance, so that the deviation from the
c     full-length approximation does not exceed 0.01 nT. (In the previous versions,
c     the upper limit NM of the order of harmonics had to be specified by users),
c
c  (2) Subroutine DIP, calculating the Earth's field in the dipole approximation:

c   (a) no longer accepts the tilt angle via the list of formal parameters. Instead,
c     the sine SPS and cosine CPS of that angle are now forwarded into DIP via the
c     first common block /GEOPACK1/.  Accordingly, there are two options: (i) to
c     calculate SPS and CPS by calling RECALC before calling DIP, or (ii) to specify
c     them explicitly. In the last case, SPS and CPS should be specified AFTER the
c     invocation of RECALC (otherwise they will be overridden by those returned by
c     RECALC).

c   (b) the Earth's dipole moment is now calculated by RECALC, based on the table of
c     the IGRF coefficients and their secular variation rates, for a given year and
c     the day of the year, and the obtained value of the moment is forwarded into DIP
c     via the second common block /GEOPACK2/. (In the previous versions, only a single
c     fixed value was provided for the geodipole moment, corresponding to the most
c     recent epoch).
c
c  (3) Subroutine RECALC now consolidates in one module all calculations needed to
c     initialize and update the values of coefficients and quantities that vary in
c     time, either due to secular changes of the main geomagnetic field or as a result
c     of Earth's diurnal rotation and orbital motion around Sun. That allowed us to
c     simplify the codes and make them more compiler-independent.
c
c  (4) Subroutine GEOMAG is now identical in its structure to other coordinate trans-
c     formation subroutines. It no longer invokes RECALC from within GEOMAG, but uses
c     precalculated values of the rotation matrix elements, obtained by a separate
c     external invocation of RECALC. This eliminates possible interference of the
c     two subroutines in the old version of the package.
c
c  (5) Subroutine TRACE (and the subsidiary modules STEP and RHAND):
c
c   (a) no longer needs to specify the highest order of spherical harmonics in the
c     main geomagnetic field expansion - it is now calculated automatically inside the
c     IGRF_GSM (or IGRF_GEO) subroutine.
c
c   (b) the internal field model can now be explicitly chosen by specifying the para-
c      meter INNAME (either IGRF_GSM or DIP).
c
c  (6) A new subroutine BCARSP was added, providing a conversion of Cartesian field
c     components into spherical ones (operation, inverse to that performed by the sub-
c     routine  BSPCAR).
c
c  (7) Two new subroutines were added, SHUETAL_MGNP and T96_MGNP, providing the position
c     of the magnetopause, according to the model of Shue et al. [1998] and the one
c     used in the T96 magnetospheric magnetic field model.
c
c
c----------------------------------------------------------------------------------
c
      SUBROUTINE IGRF_GSM (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE GEOCENTRIC SOLAR
C  MAGNETOSPHERIC COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
C  COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html Revised: 22 March, 2005)
c
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE/TIME (IYEAR,IDAY,IHOUR,MIN,ISEC)
C  WAS CHANGED, THE MODEL COEFFICIENTS AND GEO-GSM ROTATION MATRIX ELEMENTS SHOULD BE UPDATED
c  BY CALLING THE SUBROUTINE RECALC
C
C-----INPUT PARAMETERS:
C
C     XGSM,YGSM,ZGSM - CARTESIAN GSM COORDINATES (IN UNITS RE=6371.2 KM)
C
C-----OUTPUT PARAMETERS:
C
C     HXGSM,HYGSM,HZGSM - CARTESIAN GSM COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C
C     LAST MODIFICATION:  MAY 4, 2005.
C     THIS VERSION OF THE CODE ACCEPTS DATES FROM 1965 THROUGH 2010.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      CALL GEOGSM (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,-1)
      RHO2=XGEO**2+YGEO**2
      R=SQRT(RHO2+ZGEO**2)
      C=ZGEO/R
      RHO=SQRT(RHO2)
      S=RHO/R
      IF (S.LT.1.E-5) THEN
        CF=1.
        SF=0.
      ELSE
        CF=XGEO/RHO
        SF=YGEO/RHO
      ENDIF
C
      PP=1./R
      P=PP
C
C  IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.
         Y=1.
170      Q=P
         Z=D
         BI=0.
         P2=0.
         D2=0.
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.E-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BT=BBT
      IF(S.LT.1.E-5) GOTO 210
      BF=BBF/S
      GOTO 211
210   IF(C.LT.0.) BBF=-BBF
      BF=BBF

211   HE=BR*S+BT*C
      HXGEO=HE*CF-BF*SF
      HYGEO=HE*SF+BF*CF
      HZGEO=BR*C-BT*S

      CALL GEOGSM (HXGEO,HYGEO,HZGEO,HXGSM,HYGSM,HZGSM,1)

      RETURN
      END
C
c==========================================================================================
C
c
      SUBROUTINE IGRF_GEO (R,THETA,PHI,BR,BTHETA,BPHI)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE SPHERICAL GEOGRAPHIC
C  (GEOCENTRIC) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
C  COEFFICIENTS  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised: 22 March, 2005)
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE (IYEAR AND IDAY) WAS CHANGED,
C  THE MODEL COEFFICIENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC
C
C-----INPUT PARAMETERS:
C
C   R, THETA, PHI - SPHERICAL GEOGRAPHIC (GEOCENTRIC) COORDINATES:
C   RADIAL DISTANCE R IN UNITS RE=6371.2 KM, COLATITUDE THETA AND LONGITUDE PHI IN RADIANS
C
C-----OUTPUT PARAMETERS:
C
C     BR, BTHETA, BPHI - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C      (POSITIVE BR OUTWARD, BTHETA SOUTHWARD, BPHI EASTWARD)
C
C     LAST MODIFICATION:  MAY 4, 2005.
C     THIS VERSION OF THE  CODE ACCEPTS DATES FROM 1965 THROUGH 2010.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      C=COS(THETA)
      S=SIN(THETA)
      CF=COS(PHI)
      SF=SIN(PHI)
C
      PP=1./R
      P=PP
C
C  IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.
      D=0.
      BBR=0.
      BBT=0.
      BBF=0.

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.
         Y=1.
170      Q=P
         Z=D
         BI=0.
         P2=0.
         D2=0.
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.E-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BTHETA=BBT
      IF(S.LT.1.E-5) GOTO 210
      BPHI=BBF/S
      RETURN
210   IF(C.LT.0.) BBF=-BBF
      BPHI=BBF

      RETURN
      END
C
c==========================================================================================
c
       SUBROUTINE DIP (XGSM,YGSM,ZGSM,BXGSM,BYGSM,BZGSM)
C
C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH, SPECIFIED BY CALLING SUBROUTINE RECALC (SHOULD BE
C  INVOKED BEFORE THE FIRST USE OF THIS ONE AND IN CASE THE DATE/TIME WAS CHANGED).
C
C--INPUT PARAMETERS: XGSM,YGSM,ZGSM - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C--OUTPUT PARAMETERS: BXGSM,BYGSM,BZGSM - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION: MAY 4, 2005
C
C  AUTHOR: N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ AAA(10),SPS,CPS,BBB(23)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIPMOM=SQRT(G(2)**2+G(3)**2+H(3)**2)

      P=XGSM**2
      U=ZGSM**2
      V=3.*ZGSM*XGSM
      T=YGSM**2
      Q=DIPMOM/SQRT(P+T+U)**5
      BXGSM=Q*((T+U-2.*P)*SPS-V*CPS)
      BYGSM=-3.*YGSM*Q*(XGSM*SPS+ZGSM*CPS)
      BZGSM=Q*((P+T-2.*U)*CPS-V*SPS)
      RETURN
      END

C*******************************************************************
c
      SUBROUTINE SUN (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
C
      IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
C
C================================================================================
c
      SUBROUTINE SPHCAR (R,THETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
C    (THETA AND PHI IN RADIANS).
C
C                  J>0            J<0
C-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,THETA,PHI
C
C  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 (WHEN CONVERTING
C        FROM CARTESIAN TO SPHERICAL COORDS, I.E., FOR J<0)
C
C   LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES AND MORE
C                         COMMENTS ADDED)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=SQRT(SQ+Z**2)
      IF (SQ.NE.0.) GOTO 2
      PHI=0.
      IF (Z.LT.0.) GOTO 1
      THETA=0.
      RETURN
  1   THETA=3.141592654
      RETURN
  2   SQ=SQRT(SQ)
      PHI=ATAN2(Y,X)
      THETA=ATAN2(SQ,Z)
      IF (PHI.LT.0.) PHI=PHI+6.28318531
      RETURN
  3   SQ=R*SIN(THETA)
      X=SQ*COS(PHI)
      Y=SQ*SIN(PHI)
      Z=R*COS(THETA)
      RETURN
      END
C
C===========================================================================
c
      SUBROUTINE BSPCAR (THETA,PHI,BR,BTHETA,BPHI,BX,BY,BZ)
C
C   CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
C-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
C              BR,BTHETA,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
C-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
C
C   LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      S=SIN(THETA)
      C=COS(THETA)
      SF=SIN(PHI)
      CF=COS(PHI)
      BE=BR*S+BTHETA*C
      BX=BE*CF-BPHI*SF
      BY=BE*SF+BPHI*CF
      BZ=BR*C-BTHETA*S
      RETURN
      END
c
C==============================================================================
C
      SUBROUTINE BCARSP (X,Y,Z,BX,BY,BZ,BR,BTHETA,BPHI)
C
CALCULATES SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM
C
C-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
C              BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
C-----OUTPUT:  BR,BTHETA,BPHI - SPHERICAL COMPONENTS OF THE FIELD VECTOR
C
C  NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
C        AND HENCE BTHETA=BX, BPHI=BY
C
C   WRITTEN AND ADDED TO THIS PACKAGE:  APRIL 1, 2003,
C   AUTHOR:   N. A. TSYGANENKO
C
      RHO2=X**2+Y**2
      R=SQRT(RHO2+Z**2)
      RHO=SQRT(RHO2)

      IF (RHO.NE.0.) THEN
        CPHI=X/RHO
        SPHI=Y/RHO
       ELSE
        CPHI=1.
        SPHI=0.
      ENDIF

      CT=Z/R
      ST=RHO/R

      BR=(X*BX+Y*BY+Z*BZ)/R
      BTHETA=(BX*CPHI+BY*SPHI)*CT-BZ*ST
      BPHI=BY*CPHI-BX*SPHI

      RETURN
      END
C
c=====================================================================================
C
      SUBROUTINE RECALC (IYEAR,IDAY,IHOUR,MIN,ISEC)
C
C  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
C      (IGRF MODEL)
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C    IGRF_GEO, IGRF_GSM, DIP, GEOMAG, GEOGSM, MAGSM, SMGSM, GSMGSE, GEIGEO.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE AND TIME.
C
C-----INPUT PARAMETERS:
C
C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C
C-----OUTPUT PARAMETERS:   NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN
C
C    AUTHOR:  N.A. TSYGANENKO
C    DATE:    DEC.1, 1991
c
c    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
C     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
C     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
C     CASES CAUSED RUNTIME ERRORS)
c
C    REVISION OF MAY 3, 2005:
C     The table of IGRF coefficients was extended to include those for the epoch 2005
c       the maximal order of spherical harmonics was also increased up to 13
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
c
C    REVISION OF APRIL 3, 2003:
c     The code now includes preparation of the model coefficients for the subroutines
c       IGRF and GEOMAG. This eliminates the need for the SAVE statements, used in the
c        old versions, making the codes easier and more compiler-independent.
C
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * CGST,SGST,BA(6)
C
C  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
C   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
C  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
C    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
C    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.
C
      DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
     + G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
     + H95(105),G00(105),H00(105),G05(105),H05(105),DG05(45),DH05(45)
c
      DATA G65/0.,-30334.,-2119.,-1662.,2997.,1594.,1297.,-2038.,1292.,
     *856.,957.,804.,479.,-390.,252.,-219.,358.,254.,-31.,-157.,-62.,
     *45.,61.,8.,-228.,4.,1.,-111.,75.,-57.,4.,13.,-26.,-6.,13.,1.,13.,
     *5.,-4.,-14.,0.,8.,-1.,11.,4.,8.,10.,2.,-13.,10.,-1.,-1.,5.,1.,-2.,
     *-2.,-3.,2.,-5.,-2.,4.,4.,0.,2.,2.,0.,39*0./
      DATA H65/0.,0.,5776.,0.,-2016.,114.,0.,-404.,240.,-165.,0.,148.,
     *-269.,13.,-269.,0.,19.,128.,-126.,-97.,81.,0.,-11.,100.,68.,-32.,
     *-8.,-7.,0.,-61.,-27.,-2.,6.,26.,-23.,-12.,0.,7.,-12.,9.,-16.,4.,
     *24.,-3.,-17.,0.,-22.,15.,7.,-4.,-5.,10.,10.,-4.,1.,0.,2.,1.,2.,
     *6.,-4.,0.,-2.,3.,0.,-6.,39*0./
c
      DATA G70/0.,-30220.,-2068.,-1781.,3000.,1611.,1287.,-2091.,1278.,
     *838.,952.,800.,461.,-395.,234.,-216.,359.,262.,-42.,-160.,-56.,
     *43.,64.,15.,-212.,2.,3.,-112.,72.,-57.,1.,14.,-22.,-2.,13.,-2.,
     *14.,6.,-2.,-13.,-3.,5.,0.,11.,3.,8.,10.,2.,-12.,10.,-1.,0.,3.,
     *1.,-1.,-3.,-3.,2.,-5.,-1.,6.,4.,1.,0.,3.,-1.,39*0./
      DATA H70/0.,0.,5737.,0.,-2047.,25.,0.,-366.,251.,-196.,0.,167.,
     *-266.,26.,-279.,0.,26.,139.,-139.,-91.,83.,0.,-12.,100.,72.,-37.,
     *-6.,1.,0.,-70.,-27.,-4.,8.,23.,-23.,-11.,0.,7.,-15.,6.,-17.,6.,
     *21.,-6.,-16.,0.,-21.,16.,6.,-4.,-5.,10.,11.,-2.,1.,0.,1.,1.,3.,
     *4.,-4.,0.,-1.,3.,1.,-4.,39*0./
c
      DATA G75/0.,-30100.,-2013.,-1902.,3010.,1632.,1276.,-2144.,1260.,
     *830.,946.,791.,438.,-405.,216.,-218.,356.,264.,-59.,-159.,-49.,
     *45.,66.,28.,-198.,1.,6.,-111.,71.,-56.,1.,16.,-14.,0.,12.,-5.,
     *14.,6.,-1.,-12.,-8.,4.,0.,10.,1.,7.,10.,2.,-12.,10.,-1.,-1.,4.,
     *1.,-2.,-3.,-3.,2.,-5.,-2.,5.,4.,1.,0.,3.,-1.,39*0./
      DATA H75/0.,0.,5675.,0.,-2067.,-68.,0.,-333.,262.,-223.,0.,191.,
     *-265.,39.,-288.,0.,31.,148.,-152.,-83.,88.,0.,-13.,99.,75.,-41.,
     *-4.,11.,0.,-77.,-26.,-5.,10.,22.,-23.,-12.,0.,6.,-16.,4.,-19.,6.,
     *18.,-10.,-17.,0.,-21.,16.,7.,-4.,-5.,10.,11.,-3.,1.,0.,1.,1.,3.,
     *4.,-4.,-1.,-1.,3.,1.,-5.,39*0./
c
      DATA G80/0.,-29992.,-1956.,-1997.,3027.,1663.,1281.,-2180.,1251.,
     *833.,938.,782.,398.,-419.,199.,-218.,357.,261.,-74.,-162.,-48.,
     *48.,66.,42.,-192.,4.,14.,-108.,72.,-59.,2.,21.,-12.,1.,11.,-2.,
     *18.,6.,0.,-11.,-7.,4.,3.,6.,-1.,5.,10.,1.,-12.,9.,-3.,-1.,7.,2.,
     *-5.,-4.,-4.,2.,-5.,-2.,5.,3.,1.,2.,3.,0.,39*0./
      DATA H80/0.,0.,5604.,0.,-2129.,-200.,0.,-336.,271.,-252.,0.,212.,
     *-257.,53.,-297.,0.,46.,150.,-151.,-78.,92.,0.,-15.,93.,71.,-43.,
     *-2.,17.,0.,-82.,-27.,-5.,16.,18.,-23.,-10.,0.,7.,-18.,4.,-22.,9.,
     *16.,-13.,-15.,0.,-21.,16.,9.,-5.,-6.,9.,10.,-6.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6.,39*0./
c
      DATA G85/0.,-29873.,-1905.,-2072.,3044.,1687.,1296.,-2208.,1247.,
     *829.,936.,780.,361.,-424.,170.,-214.,355.,253.,-93.,-164.,-46.,
     *53.,65.,51.,-185.,4.,16.,-102.,74.,-62.,3.,24.,-6.,4.,10.,0.,21.,
     *6.,0.,-11.,-9.,4.,4.,4.,-4.,5.,10.,1.,-12.,9.,-3.,-1.,7.,1.,-5.,
     *-4.,-4.,3.,-5.,-2.,5.,3.,1.,2.,3.,0.,39*0./
      DATA H85/0.,0.,5500.,0.,-2197.,-306.,0.,-310.,284.,-297.,0.,232.,
     *-249.,69.,-297.,0.,47.,150.,-154.,-75.,95.,0.,-16.,88.,69.,-48.,
     *-1.,21.,0.,-83.,-27.,-2.,20.,17.,-23.,-7.,0.,8.,-19.,5.,-23.,11.,
     *14.,-15.,-11.,0.,-21.,15.,9.,-6.,-6.,9.,9.,-7.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6.,39*0./
c
      DATA G90/0., -29775.,  -1848.,  -2131.,   3059.,   1686.,   1314.,
     *     -2239.,   1248.,    802.,    939.,    780.,    325.,   -423.,
     *       141.,   -214.,    353.,    245.,   -109.,   -165.,    -36.,
     *        61.,     65.,     59.,   -178.,      3.,     18.,    -96.,
     *        77.,    -64.,      2.,     26.,     -1.,      5.,      9.,
     *         0.,     23.,      5.,     -1.,    -10.,    -12.,      3.,
     *         4.,      2.,     -6.,      4.,      9.,      1.,    -12.,
     *         9.,     -4.,     -2.,      7.,      1.,     -6.,     -3.,
     *        -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
     *         3.,      3.,      0.,  39*0./

      DATA H90/0.,      0.,   5406.,      0.,  -2279.,   -373.,      0.,
     *      -284.,    293.,   -352.,      0.,    247.,   -240.,     84.,
     *      -299.,      0.,     46.,    154.,   -153.,    -69.,     97.,
     *         0.,    -16.,     82.,     69.,    -52.,      1.,     24.,
     *         0.,    -80.,    -26.,      0.,     21.,     17.,    -23.,
     *        -4.,      0.,     10.,    -19.,      6.,    -22.,     12.,
     *        12.,    -16.,    -10.,      0.,    -20.,     15.,     11.,
     *        -7.,     -7.,      9.,      8.,     -7.,      2.,      0.,
     *         2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
     *         3.,     -1.,     -6.,   39*0./

      DATA G95/0., -29692.,  -1784.,  -2200.,   3070.,   1681.,   1335.,
     *     -2267.,   1249.,    759.,    940.,    780.,    290.,   -418.,
     *       122.,   -214.,    352.,    235.,   -118.,   -166.,    -17.,
     *        68.,     67.,     68.,   -170.,     -1.,     19.,    -93.,
     *        77.,    -72.,      1.,     28.,      5.,      4.,      8.,
     *        -2.,     25.,      6.,     -6.,     -9.,    -14.,      9.,
     *         6.,     -5.,     -7.,      4.,      9.,      3.,    -10.,
     *         8.,     -8.,     -1.,     10.,     -2.,     -8.,     -3.,
     *        -6.,      2.,     -4.,     -1.,      4.,      2.,      2.,
     *         5.,      1.,      0.,    39*0./

      DATA H95/0.,      0.,   5306.,      0.,  -2366.,   -413.,      0.,
     *      -262.,    302.,   -427.,      0.,    262.,   -236.,     97.,
     *      -306.,      0.,     46.,    165.,   -143.,    -55.,    107.,
     *         0.,    -17.,     72.,     67.,    -58.,      1.,     36.,
     *         0.,    -69.,    -25.,      4.,     24.,     17.,    -24.,
     *        -6.,      0.,     11.,    -21.,      8.,    -23.,     15.,
     *        11.,    -16.,    -4.,      0.,    -20.,     15.,     12.,
     *        -6.,     -8.,      8.,      5.,     -8.,      3.,      0.,
     *         1.,      0.,      4.,      5.,     -5.,     -1.,     -2.,
     *         1.,     -2.,     -7.,    39*0./

      DATA G00/0.,-29619.4, -1728.2, -2267.7,  3068.4,  1670.9,  1339.6,
     *     -2288.,  1252.1,   714.5,   932.3,   786.8,    250.,   -403.,
     *      111.3,  -218.8,   351.4,   222.3,  -130.4,  -168.6,   -12.9,
     *       72.3,    68.2,    74.2,  -160.9,    -5.9,    16.9,   -90.4,
     *       79.0,   -74.0,      0.,    33.3,     9.1,     6.9,     7.3,
     *       -1.2,    24.4,     6.6,    -9.2,    -7.9,   -16.6,     9.1,
     *        7.0,    -7.9,     -7.,      5.,     9.4,      3.,   - 8.4,
     *        6.3,    -8.9,    -1.5,     9.3,    -4.3,    -8.2,    -2.6,
     *        -6.,     1.7,    -3.1,    -0.5,     3.7,      1.,      2.,
     *        4.2,     0.3,    -1.1,     2.7,    -1.7,    -1.9,     1.5,
     *       -0.1,     0.1,    -0.7,     0.7,     1.7,     0.1,     1.2,
     *        4.0,    -2.2,    -0.3,     0.2,     0.9,    -0.2,     0.9,
     *       -0.5,     0.3,    -0.3,    -0.4,    -0.1,    -0.2,    -0.4,
     *       -0.2,    -0.9,     0.3,     0.1,    -0.4,     1.3,    -0.4,
     *        0.7,    -0.4,     0.3,    -0.1,     0.4,      0.,     0.1/


      DATA H00/0.,      0.,  5186.1,      0., -2481.6,  -458.0,      0.,
     *     -227.6,   293.4,  -491.1,      0.,   272.6,  -231.9,   119.8,
     *     -303.8,      0.,    43.8,   171.9,  -133.1,   -39.3,   106.3,
     *         0.,   -17.4,    63.7,    65.1,   -61.2,     0.7,    43.8,
     *         0.,   -64.6,   -24.2,     6.2,     24.,    14.8,   -25.4,
     *       -5.8,     0.0,    11.9,   -21.5,     8.5,   -21.5,    15.5,
     *        8.9,   -14.9,    -2.1,     0.0,   -19.7,    13.4,    12.5,
     *       -6.2,    -8.4,     8.4,     3.8,    -8.2,     4.8,     0.0,
     *        1.7,     0.0,     4.0,     4.9,    -5.9,    -1.2,    -2.9,
     *        0.2,    -2.2,    -7.4,     0.0,     0.1,     1.3,    -0.9,
     *       -2.6,     0.9,    -0.7,    -2.8,    -0.9,    -1.2,    -1.9,
     *       -0.9,     0.0,    -0.4,     0.3,     2.5,    -2.6,     0.7,
     *        0.3,     0.0,     0.0,     0.3,    -0.9,    -0.4,     0.8,
     *        0.0,    -0.9,     0.2,     1.8,    -0.4,    -1.0,    -0.1,
     *        0.7,     0.3,     0.6,     0.3,    -0.2,    -0.5,    -0.9/
     *

      DATA G05/0.,-29556.8, -1671.8, -2340.5,   3047.,  1656.9,  1335.7,
     *    -2305.3,  1246.8,   674.4,   919.8,   798.2,   211.5,  -379.5,
     *      100.2,  -227.6,   354.4,   208.8,  -136.6,  -168.3,   -14.1,
     *       72.9,    69.6,    76.6,  -151.1,   -15.0,    14.7,   -86.4,
     *       79.8,   -74.4,    -1.4,    38.6,    12.3,     9.4,     5.5,
     *        2.0,    24.8,     7.7,   -11.4,    -6.8,   -18.0,    10.0,
     *        9.4,   -11.4,    -5.0,     5.6,     9.8,     3.6,    -7.0,
     *        5.0,   -10.8,    -1.3,     8.7,    -6.7,    -9.2,    -2.2,
     *       -6.3,     1.6,    -2.5,    -0.1,     3.0,     0.3,     2.1,
     *        3.9,    -0.1,    -2.2,     2.9,    -1.6,    -1.7,     1.5,
     *       -0.2,     0.2,    -0.7,     0.5,     1.8,     0.1,     1.0,
     *        4.1,    -2.2,    -0.3,     0.3,     0.9,    -0.4,     1.0,
     *       -0.4,     0.5,    -0.3,    -0.4,     0.0,    -0.4,     0.0,
     *       -0.2,    -0.9,     0.3,     0.3,    -0.4,     1.2,    -0.4,
     *        0.7,    -0.3,     0.4,    -0.1,     0.4,    -0.1,    -0.3/

      DATA H05/0.,     0.0,  5080.0,     0.0, -2594.9,  -516.7,     0.0,
     *     -200.4,   269.3,  -524.5,     0.0,   281.4,  -225.8,   145.7,
     *     -304.7,     0.0,    42.7,   179.8,  -123.0,   -19.5,   103.6,
     *        0.0,   -20.2,    54.7,    63.7,   -63.4,     0.0,    50.3,
     *        0.0,   -61.4,   -22.5,     6.9,    25.4,    10.9,   -26.4,
     *       -4.8,     0.0,    11.2,   -21.0,     9.7,   -19.8,    16.1,
     *        7.7,   -12.8,    -0.1,     0.0,   -20.1,    12.9,    12.7,
     *       -6.7,    -8.1,     8.1,     2.9,    -7.9,     5.9,     0.0,
     *        2.4,     0.2,     4.4,     4.7,    -6.5,    -1.0,    -3.4,
     *       -0.9,    -2.3,    -8.0,     0.0,     0.3,     1.4,    -0.7,
     *       -2.4,     0.9,    -0.6,    -2.7,    -1.0,    -1.5,    -2.0,
     *       -1.4,     0.0,    -0.5,     0.3,     2.3,    -2.7,     0.6,
     *        0.4,     0.0,     0.0,     0.3,    -0.8,    -0.4,     1.0,
     *        0.0,    -0.7,     0.3,     1.7,    -0.5,    -1.0,     0.0,
     *        0.7,     0.2,     0.6,     0.4,    -0.2,    -0.5,    -1.0/

      DATA DG05/0.0,   8.8,    10.8,   -15.0,    -6.9,    -1.0,    -0.3,
     *         -3.1,  -0.9,    -6.8,    -2.5,     2.8,    -7.1,     5.9,
     *         -3.2,  -2.6,     0.4,    -3.0,    -1.2,     0.2,    -0.6,
     *         -0.8,   0.2,    -0.2,     2.1,    -2.1,    -0.4,     1.3,
     *         -0.4,   0.0,    -0.2,     1.1,     0.6,     0.4,    -0.5,
     *          0.9,  -0.2,     0.2,    -0.2,     0.2,    -0.2,     0.2,
     *          0.5,  -0.7,     0.5/

      DATA DH05/0.0,   0.0,   -21.3,     0.0,   -23.3,   -14.0,     0.0,
     *          5.4,  -6.5,    -2.0,     0.0,     2.0,     1.8,     5.6,
     *          0.0,   0.0,     0.1,     1.8,     2.0,     4.5,    -1.0,
     *          0.0,  -0.4,    -1.9,    -0.4,    -0.4,    -0.2,     0.9,
     *          0.0,   0.8,     0.4,     0.1,     0.2,    -0.9,    -0.3,
     *          0.3,   0.0,    -0.2,     0.2,     0.2,     0.4,     0.2,
     *         -0.3,   0.5,     0.4/
C
C
      IY=IYEAR
C
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2010, FOR WHICH THE IGRF COEFFICIENTS
c    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
C      NEAREST LIMITING VALUE AND PRINTS A WARNING:
C
      IF(IY.LT.1965) THEN
       IY=1965
       WRITE (*,10) IYEAR,IY
      ENDIF

      IF(IY.GT.2010) THEN
       IY=2010
       WRITE (*,10) IYEAR,IY
      ENDIF

C
C  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
C  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
C  AND THEIR DERIVATIVES:
c
      DO 20 N=1,14
         N2=2*N-1
         N2=N2*(N2-2)
         DO 20 M=1,N
            MN=N*(N-1)/2+M
20    REC(MN)=FLOAT((N-M)*(N+M-2))/FLOAT(N2)
C
      IF (IY.LT.1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
      IF (IY.LT.1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
      IF (IY.LT.1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
      IF (IY.LT.1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
      IF (IY.LT.1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
      IF (IY.LT.1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
      IF (IY.LT.2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
      IF (IY.LT.2005) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005
C
C       EXTRAPOLATE BEYOND 2005:
C

      DT=FLOAT(IY)+FLOAT(IDAY-1)/365.25-2005.
      DO 40 N=1,105
         G(N)=G05(N)
         H(N)=H05(N)
         IF (N.GT.45) GOTO 40
         G(N)=G(N)+DG05(N)*DT
         H(N)=H(N)+DH05(N)*DT
40    CONTINUE
      GOTO 300
C
C       INTERPOLATE BETWEEEN 1965 - 1970:
C
50    F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1965)/5.
      F1=1.-F2
      DO 55 N=1,105
         G(N)=G65(N)*F1+G70(N)*F2
55       H(N)=H65(N)*F1+H70(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1970 - 1975:
C
60    F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1970)/5.
      F1=1.-F2
      DO 65 N=1,105
         G(N)=G70(N)*F1+G75(N)*F2
65       H(N)=H70(N)*F1+H75(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1975 - 1980:
C
70    F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1975)/5.
      F1=1.-F2
      DO 75 N=1,105
         G(N)=G75(N)*F1+G80(N)*F2
75       H(N)=H75(N)*F1+H80(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1980 - 1985:
C
80    F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1980)/5.
      F1=1.-F2
      DO 85 N=1,105
         G(N)=G80(N)*F1+G85(N)*F2
85       H(N)=H80(N)*F1+H85(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1985 - 1990:
C
90    F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1985)/5.
      F1=1.-F2
      DO 95 N=1,105
         G(N)=G85(N)*F1+G90(N)*F2
95       H(N)=H85(N)*F1+H90(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1990 - 1995:
C
100   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1990)/5.
      F1=1.-F2
      DO 105 N=1,105
         G(N)=G90(N)*F1+G95(N)*F2
105      H(N)=H90(N)*F1+H95(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1995 - 2000:
C
110   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-1995)/5.
      F1=1.-F2
      DO 115 N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
         G(N)=G95(N)*F1+G00(N)*F2
115      H(N)=H95(N)*F1+H00(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2000 - 2005:
C
120   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2000)/5.
      F1=1.-F2
      DO 125 N=1,105
         G(N)=G00(N)*F1+G05(N)*F2
125      H(N)=H00(N)*F1+H05(N)*F2
      GOTO 300
C
C   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
C   THEM BY SCHMIDT NORMALIZATION FACTORS:
C
300   S=1.
      DO 130 N=2,14
         MN=N*(N-1)/2+1
         S=S*FLOAT(2*N-3)/FLOAT(N-1)
         G(MN)=G(MN)*S
         H(MN)=H(MN)*S
         P=S
         DO 130 M=2,N
            AA=1.
            IF (M.EQ.2) AA=2.
            P=P*SQRT(AA*FLOAT(N-M+1)/FLOAT(N+M-2))
            MNN=MN+M-1
            G(MNN)=G(MNN)*P
130         H(MNN)=H(MNN)*P

           G10=-G(2)
           G11= G(3)
           H11= H(3)
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0
C
      SQ=G11**2+H11**2
      SQQ=SQRT(SQ)
      SQR=SQRT(G10**2+SQ)
      SL0=-H11/SQQ
      CL0=-G11/SQQ
      ST0=SQQ/SQR
      CT0=G10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
C
      CALL SUN (IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
C   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
C
      S1=COS(SRASN)*COS(SDEC)
      S2=SIN(SRASN)*COS(SDEC)
      S3=SIN(SDEC)
      CGST=COS(GST)
      SGST=SIN(GST)
C
C  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
C   IN THE SYSTEM GEI:
C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM GEI
C   BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT LENGTH:
C
      Y1=DIP2*S3-DIP3*S2
      Y2=DIP3*S1-DIP1*S3
      Y3=DIP1*S2-DIP2*S1
      Y=SQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
C
C   THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
C    HAS THE COMPONENTS:
C
      Z1=S2*Y3-S3*Y2
      Z2=S3*Y1-S1*Y3
      Z3=S1*Y2-S2*Y1
C
C    THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
C     COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
C   THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL HANDBOOKS).
C    HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:
C
      DJ=FLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY)
     * -0.5+FLOAT(IHOUR*3600+MIN*60+ISEC)/86400.
      T=DJ/36525.
      OBLIQ=(23.45229-0.0130125*T)/57.2957795
      DZ1=0.
      DZ2=-SIN(OBLIQ)
      DZ3=COS(OBLIQ)
C
C  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S :
C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
C
C   THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
C   CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE), EM32=(EZGSM,EYGSE)=-EM23,
C     AND EM33=(EZGSM,EZGSE)=EM22
C
      CHI=Y1*DY1+Y2*DY2+Y3*DY3
      SHI=Y1*DZ1+Y2*DZ2+Y3*DZ3
      HI=ASIN(SHI)
C
C    TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)
C
      SPS=DIP1*S1+DIP2*S2+DIP3*S3
      CPS=SQRT(1.-SPS**2)
      PSI=ASIN(SPS)
C
C    THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
C CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED AS FOLLOWS:
C
C IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THE COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
C  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:
C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
C
      XMUT=(ATAN2(SFI,CFI)+3.1415926536)*3.8197186342
C
C  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:
C
C   A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
C   A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
C   A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),
C
C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
C                                                           AND  THEREFORE:
C
      A11=S1*CGST+S2*SGST
      A12=-S1*SGST+S2*CGST
      A13=S3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
C
 10   FORMAT(//1X,
     *'**** RECALC WARNS: YEAR IS OUT OF INTERVAL 1965-2010: IYEAR=',I4,
     * /,6X,'CALCULATIONS WILL BE DONE FOR IYEAR=',I4,/)
      RETURN
      END
c
c====================================================================
C
      SUBROUTINE GEOMAG (XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J)
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
C
C                    J>0                       J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
C-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO
C
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOMAG IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR AND/OR IDAY HAVE BEEN CHANGED
C
C
C   LAST MOFIFICATION:  MARCH 30, 2003 (INVOCATION OF RECALC INSIDE THIS S/R WAS REMOVED)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB(19),BB(8)

      IF(J.GT.0) THEN
       XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
       YMAG=YGEO*CL0-XGEO*SL0
       ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      ELSE
       XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
       YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
       ZGEO=ZMAG*CT0-XMAG*ST0
      ENDIF

      RETURN
      END
c
c=========================================================================================
c
      SUBROUTINE GEIGEO (XGEI,YGEI,ZGEI,XGEO,YGEO,ZGEO,J)
C
C   CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
C   OR VICA VERSA.
C                    J>0                J<0
C----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
C----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEIGEO IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  MARCH 31, 2003

C     AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ A(27),CGST,SGST,B(6)
C
      IF(J.GT.0) THEN
       XGEO=XGEI*CGST+YGEI*SGST
       YGEO=YGEI*CGST-XGEI*SGST
       ZGEO=ZGEI
      ELSE
       XGEI=XGEO*CGST-YGEO*SGST
       YGEI=YGEO*CGST+XGEO*SGST
       ZGEI=ZGEO
      ENDIF

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE MAGSM (XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J)
C
C  CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICA VERSA
C
C                    J>0              J<0
C----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
C----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE MAGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ A(8),SFI,CFI,B(7),AB(10),BA(8)
C
      IF (J.GT.0) THEN
       XSM=XMAG*CFI-YMAG*SFI
       YSM=XMAG*SFI+YMAG*CFI
       ZSM=ZMAG
      ELSE
       XMAG=XSM*CFI+YSM*SFI
       YMAG=YSM*CFI-XSM*SFI
       ZMAG=ZSM
      ENDIF

      RETURN
      END
C
C=======================================================================================
C
       SUBROUTINE GSMGSE (XGSM,YGSM,ZGSM,XGSE,YGSE,ZGSE,J)
C
C CONVERTS GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDS TO SOLAR ECLIPTIC (GSE) ONES
C   OR VICA VERSA.
C                    J>0                J<0
C-----INPUT: J,XGSM,YGSM,ZGSM    J,XGSE,YGSE,ZGSE
C----OUTPUT:   XGSE,YGSE,ZGSE      XGSM,YGSM,ZGSM
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GSMGSE IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ A(12),SHI,CHI,AB(13),BA(8)
C
      IF (J.GT.0) THEN
       XGSE=XGSM
       YGSE=YGSM*CHI-ZGSM*SHI
       ZGSE=YGSM*SHI+ZGSM*CHI
      ELSE
       XGSM=XGSE
       YGSM=YGSE*CHI+ZGSE*SHI
       ZGSM=ZGSE*CHI-YGSE*SHI
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
       SUBROUTINE SMGSM (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,J)
C
C CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR MAGNETOSPHERIC
C   (GSM) COORDINATES OR VICA VERSA.
C                  J>0                 J<0
C-----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
C----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM
C
C  ATTENTION:  SUBROUTINE RECALC  MUST BE INVOKED BEFORE SMGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C     LAST MODIFICATION:  MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/ A(10),SPS,CPS,B(15),AB(8)

      IF (J.GT.0) THEN
       XGSM=XSM*CPS+ZSM*SPS
       YGSM=YSM
       ZGSM=ZSM*CPS-XSM*SPS
      ELSE
       XSM=XGSM*CPS-ZGSM*SPS
       YSM=YGSM
       ZSM=XGSM*SPS+ZGSM*CPS
      ENDIF

      RETURN
      END
C
C==========================================================================================
C
      SUBROUTINE GEOGSM (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,J)
C
C CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR MAGNETOSPHERIC (GSM) COORDINATES
C   OR VICA VERSA.
C
C                   J>0                   J<0
C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSM,YGSM,ZGSM
C---- OUTPUT:    XGSM,YGSM,ZGSM      XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC  MUST BE INVOKED BEFORE GEOGSM IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE BEEN CHANGED
C
C     LAST MODIFICATION: MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      COMMON /GEOPACK1/AA(17),A11,A21,A31,A12,A22,A32,A13,A23,A33,D,B(8)
C
      IF (J.GT.0) THEN
       XGSM=A11*XGEO+A12*YGEO+A13*ZGEO
       YGSM=A21*XGEO+A22*YGEO+A23*ZGEO
       ZGSM=A31*XGEO+A32*YGEO+A33*ZGEO
      ELSE
       XGEO=A11*XGSM+A21*YGSM+A31*ZGSM
       YGEO=A12*XGSM+A22*YGSM+A32*ZGSM
       ZGEO=A13*XGSM+A23*YGSM+A33*ZGSM
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
      SUBROUTINE RHAND (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
C
C  CALCULATES THE COMPONENTS OF THE RIGHT HAND SIDE VECTOR IN THE GEOMAGNETIC FIELD
C    LINE EQUATION  (a subsidiary subroutine for the subroutine STEP)
C
C     LAST MODIFICATION:  MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      DIMENSION PARMOD(10)
C
C     EXNAME AND INNAME ARE NAMES OF SUBROUTINES FOR THE EXTERNAL AND INTERNAL
C     PARTS OF THE TOTAL FIELD
C
      COMMON /GEOPACK1/ A(15),PSI,AA(10),DS3,BB(8)

      CALL EXNAME (IOPT,PARMOD,PSI,X,Y,Z,BXGSM,BYGSM,BZGSM)
      CALL INNAME (X,Y,Z,HXGSM,HYGSM,HZGSM)

      BX=BXGSM+HXGSM
      BY=BYGSM+HYGSM
      BZ=BZGSM+HZGSM
      B=DS3/SQRT(BX**2+BY**2+BZ**2)
      R1=BX*B
      R2=BY*B
      R3=BZ*B
      RETURN
      END
C
C===================================================================================
C
      SUBROUTINE STEP (X,Y,Z,DS,ERRIN,IOPT,PARMOD,EXNAME,INNAME)
C
C   RE-CALCULATES {X,Y,Z}, MAKING A STEP ALONG A FIELD LINE.
C   DS IS THE STEP SIZE, ERRIN IS PERMISSIBLE ERROR VALUE, IOPT SPECIFIES THE EXTERNAL
C   MODEL VERSION, THE ARRAY PARMOD CONTAINS INPUT PARAMETERS FOR THAT MODEL
C   EXNAME IS THE NAME OF THE EXTERNAL FIELD SUBROUTINE
C   INNAME IS THE NAME OF THE INTERNAL FIELD SUBROUTINE (EITHER DIP OR IGRF)
C
C   ALL THE PARAMETERS ARE INPUT ONES; OUTPUT IS THE RENEWED TRIPLET X,Y,Z
C
C     LAST MODIFICATION:  MARCH 31, 2003
C
C     AUTHOR:  N. A. TSYGANENKO
C
      DIMENSION PARMOD(10)
      COMMON /GEOPACK1/ A(26),DS3,B(8)
      EXTERNAL EXNAME,INNAME

  1   DS3=-DS/3.
      CALL RHAND (X,Y,Z,R11,R12,R13,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND (X+R11,Y+R12,Z+R13,R21,R22,R23,IOPT,PARMOD,EXNAME,
     * INNAME)
      CALL RHAND (X+.5*(R11+R21),Y+.5*(R12+R22),Z+.5*
     *(R13+R23),R31,R32,R33,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND (X+.375*(R11+3.*R31),Y+.375*(R12+3.*R32
     *),Z+.375*(R13+3.*R33),R41,R42,R43,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND (X+1.5*(R11-3.*R31+4.*R41),Y+1.5*(R12-
     *3.*R32+4.*R42),Z+1.5*(R13-3.*R33+4.*R43),
     *R51,R52,R53,IOPT,PARMOD,EXNAME,INNAME)
      ERRCUR=ABS(R11-4.5*R31+4.*R41-.5*R51)+ABS(R12-4.5*R32+4.*R42-.5*
     *R52)+ABS(R13-4.5*R33+4.*R43-.5*R53)
      IF (ERRCUR.LT.ERRIN) GOTO 2
      DS=DS*.5
      GOTO 1
  2   X=X+.5*(R11+4.*R41+R51)
      Y=Y+.5*(R12+4.*R42+R52)
      Z=Z+.5*(R13+4.*R43+R53)
      IF(ERRCUR.LT.ERRIN*.04.AND.ABS(DS).LT.1.33) DS=DS*1.5
      RETURN
      END
C
C==============================================================================
C
      SUBROUTINE TRACE (XI,YI,ZI,DIR,RLIM,R0,IOPT,PARMOD,EXNAME,INNAME,
     *XF,YF,ZF,XX,YY,ZZ,L)
C
C  TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
C  SURFACE OR TO A MODEL LIMITING BOUNDARY.
C
C  THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION USED
C  IN THE MAPPING IS CALCULATED AUTOMATICALLY. IF INNAME=IGRF_GSM, THEN AN IGRF MODEL
C  FIELD WILL BE USED, AND IF INNAME=DIP, A PURE DIPOLE FIELD WILL BE USED.

C  IN ANY CASE, BEFORE CALLING TRACE, ONE SHOULD INVOKE RECALC, TO CALCULATE CORRECT
C  VALUES OF THE IGRF COEFFICIENTS AND ALL QUANTITIES NEEDED FOR TRANSFORMATIONS
C  BETWEEN COORDINATE SYSTEMS INVOLVED IN THIS CALCULATIONS.
C
C  ALTERNATIVELY, THE SUBROUTINE RECALC CAN BE INVOKED WITH THE DESIRED VALUES OF
C  IYEAR AND IDAY (TO SPECIFY THE DIPOLE MOMENT), WHILE THE VALUES OF THE DIPOLE
C  TILT ANGLE PSI (IN RADIANS) AND ITS SINE (SPS) AND COSINE (CPS) CAN BE EXPLICITLY
C  SPECIFIED AND FORWARDED TO THE COMMON BLOCK GEOPACK1 (11th, 12th, AND 16th ELEMENTS, RESP.)
C
C------------- INPUT PARAMETERS:
C
C   XI,YI,ZI - GSM COORDS OF INITIAL POINT (IN EARTH RADII, 1 RE = 6371.2 km),
C
C   DIR - SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN WE MOVE ANTIPARALLEL TO THE
C     FIELD VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
C     AND IF DIR=-1.0 THEN THE TRACING GOES IN THE OPPOSITE DIRECTION.
C
C   R0 -  RADIUS OF A SPHERE (IN RE) FOR WHICH THE FIELD LINE ENDPOINT COORDINATES
C     XF,YF,ZF  SHOULD BE CALCULATED.
C
C   RLIM - UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.
C
C   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING AN OPTION OF THE EXTERNAL FIELD
C       MODEL (E.G., INTERVAL OF THE KP-INDEX). ALTERNATIVELY, ONE CAN USE THE ARRAY
C       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT IS JUST A DUMMY PARAMETER.
C
C   PARMOD -  A 10-ELEMENT ARRAY CONTAINING MODEL PARAMETERS, NEEDED FOR A UNIQUE
C      SPECIFICATION OF THE EXTERNAL FIELD. THE CONCRETE MEANING OF THE COMPONENTS
C      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THE EXTERNAL FIELD MODEL.
C
C   EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD
C    (E.G., T96_01).
C   INNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE INTERNAL MAGNETIC FIELD
C    (EITHER DIP OR IGRF_GSM).
C
C-------------- OUTPUT PARAMETERS:
C
C   XF,YF,ZF - GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
C   XX,YY,ZZ - ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. HERE THEIR MAXIMAL LENGTH WAS
C      ASSUMED EQUAL TO 999.
C   L - ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS 999, TRACING
C     TERMINATES, AND A WARNING IS DISPLAYED.
C
C
C     LAST MODIFICATION:  MARCH 31, 2003.
C
C     AUTHOR:  N. A. TSYGANENKO
C
      DIMENSION XX(1000),YY(1000),ZZ(1000), PARMOD(10)
      COMMON /GEOPACK1/ AA(26),DD,BB(8)
      EXTERNAL EXNAME,INNAME
C
      ERR=0.0001
      L=0
      DS=0.5*DIR
      X=XI
      Y=YI
      Z=ZI
      DD=DIR
      AL=0.
c
c  here we call RHAND just to find out the sign of the radial component of the field
c   vector, and to determine the initial direction of the tracing (i.e., either away
c   or towards Earth):
c
      CALL RHAND (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
      AD=0.01
      IF (X*R1+Y*R2+Z*R3.LT.0.) AD=-0.01
C
c     |AD|=0.01 and its sign follows the rule:
c (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
c (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
c     AD is defined in order to initialize the value of RR (radial distance at previous step):

      RR=SQRT(X**2+Y**2+Z**2)+AD
  1   L=L+1
      IF(L.GT.999) GOTO 7
      XX(L)=X
      YY(L)=Y
      ZZ(L)=Z
      RYZ=Y**2+Z**2
      R2=X**2+RYZ
      R=SQRT(R2)

c  check if the line hit the outer tracing boundary; if yes, then terminate
c   the tracing (label 8):

      IF (R.GT.RLIM.OR.RYZ.GT.1600.D0.OR.X.GT.20.D0) GOTO 8
c
c  check whether or not the inner tracing boundary was crossed from outside,
c  if yes, then calculate the footpoint position by interpolation (go to label 6):
c
      IF (R.LT.R0.AND.RR.GT.R) GOTO 6

c  check if (i) we are moving outward, or (ii) we are still sufficiently
c    far from Earth (beyond R=5Re); if yes, proceed further:
c
      IF (R.GE.RR.OR.R.GT.5.) GOTO 5

c  now we moved closer inward (between R=3 and R=5); go to 3 and begin logging
c  previous values of X,Y,Z, to be used in the interpolation (after having
c  crossed the inner tracing boundary):

      IF (R.GE.3.) GOTO 3
c
c  we entered inside the sphere R=3: to avoid too large steps (and hence inaccurate
c  interpolated position of the footpoint), enforce the progressively smaller
c  stepsize values as we approach the inner boundary R=R0:
c
      FC=0.2
      IF(R-R0.LT.0.05) FC=0.05
      AL=FC*(R-R0+0.2)
      DS=DIR*AL
      GOTO 4
  3   DS=DIR
  4   XR=X
      YR=Y
      ZR=Z
  5   RR=R
      CALL STEP (X,Y,Z,DS,ERR,IOPT,PARMOD,EXNAME,INNAME)
      GOTO 1
c
c  find the footpoint position by interpolating between the current and previous
c   field line points:
c
  6   R1=(R0-R)/(RR-R)
      X=X-(X-XR)*R1
      Y=Y-(Y-YR)*R1
      Z=Z-(Z-ZR)*R1
      GOTO 8
  7   WRITE (*,10)
      L=999
  8   XF=X
      YF=Y
      ZF=Z
      RETURN
 10   FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE ARE',
     *' TERMINATED: THE CURRENT NUMBER OF POINTS EXCEEDED 1000 ****'//)
      END
c
C====================================================================================
C
      SUBROUTINE SHUETAL_MGNP(XN_PD,VEL,BZIMF,XGSM,YGSM,ZGSM,
     *  XMGNP,YMGNP,ZMGNP,DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH COORDINATES (XGSM,YGSM,ZGSM) AND SPECIFIED CONDITIONS
C  IN THE INCOMING SOLAR WIND, THIS SUBROUTINE:
C
C (1) DETERMINES IF THE POINT (XGSM,YGSM,ZGSM) LIES INSIDE OR OUTSIDE THE
C      MODEL MAGNETOPAUSE OF SHUE ET AL. (JGR-A, V.103, P. 17691, 1998).
C
C (2) CALCULATES THE GSM POSITION OF A POINT {XMGNP,YMGNP,ZMGNP}, LYING AT THE MODEL
C      MAGNETOPAUSE AND ASYMPTOTICALLY TENDING TO THE NEAREST BOUNDARY POINT WITH
C      RESPECT TO THE OBSERVATION POINT {XGSM,YGSM,ZGSM}, AS IT APPROACHES THE MAGNETO-
C      PAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         BZIMF - IMF BZ IN NANOTESLAS
C
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSM,YGSM,ZGSM - GSM POSITION OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT
C          DIST - DISTANCE (IN RE) BETWEEN THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
C                 AND THE MODEL NAGNETOPAUSE
C          ID -  POSITION FLAG:  ID=+1 (-1) MEANS THAT THE OBSERVATION POINT
C          LIES INSIDE (OUTSIDE) OF THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  OTHER SUBROUTINES USED: T96_MGNP
C
c          AUTHOR:  N.A. TSYGANENKO,
C          DATE:    APRIL 4, 2003.
C
      IF (VEL.LT.0.) THEN
        PD=XN_PD
      ELSE
        PD=1.94E-6*XN_PD*VEL**2  ! PD IS THE SOLAR WIND DYNAMIC PRESSURE (IN nPa)
      ENDIF

c
c  DEFINE THE ANGLE PHI, MEASURED DUSKWARD FROM THE NOON-MIDNIGHT MERIDIAN PLANE;
C  IF THE OBSERVATION POINT LIES ON THE X AXIS, THE ANGLE PHI CANNOT BE UNIQUELY
C  DEFINED, AND WE SET IT AT ZERO:
c
      IF (YGSM.NE.0..OR.ZGSM.NE.0.) THEN
         PHI=ATAN2(YGSM,ZGSM)
      ELSE
         PHI=0.
      ENDIF
C
C  FIRST, FIND OUT IF THE OBSERVATION POINT LIES INSIDE THE SHUE ET AL BDRY
C  AND SET THE VALUE OF THE ID FLAG:
C
      ID=-1
      R0=(10.22+1.29*TANH(0.184*(BZIMF+8.14)))*PD**(-.15151515)
      ALPHA=(0.58-0.007*BZIMF)*(1.+0.024*ALOG(PD))
      R=SQRT(XGSM**2+YGSM**2+ZGSM**2)
      RM=R0*(2./(1.+XGSM/R))**ALPHA
      IF (R.LE.RM) ID=+1
C
C  NOW, FIND THE CORRESPONDING T96 MAGNETOPAUSE POSITION, TO BE USED AS
C  A STARTING APPROXIMATION IN THE SEARCH OF A CORRESPONDING SHUE ET AL.
C  BOUNDARY POINT:
C
      CALL T96_MGNP (PD,-1.,XGSM,YGSM,ZGSM,XMT96,YMT96,ZMT96,DIST,ID96)
C
      RHO2=YMT96**2+ZMT96**2
      R=SQRT(RHO2+XMT96**2)
      ST=SQRT(RHO2)/R
      CT=XMT96/R
C
C  NOW, USE NEWTON'S ITERATIVE METHOD TO FIND THE NEAREST POINT AT THE
C   SHUE ET AL.'S BOUNDARY:
C
      NIT=0

  1   T=ATAN2(ST,CT)
      RM=R0*(2./(1.+CT))**ALPHA

      F=R-RM
      GRADF_R=1.
      GRADF_T=-ALPHA/R*RM*ST/(1.+CT)
      GRADF=SQRT(GRADF_R**2+GRADF_T**2)

      DR=-F/GRADF**2
      DT= DR/R*GRADF_T

      R=R+DR
      T=T+DT
      ST=SIN(T)
      CT=COS(T)

      DS=SQRT(DR**2+(R*DT)**2)

      NIT=NIT+1

      IF (NIT.GT.1000) THEN
         PRINT *,
     *' BOUNDARY POINT COULD NOT BE FOUND; ITERATIONS DO NOT CONVERGE'
      ENDIF

      IF (DS.GT.1.E-4) GOTO 1

      XMGNP=R*COS(T)
      RHO=  R*SIN(T)

      YMGNP=RHO*SIN(PHI)
      ZMGNP=RHO*COS(PHI)

      DIST=SQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE T96_MGNP (XN_PD,VEL,XGSM,YGSM,ZGSM,XMGNP,YMGNP,ZMGNP,
     * DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSM,YGSM,ZGSM), THIS SUBROUTINE DEFINES
C  THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE, HAVING THE
C  SAME VALUE OF THE ELLIPSOIDAL TAU-COORDINATE, AND THE DISTANCE BETWEEN THEM.  THIS IS
C  NOT THE SHORTEST DISTANCE D_MIN TO THE BOUNDARY, BUT DIST ASYMPTOTICALLY TENDS TO D_MIN,
C  AS THE OBSERVATION POINT GETS CLOSER TO THE MAGNETOPAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSM,YGSM,ZGSM - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT, HAVING THE SAME
C          VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
C          DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
C          ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSM,YGSM,ZGSM)
C          LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C  (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)
C
c   AUTHOR:  N.A. TSYGANENKO
C   DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.
C
C
C  DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
C   IF NOT EXPLICITLY SPECIFIED IN THE INPUT:

      IF (VEL.LT.0.) THEN
       PD=XN_PD
      ELSE
       PD=1.94E-6*XN_PD*VEL**2
C
      ENDIF
C
C  RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED EQUAL TO 2 nPa:

      RAT=PD/2.0
      RAT16=RAT**0.14

C (THE POWER INDEX 0.14 IN THE SCALING FACTOR IS THE BEST-FIT VALUE OBTAINED FROM DATA
C    AND USED IN THE T96_01 VERSION)
C
C  VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa:
C
      A0=70.
      S00=1.08
      X00=5.48
C
C   VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED BY THE ACTUAL PRESSURE:
C
      A=A0/RAT16
      S0=S00
      X0=X00/RAT16
      XM=X0-A
C
C  (XM IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE ELLIPSOID AND THE CYLINDER)
C
C     (FOR DETAILS ON THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C      N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C      ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
       IF (YGSM.NE.0..OR.ZGSM.NE.0.) THEN
          PHI=ATAN2(YGSM,ZGSM)
       ELSE
          PHI=0.
       ENDIF
C
       RHO=SQRT(YGSM**2+ZGSM**2)
C
       IF (XGSM.LT.XM) THEN
           XMGNP=XGSM
           RHOMGNP=A*SQRT(S0**2-1)
           YMGNP=RHOMGNP*SIN(PHI)
           ZMGNP=RHOMGNP*COS(PHI)
           DIST=SQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LE.RHO) ID=-1
           RETURN
       ENDIF
C
          XKSI=(XGSM-X0)/A+1.
          XDZT=RHO/A
          SQ1=SQRT((1.+XKSI)**2+XDZT**2)
          SQ2=SQRT((1.-XKSI)**2+XDZT**2)
          SIGMA=0.5*(SQ1+SQ2)
          TAU=0.5*(SQ1-SQ2)
C
C  NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE
C
          XMGNP=X0-A*(1.-S0*TAU)
          ARG=(S0**2-1.)*(1.-TAU**2)
          IF (ARG.LT.0.) ARG=0.
          RHOMGNP=A*SQRT(ARG)
          YMGNP=RHOMGNP*SIN(PHI)
          ZMGNP=RHOMGNP*COS(PHI)
C
C  NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSM,YGSM,ZGSM} AND {XMGNP,YMGNP,ZMGNP}:
C   (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
C    TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):
C
      DIST=SQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LE.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END
C
C===================================================================================
C
c
