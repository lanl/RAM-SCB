! Subroutine GEOMAG
!
! Excerpted from the IRI package, Stan Solomon, 6/1988.  Original name was
! GMAG, which was changed to avoid conflict when used at the same time
! as IRI.  Labeled common /CONST/ which contained only the degree/radian
! conversion factor was changed to a data statement.  Longitude range
! from -180 to 180.
!
! Minor refactor to f90, SCS, 2016.
!
! CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT)
! FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
! AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
! LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.
!
    SUBROUTINE GEO2MAG(ART,LONG,LATI,MLONG,MLAT)
      INTEGER ART
      REAL MLONG,MLAT,LONG,LATI
      DATA FAKTOR/.0174532952/
      ZPI=FAKTOR*360.
      CBG=11.4*FAKTOR
      CI=COS(CBG)
      SI=SIN(CBG)
      IF(ART.NE.0) THEN
        CBM=COS(MLAT*FAKTOR)
        SBM=SIN(MLAT*FAKTOR)
        CLM=COS(MLONG*FAKTOR)
        SLM=SIN(MLONG*FAKTOR)
        SBG=SBM*CI-CBM*CLM*SI
        LATI=ASIN(SBG)
        CBG=COS(LATI)
        SLG=(CBM*SLM)/CBG
        CLG=(SBM*SI+CBM*CLM*CI)/CBG
        IF(CLG.GT.1.) CLG=1.
        LONG=ACOS(CLG)
        IF(SLG.LT.0.0) LONG=ZPI-ACOS(CLG)
        LATI=LATI/FAKTOR
        LONG=LONG/FAKTOR
        LONG=LONG-69.8
        IF(LONG.LT.-180.0) LONG=LONG+360.0
        IF(LONG.GT. 180.0) LONG=LONG-360.0
      ELSE
        YLG=LONG+69.8
        CBG=COS(LATI*FAKTOR)
        SBG=SIN(LATI*FAKTOR)
        CLG=COS(YLG*FAKTOR)
        SLG=SIN(YLG*FAKTOR)
        SBM=SBG*CI+CBG*CLG*SI
        MLAT=ASIN(SBM)
        CBM=COS(MLAT)
        SLM=(CBG*SLG)/CBM
        CLM=(-SBG*SI+CBG*CLG*CI)/CBM
        IF(CLM.GT.1.) CLM=1.
        MLONG=ACOS(CLM)
        IF(SLM.LT..0) MLONG=ZPI-ACOS(CLM)
        MLAT=MLAT/FAKTOR
        MLONG=MLONG/FAKTOR
        IF(MLONG.LT.-180.0) MLONG=MLONG+360.0
        IF(MLONG.GT. 180.0) MLONG=MLONG-360.0
      ENDIF
      RETURN
    END SUBROUTINE GEO2MAG
