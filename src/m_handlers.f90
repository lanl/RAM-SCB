MODULE m_handlers

  USE FoX_common
  USE FoX_sax
  USE Module1, ONLY : isotropy, iReduceAnisotropy, iOuterMethod, iPressureChoice, method, &
       iWantAlphaExtrapolation, numit, decreaseConv, blendInitial, iSm, iSm2, nimax, &
       iAMR, iAzimOffset, isSorDetailNeeded, isEnergDetailNeeded, iLossCone, iOutput, iTilt, &
       iConvectionChoice, iHourChoice, iCompDomain, prefixIn, prefixOut, prefixPres
  IMPLICIT NONE
  PRIVATE

  ! This shows how you might use a SAX processor
  ! to extract character data from within certain nodes.

  ! It will print out the character content of all scalar
  ! nodes inside parameter nodes.

  ! It will also remember the value of the scalar within the 
  ! parameter node called DM.EnergyTolerance

  ! Note that when dealing with character data, we need to 
  ! concatenate all the character handling events in question:
  ! we do this using some of the character handling routines 
  ! provided as a convenience by FoX.

  PUBLIC :: characters_handler
  PUBLIC :: endElement_handler
  PUBLIC :: startDocument_handler
  PUBLIC :: startElement_handler

  LOGICAL :: inScalar
  LOGICAL :: inParameter
  LOGICAL :: isotropyFound
  ! REAL :: isotropy
  CHARACTER(200) :: c

  ! PUBLIC :: etol

CONTAINS

  SUBROUTINE characters_handler(chunk)
    CHARACTER(len=*), INTENT(in) :: chunk

    IF (inScalar.AND.inParameter) THEN
       PRINT*, "Found some scalar parameter data:"
       PRINT*, TRIM(chunk)
       IF (isotropyFound) THEN
          c = TRIM(c)//" "//chunk
       ENDIF
    ENDIF

  END SUBROUTINE characters_handler

  SUBROUTINE endElement_handler(URI, localname, name)
    CHARACTER(len=*), INTENT(in)     :: URI
    CHARACTER(len=*), INTENT(in)     :: localname
    CHARACTER(len=*), INTENT(in)     :: name

    ! Note if we are leaving a scalar or parameter element.
    IF (localName=="scalar") THEN
      inScalar = .FALSE.
      IF (isotropyFound) THEN
        ! pull the data out of the concatenated string:
        CALL rts(c, isotropy)
      ENDIF
    ELSEIF (localName=="parameter") THEN
      inParameter = .FALSE.
      isotropyFound = .FALSE.
    ENDIF
  END SUBROUTINE endElement_handler
  
  SUBROUTINE startDocument_handler
    ! Initialize state variables
    inScalar = .FALSE.
    inParameter = .FALSE.
    ! Initalize other module variables
    c = ''
  END SUBROUTINE startDocument_handler
  
  SUBROUTINE startElement_handler(URI, localname, name, attributes)
    CHARACTER(len=*), INTENT(in)   :: URI
    CHARACTER(len=*), INTENT(in)   :: localname
    CHARACTER(len=*), INTENT(in)   :: name
    TYPE(dictionary_t), INTENT(in) :: attributes

    INTEGER :: i

    ! Note if we are entering a scalar or parameter element.
    IF (localName=="scalar") THEN
      inScalar = .TRUE.
    ELSEIF (localName=="parameter") THEN
      inParameter = .TRUE.
      ! Loop over the attributes, looking for the name ...
      DO i = 1, getLength(attributes)
        IF (getQName(attributes, i)=="name") THEN
          ! And if the attribute is called name, check if the name
          ! is what we are interested in.
          PRINT*, "name=", getValue(attributes, i)
          isotropyFound = getValue(attributes, i)=="Isotropic"
         ! reduceAnisotropyFound = getValue(attributes, i)=="ReduceAnisotropy"
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE startElement_handler

END MODULE m_handlers
