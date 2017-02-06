FUNCTION fr(r)
  !cc   function that defines the flux surface boundary
  !--------------------
  !      use clich1
  !     use clich6
  !--------------------

  USE nrtype
  USE Module1

  IMPLICIT NONE

  REAL(DP), INTENT(IN)                         :: r
  REAL(DP)                                     :: fr
  REAL(DP)                                     :: pdip, pcom, pimf
  REAL(DP)                                     :: bi, bc2, bs3, psibc, thet, xcent

  COMMON/param0/bi,bc2,bs3,psibc,thet,xcent
  pdip=-xzero3*COS(thet)**2/r
  pimf=0.5*bi*r**2*COS(thet)**2
  pcom=bs3*(r-re)**2*SIN(3.*thet)**2
  !      pcom=pcom+bc2*(r-re)*cos(2.*thet)
  pcom=pcom+bc2*(r-re)*(COS(2.*thet)+2.*COS(thet))/3.
  fr=1.0-(pdip+pimf+pcom)/psibc
  RETURN
END FUNCTION fr
