!****************************************************************************
! Interface between RAM and 3D equilibrium code (computes h, I integrals, 
! as well as RAM flux along field lines)
! Author: S. Zaharia, 2006-2009
! Latest version: Nov. 2009
! Copyright (c) 2016, Los Alamos National Security, LLC
! All rights reserved.
!****************************************************************************


MODULE Module_RAM
  USE nrtype, ONLY : SP, DP
  USE Module1, ONLY : iHourChoice, iOutput, nthe, npsi, nzeta, bZ, rank, numProc, prefixOut, nXRaw, nAzimRAM, iDumpRAMFlux
  USE ModRamMain,  ONLY: PathRamIn, IsComponent, FLUX, NE, NPA, NT, RadiusMax
  USE ModRamCouple,ONLY: SwmfPot_II
  USE mpi

  IMPLICIT NONE

INTERFACE DQAGI
   SUBROUTINE DQAGI(F, BOUND, INF, EPSABS, EPSREL, RESULT, ABSERR, &
        NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
     USE nrtype, ONLY : DP
     REAL(DP), EXTERNAL :: F
     REAL(DP), INTENT(IN) :: BOUND
     INTEGER, INTENT(IN) :: INF
     REAL(DP), INTENT(IN) :: EPSABS, EPSREL
     REAL(DP), INTENT(OUT) :: RESULT, ABSERR
     INTEGER, INTENT(OUT) :: NEVAL, IER
     INTEGER, INTENT(IN) :: LIMIT, LENW
     INTEGER, INTENT(OUT) :: LAST
     INTEGER, INTENT(IN OUT) :: IWORK(LIMIT)
     REAL(DP), INTENT(IN OUT) :: WORK(LENW)
   END SUBROUTINE DQAGI
END INTERFACE

INTERFACE DQAG
   SUBROUTINE DQAG(F, A, B, EPSABS, EPSREL, KEY, RESULT, ABSERR, &
        NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
     USE nrtype, ONLY : DP
     REAL(DP), EXTERNAL :: F
     REAL(DP), INTENT(IN) :: A, B 
     REAL(DP), INTENT(IN) :: EPSABS, EPSREL
     REAL(DP), INTENT(OUT) :: RESULT, ABSERR
     INTEGER, INTENT(OUT) :: NEVAL, IER
     INTEGER, INTENT(IN) :: KEY, LIMIT, LENW
     INTEGER, INTENT(OUT) :: LAST
     INTEGER, INTENT(IN OUT) :: IWORK(LIMIT)
     REAL(DP), INTENT(IN OUT) :: WORK(LENW)
   END SUBROUTINE DQAG
END INTERFACE

  INTERFACE qromb
     FUNCTION qromb(func,a,b)
       USE nrtype; USE nrutil, ONLY : nrerror
       USE nr, ONLY : polint,trapzd
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP) :: qromb
       INTERFACE
          FUNCTION func(x)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(SIZE(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION qromb
  END INTERFACE

  INTERFACE splint_interface
     FUNCTION splint(xa,ya,y2a,x)
       USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
       USE nr, ONLY: locate
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: splint
     END FUNCTION splint
  END INTERFACE

  INTERFACE spline_interface
     SUBROUTINE spline(x,y,yp1,ypn,y2)
       USE nrtype; USE nrutil, ONLY : assert_eq
       USE nr, ONLY : tridag
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(DP), INTENT(IN) :: yp1,ypn
       REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
     END SUBROUTINE spline
  END INTERFACE

  INTERFACE polint
     SUBROUTINE polint(xa,ya,x,y,dy)
       USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y, dy
     END SUBROUTINE polint
  END INTERFACE

  INTERFACE hdens_rairden
     FUNCTION hdens_rairden(radius)
       USE nrtype, ONLY : DP
       REAL(DP), INTENT(IN) :: radius
       REAL(DP) :: hdens_rairden
     END FUNCTION hdens_rairden
  END INTERFACE

 INTERFACE Spline_2D_point
     SUBROUTINE Spline_2D_point(x_1, x_2, f, x, y, func, ierrDomain)
       USE EZspline_obj ! import the modules
       USE EZspline  
       USE nrtype, ONLY : DP
       IMPLICIT NONE
       INTEGER, PARAMETER :: r8 = DP
       REAL(r8), DIMENSION(:), INTENT(IN) :: x_1, x_2 ! independent variable
       REAL(r8), DIMENSION(:,:), INTENT(IN) :: f
       REAL(r8), DIMENSION(:,:), INTENT(OUT) :: func   ! interpolated values 
       REAL(r8), DIMENSION(:,:), INTENT(IN)    :: x, y  ! grid of points for output
       INTEGER, INTENT(out) :: ierrDomain
     END SUBROUTINE Spline_2D_point
  END INTERFACE

  INTEGER, PARAMETER :: NRAD = nXRaw+1
  INTEGER :: INDEXPA(nthe,npsi,nzeta,NPA) ! Index that maps the pitch angle on each point of a field line to the equatorial (RAM) one
  INTEGER :: INCFD, j, k, L
  REAL(DP) :: bfMirror(NPA), bfInterm(NPA)
  REAL(DP) :: chiMirror(NPA)
  REAL(DP) :: derivs2(1:nthe), dBdTheta(1:nthe), dThetadB(1:nthe), derivs4(1:nthe), distance(nthe), dyDummyDens(nthe)
  REAL(DP) :: flux3DEQ(4, npsi,nzeta,NE,NPA)
END MODULE Module_RAM

