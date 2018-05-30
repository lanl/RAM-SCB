!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModScbEquation
! Contains subroutines for calculating the left hand side and
! right hand side of the force balance equation


  implicit none

contains

!====================================================!
!================ LHS EQUATIONS =====================!
!====================================================!
!==============================================================================
SUBROUTINE metrica

  USE ModRamParams,    ONLY: verbose
  USE ModScbMain,      ONLY: DP
  USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, ny, rdr, &
                             rdt, rdp, rdt4, rdp4, rdr2, rdp2, rdt2, rdpsq, &
                             rdtsq, rdpdt4, rdr4, rdrsq
  use ModScbVariables, ONLY: f, fzet, x, y, z, sumb, sumdb, bf, bsq, jacobian, &
                             thetaVal, rhoVal, zetaVal, vecd, vec1, vec2, vec3, &
                             vec4, vec6, vec7, vec8, vec9, left, right

  use ModRamGSL, ONLY: GSL_Derivs


  implicit none

  INTEGER :: i, j, k, GSLerr, nl0

  REAL(DP) :: xa, xb, xc, xd, xe, xta, xtb, xtc, xtd, xte, xpa, xpb, xpc, xpd, &
              xpe, xra, xrb, xrc, xrd, xre, yta, ytb, ytc, ytd, yte, ypa, ypb, &
              ypc, ypd, ype, yra, yrb, yrc, yrd, yre, zta, ztb, ztc, ztd, zte, &
              zpa, zpb, zpc, zpd, zpe, zra, zrb, zrc, zrd, zre, aja, ajb, ajc, &
              ajd, aje, grxa, grya, grza, grxb, gryb, grzb, grxc, gryc, grzc, &
              grxd, gryd, grzd, grxe, grye, grze, gpxa, gpya, gpza, gpxb, gpyb, &
              gpzb, gpxc, gpyc, gpzc, gpxd, gpyd, gpzd, gpxe, gpye, gpze, gtxa, &
              gtya, gtza, gtxb, gtyb, gtzb, gtxc, gtyc, gtzc, gtxd, gtyd, gtzd, &
              gtxe, gtye, gtze, grsa, grsb, grsc, grsd, grse, gpsa, gpsb, gpsc, &
              gpsd, gpse, gtsa, gtsb, gtsc, gtsd, gtse, grgpa, grgpb, grgpc, &
              grgpd, grgpe, gpgta, gpgtb, gpgtc, gpgtd, gpgte, gtgra, gtgrb, &
              gtgrc, gtgrd, gtgre, v1a, v1b, v1c, v2a, v2b, v2c, v2d, v3a, v3b, &
              v3c, v3d, bsqa, bsqb, bsqc, bsqd

  !jacobian = 0._dp
  vecd = 0._dp
  vec1 = 0._dp
  vec2 = 0._dp
  vec3 = 0._dp
  vec4 = 0._dp
  vec6 = 0._dp
  vec7 = 0._dp
  vec8 = 0._dp
  vec9 = 0._dp
  nl0 = 0

  DO j=left+1,right-1
     DO k=2,nzeta
        DO i=2,nthem

           xa = .5*(x(i,j,k)+x(i+1,j,k))
           xb = .5*(x(i,j,k)+x(i,j,k+1))
           xc = .5*(x(i,j,k)+x(i-1,j,k))
           xd = .5*(x(i,j,k)+x(i,j,k-1))
           xe = x(i,j,k)

           xta=(x(i+1,j,k)-x(i,j,k))*rdt
           xtb=(x(i+1,j,k+1)+x(i+1,j,k)-x(i-1,j,k+1)-x(i-1,j,k))*rdt4
           xtc=(x(i,j,k)-x(i-1,j,k))*rdt
           xtd=(x(i+1,j,k)+x(i+1,j,k-1)-x(i-1,j,k)-x(i-1,j,k-1))*rdt4
           xte=(x(i+1,j,k)-x(i-1,j,k))*rdt2
           
           xpa=(x(i+1,j,k+1)+x(i,j,k+1)-x(i+1,j,k-1)-x(i,j,k-1))*rdp4
           xpb=(x(i,j,k+1)-x(i,j,k))*rdp
           xpc=(x(i,j,k+1)+x(i-1,j,k+1)-x(i,j,k-1)-x(i-1,j,k-1))*rdp4
           xpd=(x(i,j,k)-x(i,j,k-1))*rdp
           xpe=(x(i,j,k+1)-x(i,j,k-1))*rdp2
           
           xra=(x(i+1,j+1,k)+x(i,j+1,k)-x(i+1,j-1,k)-x(i,j-1,k))*rdr4
           xrb=(x(i,j+1,k+1)+x(i,j+1,k)-x(i,j-1,k+1)-x(i,j-1,k))*rdr4
           xrc=(x(i,j+1,k)+x(i-1,j+1,k)-x(i,j-1,k)-x(i-1,j-1,k))*rdr4
           xrd=(x(i,j+1,k-1)+x(i,j+1,k)-x(i,j-1,k-1)-x(i,j-1,k))*rdr4
           xre=(x(i,j+1,k)-x(i,j-1,k))*rdr2

           yta=(y(i+1,j,k)-y(i,j,k))*rdt
           ytb=(y(i+1,j,k+1)+y(i+1,j,k)-y(i-1,j,k+1)-y(i-1,j,k))*rdt4
           ytc=(y(i,j,k)-y(i-1,j,k))*rdt
           ytd=(y(i+1,j,k)+y(i+1,j,k-1)-y(i-1,j,k)-y(i-1,j,k-1))*rdt4
           yte=(y(i+1,j,k)-y(i-1,j,k))*rdt2
           
           ypa=(y(i+1,j,k+1)+y(i,j,k+1)-y(i+1,j,k-1)-y(i,j,k-1))*rdp4
           ypb=(y(i,j,k+1)-y(i,j,k))*rdp
           ypc=(y(i,j,k+1)+y(i-1,j,k+1)-y(i,j,k-1)-y(i-1,j,k-1))*rdp4
           ypd=(y(i,j,k)-y(i,j,k-1))*rdp
           ype=(y(i,j,k+1)-y(i,j,k-1))*rdp2
           
           yra=(y(i+1,j+1,k)+y(i,j+1,k)-y(i+1,j-1,k)-y(i,j-1,k))*rdr4
           yrb=(y(i,j+1,k+1)+y(i,j+1,k)-y(i,j-1,k+1)-y(i,j-1,k))*rdr4
           yrc=(y(i,j+1,k)+y(i-1,j+1,k)-y(i,j-1,k)-y(i-1,j-1,k))*rdr4
           yrd=(y(i,j+1,k-1)+y(i,j+1,k)-y(i,j-1,k-1)-y(i,j-1,k))*rdr4
           yre=(y(i,j+1,k)-y(i,j-1,k))*rdr2

           zta=(z(i+1,j,k)-z(i,j,k))*rdt
           ztb=(z(i+1,j,k+1)+z(i+1,j,k)-z(i-1,j,k+1)-z(i-1,j,k))*rdt4
           ztc=(z(i,j,k)-z(i-1,j,k))*rdt
           ztd=(z(i+1,j,k)+z(i+1,j,k-1)-z(i-1,j,k)-z(i-1,j,k-1))*rdt4
           zte=(z(i+1,j,k)-z(i-1,j,k))*rdt2
           
           zpa=(z(i+1,j,k+1)+z(i,j,k+1)-z(i+1,j,k-1)-z(i,j,k-1))*rdp4
           zpb=(z(i,j,k+1)-z(i,j,k))*rdp
           zpc=(z(i,j,k+1)+z(i-1,j,k+1)-z(i,j,k-1)-z(i-1,j,k-1))*rdp4
           zpd=(z(i,j,k)-z(i,j,k-1))*rdp
           zpe=(z(i,j,k+1)-z(i,j,k-1))*rdp2
           
           zra=(z(i+1,j+1,k)+z(i,j+1,k)-z(i+1,j-1,k)-z(i,j-1,k))*rdr4
           zrb=(z(i,j+1,k+1)+z(i,j+1,k)-z(i,j-1,k+1)-z(i,j-1,k))*rdr4
           zrc=(z(i,j+1,k)+z(i-1,j+1,k)-z(i,j-1,k)-z(i-1,j-1,k))*rdr4
           zrd=(z(i,j+1,k-1)+z(i,j+1,k)-z(i,j-1,k-1)-z(i,j-1,k))*rdr4
           zre=(z(i,j+1,k)-z(i,j-1,k))*rdr2

           aja = xra*(ypa*zta-yta*zpa)  &
                +xpa*(yta*zra-yra*zta)  &
                +xta*(yra*zpa-ypa*zra)
           ajb = xrb*(ypb*ztb-ytb*zpb)  &
                +xpb*(ytb*zrb-yrb*ztb)  &
                +xtb*(yrb*zpb-ypb*zrb)
           ajc = xrc*(ypc*ztc-ytc*zpc)  &
                +xpc*(ytc*zrc-yrc*ztc)  &
                +xtc*(yrc*zpc-ypc*zrc)
           ajd = xrd*(ypd*ztd-ytd*zpd)  &
                +xpd*(ytd*zrd-yrd*ztd)  &
                +xtd*(yrd*zpd-ypd*zrd)
           aje = xre*(ype*zte-yte*zpe)  &
                +xpe*(yte*zre-yre*zte)  &
                +xte*(yre*zpe-ype*zre)

           IF (aje < 0._dp) THEN
              jacobian(i,j,k) = aje
              if (verbose) PRINT*, 'metrica: J < 0; i, j, k = ', i, j, k
              if (verbose) write(*,*) jacobian(i,j,k)
              return
           END IF

           grxa = (ypa*zta-yta*zpa)/aja
           grya = (zpa*xta-zta*xpa)/aja
           grza = (xpa*yta-xta*ypa)/aja
           grxb = (ypb*ztb-ytb*zpb)/ajb
           gryb = (zpb*xtb-ztb*xpb)/ajb
           grzb = (xpb*ytb-xtb*ypb)/ajb
           grxc = (ypc*ztc-ytc*zpc)/ajc
           gryc = (zpc*xtc-ztc*xpc)/ajc
           grzc = (xpc*ytc-xtc*ypc)/ajc
           grxd = (ypd*ztd-ytd*zpd)/ajd
           gryd = (zpd*xtd-ztd*xpd)/ajd
           grzd = (xpd*ytd-xtd*ypd)/ajd
           grxe = (ype*zte-yte*zpe)/aje
           grye = (zpe*xte-zte*xpe)/aje
           grze = (xpe*yte-xte*ype)/aje


           gpxa = (yta*zra-yra*zta)/aja
           gpya = (zta*xra-zra*xta)/aja
           gpza = (xta*yra-xra*yta)/aja
           gpxb = (ytb*zrb-yrb*ztb)/ajb
           gpyb = (ztb*xrb-zrb*xtb)/ajb
           gpzb = (xtb*yrb-xrb*ytb)/ajb
           gpxc = (ytc*zrc-yrc*ztc)/ajc
           gpyc = (ztc*xrc-zrc*xtc)/ajc
           gpzc = (xtc*yrc-xrc*ytc)/ajc
           gpxd = (ytd*zrd-yrd*ztd)/ajd
           gpyd = (ztd*xrd-zrd*xtd)/ajd
           gpzd = (xtd*yrd-xrd*ytd)/ajd
           gpxe = (yte*zre-yre*zte)/aje
           gpye = (zte*xre-zre*xte)/aje
           gpze = (xte*yre-xre*yte)/aje

           gtxa = (yra*zpa-ypa*zra)/aja
           gtya = (zra*xpa-zpa*xra)/aja
           gtza = (xra*ypa-xpa*yra)/aja
           gtxb = (yrb*zpb-ypb*zrb)/ajb
           gtyb = (zrb*xpb-zpb*xrb)/ajb
           gtzb = (xrb*ypb-xpb*yrb)/ajb
           gtxc = (yrc*zpc-ypc*zrc)/ajc
           gtyc = (zrc*xpc-zpc*xrc)/ajc
           gtzc = (xrc*ypc-xpc*yrc)/ajc
           gtxd = (yrd*zpd-ypd*zrd)/ajd
           gtyd = (zrd*xpd-zpd*xrd)/ajd
           gtzd = (xrd*ypd-xpd*yrd)/ajd
           gtxe = (yre*zpe-ype*zre)/aje
           gtye = (zre*xpe-zpe*xre)/aje
           gtze = (xre*ype-xpe*yre)/aje

           grsa = ( grxa**2+grya**2+grza**2 )
           grsb = ( grxb**2+gryb**2+grzb**2 )
           grsc = ( grxc**2+gryc**2+grzc**2 )
           grsd = ( grxd**2+gryd**2+grzd**2 )
           grse = ( grxe**2+grye**2+grze**2 )
           gpsa = ( gpxa**2+gpya**2+gpza**2 )
           gpsb = ( gpxb**2+gpyb**2+gpzb**2 )
           gpsc = ( gpxc**2+gpyc**2+gpzc**2 )
           gpsd = ( gpxd**2+gpyd**2+gpzd**2 )
           gpse = ( gpxe**2+gpye**2+gpze**2 )
           gtsa = ( gtxa**2+gtya**2+gtza**2 )
           gtsb = ( gtxb**2+gtyb**2+gtzb**2 )
           gtsc = ( gtxc**2+gtyc**2+gtzc**2 )
           gtsd = ( gtxd**2+gtyd**2+gtzd**2 )
           gtse = ( gtxe**2+gtye**2+gtze**2 )

           grgpa = (gpxa*grxa+gpya*grya+gpza*grza)
           grgpb = (gpxb*grxb+gpyb*gryb+gpzb*grzb)
           grgpc = (gpxc*grxc+gpyc*gryc+gpzc*grzc)
           grgpd = (gpxd*grxd+gpyd*gryd+gpzd*grzd)
           grgpe = (gpxe*grxe+gpye*grye+gpze*grze)

           gpgta = (gpxa*gtxa+gpya*gtya+gpza*gtza)
           gpgtb = (gpxb*gtxb+gpyb*gtyb+gpzb*gtzb)
           gpgtc = (gpxc*gtxc+gpyc*gtyc+gpzc*gtzc)
           gpgtd = (gpxd*gtxd+gpyd*gtyd+gpzd*gtzd)
           gpgte = (gpxe*gtxe+gpye*gtye+gpze*gtze)

           gtgra = (gtxa*grxa+gtya*grya+gtza*grza)
           gtgrb = (gtxb*grxb+gtyb*gryb+gtzb*grzb)
           gtgrc = (gtxc*grxc+gtyc*gryc+gtzc*grzc)
           gtgrd = (gtxd*grxd+gtyd*gryd+gtzd*grzd)
           gtgre = (gtxe*grxe+gtye*grye+gtze*grze)

           v1a=(grsa*gtsa-gtgra**2)*aja*rdtsq
           v1b=(grsb*gtsb-gtgrb**2)*ajb*rdtsq
           v1c=(grsc*gtsc-gtgrc**2)*ajc*rdtsq
           v2a=(grsa*gpgta-grgpa*gtgra)*aja*rdpdt4
           v2b=(grsb*gpgtb-grgpb*gtgrb)*ajb*rdpdt4
           v2c=(grsc*gpgtc-grgpc*gtgrc)*ajc*rdpdt4
           v2d=(grsd*gpgtd-grgpd*gtgrd)*ajd*rdpdt4
           v3b=(grsb*gpsb-grgpb**2)*ajb*rdpsq
           v3c=(grsc*gpsc-grgpc**2)*ajc*rdpsq
           v3d=(grsd*gpsd-grgpd**2)*ajd*rdpsq

           vecd(i,j,k) = (v1a+v1c) + (v3b+v3d)
           vec1(i,j,k) = (v2c+v2d)
           vec2(i,j,k) = (v2c-v2a) + v3d
           vec3(i,j,k) = -(v2a+v2d)
           vec4(i,j,k) = v1c + (v2d-v2b)
           vec6(i,j,k) = v1a + (v2b-v2d)
           vec7(i,j,k) = -(v2c+v2b)
           vec8(i,j,k) = v3b + (v2a-v2c)
           vec9(i,j,k) = (v2a + v2b)

        END DO
     END DO
  END DO

  RETURN

END SUBROUTINE metrica

!------------------------------------------------------------------------------
SUBROUTINE metric

  USE ModRamParams,    ONLY: verbose
  USE ModScbMain,      ONLY: DP
  USE ModScbGrids,     ONLY: nthe, nthem, npsi, npsim, nzeta, nzetap, ny, rdr, &
                             rdt, rdp, rdt4, rdp4, rdr2, rdr4, rdp2, rdt2, rdpsq, &
                             rdtsq, rdpdt4, rdrsq, rdtdr4
  use ModScbVariables, ONLY: f, fzet, x, y, z, sumb, sumdb, bf, bsq, jacobian, &
                             thetaVal, rhoVal, zetaVal, vecd, vec1, vec2, vec3, &
                             vec4, vec6, vec7, vec8, vec9, left, right

  use ModRamGSL, ONLY: GSL_Derivs


  implicit none

  INTEGER :: i, j, k, GSLerr, nl0

  REAL(DP) :: xa, xb, xc, xd, xe, xta, xtb, xtc, xtd, xte, xpa, xpb, xpc, xpd, &
              xpe, xra, xrb, xrc, xrd, xre, yta, ytb, ytc, ytd, yte, ypa, ypb, &
              ypc, ypd, ype, yra, yrb, yrc, yrd, yre, zta, ztb, ztc, ztd, zte, &
              zpa, zpb, zpc, zpd, zpe, zra, zrb, zrc, zrd, zre, aja, ajb, ajc, &
              ajd, aje, grxa, grya, grza, grxb, gryb, grzb, grxc, gryc, grzc, &
              grxd, gryd, grzd, grxe, grye, grze, gpxa, gpya, gpza, gpxb, gpyb, &
              gpzb, gpxc, gpyc, gpzc, gpxd, gpyd, gpzd, gpxe, gpye, gpze, gtxa, &
              gtya, gtza, gtxb, gtyb, gtzb, gtxc, gtyc, gtzc, gtxd, gtyd, gtzd, &
              gtxe, gtye, gtze, grsa, grsb, grsc, grsd, grse, gpsa, gpsb, gpsc, &
              gpsd, gpse, gtsa, gtsb, gtsc, gtsd, gtse, grgpa, grgpb, grgpc, &
              grgpd, grgpe, gpgta, gpgtb, gpgtc, gpgtd, gpgte, gtgra, gtgrb, &
              gtgrc, gtgrd, gtgre, v1a, v1b, v1c, v2a, v2b, v2c, v2d, v3a, v3b, &
              v3c, v3d, bsqa, bsqb, bsqc, bsqd

  !jacobian = 0._dp
  vecd = 0._dp
  vec1 = 0._dp
  vec2 = 0._dp
  vec3 = 0._dp
  vec4 = 0._dp
  vec6 = 0._dp
  vec7 = 0._dp
  vec8 = 0._dp
  vec9 = 0._dp
  nl0 = 0

  DO j=left+1,right-1
     DO k=2,nzeta
        DO i=2,nthem

           xa = .5*(x(i,j,k)+x(i+1,j,k))
           xb = .5*(x(i,j,k)+x(i,j+1,k))
           xc = .5*(x(i,j,k)+x(i-1,j,k))
           xd = .5*(x(i,j,k)+x(i,j-1,k))
           xe = x(i,j,k)

           xta=(x(i+1,j,k)-x(i,j,k))*rdt
           xtb=(x(i+1,j+1,k)+x(i+1,j,k)-x(i-1,j+1,k)-x(i-1,j,k))*rdt4
           xtc=(x(i,j,k)-x(i-1,j,k))*rdt
           xtd=(x(i+1,j,k)+x(i+1,j-1,k)-x(i-1,j,k)-x(i-1,j-1,k))*rdt4
           xte=(x(i+1,j,k)-x(i-1,j,k))*rdt2
           xpa=(x(i+1,j,k+1)+x(i,j,k+1)-x(i+1,j,k-1)-x(i,j,k-1))*rdp4
           xpb=(x(i,j+1,k+1)+x(i,j,k+1)-x(i,j+1,k-1)-x(i,j,k-1))*rdp4
           xpc=(x(i,j,k+1)+x(i-1,j,k+1)-x(i,j,k-1)-x(i-1,j,k-1))*rdp4
           xpd=(x(i,j,k+1)+x(i,j-1,k+1)-x(i,j,k-1)-x(i,j-1,k-1))*rdp4
           xpe=(x(i,j,k+1)-x(i,j,k-1))*rdp2

           xra=(x(i+1,j+1,k)+x(i,j+1,k)-x(i+1,j-1,k)-x(i,j-1,k))*rdr4
           xrb=(x(i,j+1,k)-x(i,j,k))*rdr
           xrc=(x(i,j+1,k)+x(i-1,j+1,k)-x(i,j-1,k)-x(i-1,j-1,k))*rdr4
           xrd=(x(i,j,k)-x(i,j-1,k))*rdr
           xre=(x(i,j+1,k)-x(i,j-1,k))*rdr2

           yta=(y(i+1,j,k)-y(i,j,k))*rdt
           ytb=(y(i+1,j+1,k)+y(i+1,j,k)-y(i-1,j+1,k)-y(i-1,j,k))*rdt4
           ytc=(y(i,j,k)-y(i-1,j,k))*rdt
           ytd=(y(i+1,j,k)+y(i+1,j-1,k)-y(i-1,j,k)-y(i-1,j-1,k))*rdt4
           yte=(y(i+1,j,k)-y(i-1,j,k))*rdt2
           ypa=(y(i+1,j,k+1)+y(i,j,k+1)-y(i+1,j,k-1)-y(i,j,k-1))*rdp4
           ypb=(y(i,j+1,k+1)+y(i,j,k+1)-y(i,j+1,k-1)-y(i,j,k-1))*rdp4
           ypc=(y(i,j,k+1)+y(i-1,j,k+1)-y(i,j,k-1)-y(i-1,j,k-1))*rdp4
           ypd=(y(i,j,k+1)+y(i,j-1,k+1)-y(i,j,k-1)-y(i,j-1,k-1))*rdp4
           ype=(y(i,j,k+1)-y(i,j,k-1))*rdp2
           yra=(y(i+1,j+1,k)+y(i,j+1,k)-y(i+1,j-1,k)-y(i,j-1,k))*rdr4
           yrb=(y(i,j+1,k)-y(i,j,k))*rdr
           yrc=(y(i,j+1,k)+y(i-1,j+1,k)-y(i,j-1,k)-y(i-1,j-1,k))*rdr4
           yrd=(y(i,j,k)-y(i,j-1,k))*rdr
           yre=(y(i,j+1,k)-y(i,j-1,k))*rdr2

           zta=(z(i+1,j,k)-z(i,j,k))*rdt
           ztb=(z(i+1,j+1,k)+z(i+1,j,k)-z(i-1,j+1,k)-z(i-1,j,k))*rdt4
           ztc=(z(i,j,k)-z(i-1,j,k))*rdt
           ztd=(z(i+1,j,k)+z(i+1,j-1,k)-z(i-1,j,k)-z(i-1,j-1,k))*rdt4
           zte=(z(i+1,j,k)-z(i-1,j,k))*rdt2
           zpa=(z(i+1,j,k+1)+z(i,j,k+1)-z(i+1,j,k-1)-z(i,j,k-1))*rdp4
           zpb=(z(i,j+1,k+1)+z(i,j,k+1)-z(i,j+1,k-1)-z(i,j,k-1))*rdp4
           zpc=(z(i,j,k+1)+z(i-1,j,k+1)-z(i,j,k-1)-z(i-1,j,k-1))*rdp4
           zpd=(z(i,j,k+1)+z(i,j-1,k+1)-z(i,j,k-1)-z(i,j-1,k-1))*rdp4
           zpe=(z(i,j,k+1)-z(i,j,k-1))*rdp2
           zra=(z(i+1,j+1,k)+z(i,j+1,k)-z(i+1,j-1,k)-z(i,j-1,k))*rdr4
           zrb=(z(i,j+1,k)-z(i,j,k))*rdr
           zrc=(z(i,j+1,k)+z(i-1,j+1,k)-z(i,j-1,k)-z(i-1,j-1,k))*rdr4
           zrd=(z(i,j,k)-z(i,j-1,k))*rdr
           zre=(z(i,j+1,k)-z(i,j-1,k))*rdr2

           aja = xra*(ypa*zta-yta*zpa)  &
                +xpa*(yta*zra-yra*zta)  &
                +xta*(yra*zpa-ypa*zra)
           ajb = xrb*(ypb*ztb-ytb*zpb)  &
                +xpb*(ytb*zrb-yrb*ztb)  &
                +xtb*(yrb*zpb-ypb*zrb)
           ajc = xrc*(ypc*ztc-ytc*zpc)  &
                +xpc*(ytc*zrc-yrc*ztc)  &
                +xtc*(yrc*zpc-ypc*zrc)
           ajd = xrd*(ypd*ztd-ytd*zpd)  &
                +xpd*(ytd*zrd-yrd*ztd)  &
                +xtd*(yrd*zpd-ypd*zrd)
           aje = xre*(ype*zte-yte*zpe)  &
                +xpe*(yte*zre-yre*zte)  &
                +xte*(yre*zpe-ype*zre)

           IF (aje < 0._dp) THEN
              jacobian(i,j,k) = aje
              if (verbose) PRINT*, 'metric: J < 0; i, j, k = ', i, j, k
              if (verbose) PRINT*, jacobian(i,j,k)
              return
           END IF

           grxa = (ypa*zta-yta*zpa)/aja
           grya = (zpa*xta-zta*xpa)/aja
           grza = (xpa*yta-xta*ypa)/aja
           grxb = (ypb*ztb-ytb*zpb)/ajb
           gryb = (zpb*xtb-ztb*xpb)/ajb
           grzb = (xpb*ytb-xtb*ypb)/ajb
           grxc = (ypc*ztc-ytc*zpc)/ajc
           gryc = (zpc*xtc-ztc*xpc)/ajc
           grzc = (xpc*ytc-xtc*ypc)/ajc
           grxd = (ypd*ztd-ytd*zpd)/ajd
           gryd = (zpd*xtd-ztd*xpd)/ajd
           grzd = (xpd*ytd-xtd*ypd)/ajd
           grxe = (ype*zte-yte*zpe)/aje
           grye = (zpe*xte-zte*xpe)/aje
           grze = (xpe*yte-xte*ype)/aje

           gpxa = (yta*zra-yra*zta)/aja
           gpya = (zta*xra-zra*xta)/aja
           gpza = (xta*yra-xra*yta)/aja
           gpxb = (ytb*zrb-yrb*ztb)/ajb
           gpyb = (ztb*xrb-zrb*xtb)/ajb
           gpzb = (xtb*yrb-xrb*ytb)/ajb
           gpxc = (ytc*zrc-yrc*ztc)/ajc
           gpyc = (ztc*xrc-zrc*xtc)/ajc
           gpzc = (xtc*yrc-xrc*ytc)/ajc
           gpxd = (ytd*zrd-yrd*ztd)/ajd
           gpyd = (ztd*xrd-zrd*xtd)/ajd
           gpzd = (xtd*yrd-xrd*ytd)/ajd
           gpxe = (yte*zre-yre*zte)/aje
           gpye = (zte*xre-zre*xte)/aje
           gpze = (xte*yre-xre*yte)/aje

           gtxa = (yra*zpa-ypa*zra)/aja
           gtya = (zra*xpa-zpa*xra)/aja
           gtza = (xra*ypa-xpa*yra)/aja
           gtxb = (yrb*zpb-ypb*zrb)/ajb
           gtyb = (zrb*xpb-zpb*xrb)/ajb
           gtzb = (xrb*ypb-xpb*yrb)/ajb
           gtxc = (yrc*zpc-ypc*zrc)/ajc
           gtyc = (zrc*xpc-zpc*xrc)/ajc
           gtzc = (xrc*ypc-xpc*yrc)/ajc
           gtxd = (yrd*zpd-ypd*zrd)/ajd
           gtyd = (zrd*xpd-zpd*xrd)/ajd
           gtzd = (xrd*ypd-xpd*yrd)/ajd
           gtxe = (yre*zpe-ype*zre)/aje
           gtye = (zre*xpe-zpe*xre)/aje
           gtze = (xre*ype-xpe*yre)/aje

           grsa = ( grxa**2+grya**2+grza**2 )
           grsb = ( grxb**2+gryb**2+grzb**2 )
           grsc = ( grxc**2+gryc**2+grzc**2 )
           grsd = ( grxd**2+gryd**2+grzd**2 )
           grse = ( grxe**2+grye**2+grze**2 )
           gpsa = ( gpxa**2+gpya**2+gpza**2 )
           gpsb = ( gpxb**2+gpyb**2+gpzb**2 )
           gpsc = ( gpxc**2+gpyc**2+gpzc**2 )
           gpsd = ( gpxd**2+gpyd**2+gpzd**2 )
           gpse = ( gpxe**2+gpye**2+gpze**2 )
           gtsa = ( gtxa**2+gtya**2+gtza**2 )
           gtsb = ( gtxb**2+gtyb**2+gtzb**2 )
           gtsc = ( gtxc**2+gtyc**2+gtzc**2 )
           gtsd = ( gtxd**2+gtyd**2+gtzd**2 )
           gtse = ( gtxe**2+gtye**2+gtze**2 )

           grgpa = (gpxa*grxa+gpya*grya+gpza*grza)
           grgpb = (gpxb*grxb+gpyb*gryb+gpzb*grzb)
           grgpc = (gpxc*grxc+gpyc*gryc+gpzc*grzc)
           grgpd = (gpxd*grxd+gpyd*gryd+gpzd*grzd)
           grgpe = (gpxe*grxe+gpye*grye+gpze*grze)

           gpgta = (gpxa*gtxa+gpya*gtya+gpza*gtza)
           gpgtb = (gpxb*gtxb+gpyb*gtyb+gpzb*gtzb)
           gpgtc = (gpxc*gtxc+gpyc*gtyc+gpzc*gtzc)
           gpgtd = (gpxd*gtxd+gpyd*gtyd+gpzd*gtzd)
           gpgte = (gpxe*gtxe+gpye*gtye+gpze*gtze)

           gtgra = (gtxa*grxa+gtya*grya+gtza*grza)
           gtgrb = (gtxb*grxb+gtyb*gryb+gtzb*grzb)
           gtgrc = (gtxc*grxc+gtyc*gryc+gtzc*grzc)
           gtgrd = (gtxd*grxd+gtyd*gryd+gtzd*grzd)
           gtgre = (gtxe*grxe+gtye*grye+gtze*grze)

           v1a=(gpgta**2-gpsa*gtsa)*aja*rdtsq
           v1b=(gpgtb**2-gpsb*gtsb)*ajb*rdtsq
           v1c=(gpgtc**2-gpsc*gtsc)*ajc*rdtsq
           v2a=(grgpa*gpgta-gpsa*gtgra)*aja*rdtdr4
           v2b=(grgpb*gpgtb-gpsb*gtgrb)*ajb*rdtdr4
           v2c=(grgpc*gpgtc-gpsc*gtgrc)*ajc*rdtdr4
           v2d=(grgpd*gpgtd-gpsd*gtgrd)*ajd*rdtdr4
           v3b=(grgpb**2-grsb*gpsb)*ajb*rdrsq
           v3c=(grgpc**2-grsc*gpsc)*ajc*rdrsq
           v3d=(grgpd**2-grsd*gpsd)*ajd*rdrsq

           vecd(i,j,k) = (v1a+v1c+v3b+v3d)
           vec1(i,j,k) = (v2c+v2d)
           vec2(i,j,k) = ((v2c-v2a)+v3d)
           vec3(i,j,k) = -(v2a+v2d)
           vec4(i,j,k) = (v1c+(v2d-v2b))
           vec6(i,j,k) = (v1a+(v2b-v2d))
           vec7(i,j,k) = -(v2c+v2b)
           vec8(i,j,k) = (v3b+(v2a-v2c))
           vec9(i,j,k) = (v2b+v2a)
        END DO
     END DO
  END DO

  RETURN

END SUBROUTINE metric

!====================================================!
!================ RHS EQUATIONS =====================!
!====================================================!
!==============================================================================
SUBROUTINE newk
  !     define right hand side of alpha equation  

  USE ModScbMain,      ONLY: DP
  use ModScbParams,    ONLY: iOuterMethod, isotropy
  USE ModScbGrids,     ONLY: nthe, npsi, nzeta
  use ModScbVariables, ONLY: jacobian, bsq, f, fzet, alphaVal, GradRhoGradTheta, &
                             GradRhoGradZeta, GradThetaGradZeta, GradRhoSq, &
                             GradZetaSq, GradThetaSq, dPdPsi, dSqPdAlphaSq, &
                             dPPerdAlpha, dBsqdAlpha, dBsqdTheta, dPPerdTheta, &
                             dBBdAlpha, dPdAlpha, sigma, dBsqdRho, dPPerdRho, &
                             GradRhoSq, dPPerdZeta, dBsqdZeta, vecx, vecd, &
                             left, right


  implicit none

  INTEGER         :: i, j, k
  REAL(DP) :: xpz, xpt, c0, tz, tt

  IF (isotropy == 1) THEN      ! Isotropic case
     DO j = left, right
        IF (iOuterMethod == 2) THEN   ! Newton method
           DO k = 2, nzeta
              vecx(1:nthe, j, k) = -jacobian(1:nthe,j,k)/f(j)**2 * (dpdAlpha(1:nthe,j,k) - dSqpdAlphaSq(1:nthe,j,k) * &
                                     alphaVal(k))
              vecd(1:nthe, j, k) = vecd(1:nthe, j, k) - jacobian(1:nthe,j,k)/f(j)**2 * dSqpdAlphaSq(1:nthe,j,k)
           END DO
        ELSE                          ! Picard method
           vecx(1:nthe, j, 1:nzeta) = - dpdAlpha(1:nthe, j, 1:nzeta) * jacobian(1:nthe, j, 1:nzeta) / f(j)**2
        END IF
     END DO

  ELSE                    ! Anisotropic case
     DO i = 1, nthe
        DO j = left, right
           DO k = 1, nzeta
              IF (iOuterMethod == 1) THEN ! Picard method
                 xpz = (GradRhoSq(i,j,k)*GradZetaSq(i,j,k)-GradRhoGradZeta(i,j,k)**2)
                 xpt = (GradRhoSq(i,j,k)*GradThetaGradZeta(i,j,k) &
                        - GradRhoGradZeta(i,j,k)*GradRhoGradTheta(i,j,k))
                 c0  = -(f(j)**2*fzet(k))/sigma(i,j,k)/bsq(i,j,k)
                 tz  = dPPerdZeta(i,j,k) + 0.5*(1.-sigma(i,j,k))*dBsqdZeta(i,j,k)
                 tt  = dPPerdTheta(i,j,k) + 0.5*(1.-sigma(i,j,k))*dBsqdTheta(i,j,k)
                 vecx(i,j,k) = jacobian(i,j,k)/f(j)**2*c0*(tz*xpz+tt*xpt)
              ELSE ! Newton method
!                 vecx(i,j,k) = jacobian(i,j,k) / f(j)**2 * (-1./sigma(i,j,k) * dPperdAlpha(i,j,k) &
!                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j)**2 * fzet(k) * (gradRhoSq(i,j,k)* &
!                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradRhoGradZeta(i,j,k)) * &
!                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k))*0.5*dBsqdTheta(i,j,k)) - &
!                      (1. - sigma(i,j,k)) / sigma(i,j,k) * 0.5 * dBsqdAlpha(i,j,k))
!                 vecx(i,j,k) = vecx(i,j,k) - jacobian(i,j,k)/((f(j))**2) * dBBdAlpha(i,j,k) * alphaVal(k)
!                 vecd(i,j,k) = vecd(i,j,k) + jacobian(i,j,k)/(f(j)**2) * dBBdAlpha(i,j,k)
              END IF
           END DO
        END DO
     END DO
  END IF

  RETURN
END SUBROUTINE newk

!------------------------------------------------------------------------------
SUBROUTINE newj
  !     define right hand side of J dot grad(alfa) equilibrium equation

  USE ModScbMain,      ONLY: DP
  use ModScbParams,    ONLY: iOuterMethod, isotropy
  USE ModScbGrids,     ONLY: nthe, npsi, nzeta
  use ModScbVariables, ONLY: f, fzet, jacobian, bsq, psiVal, GradRhoGradTheta, &
                             GradRhoGradZeta, GradThetaGradZeta, GradRhoSq, &
                             dPdPsi, dSqPdPsiSq, dPPerdPsi, dBsqdPsi, dBsqdTheta, &
                             dPPerdTheta, dBBdPsi, sigma, dPPerdZeta, dBsqdZeta, &
                             dPPerdRho, dBsqdRho, GradZetaSq, vecr, vecd, &
                             nThetaEquator, nZetaMidnight, left, right
  

  implicit none

  INTEGER         :: i, j, k
  REAL(DP) :: xpr, xpt, c0, tr, tt

  Isotropic_case: IF (isotropy == 1) THEN
     DO k = 1, nzeta
        Newton_method:   IF (iOuterMethod == 2) THEN
           DO j = left, right
              vecr(1:nthe, j, k) = jacobian(1:nthe, j, k)/fzet(k)**2 * (dpdPsi(1:nthe, j, k) - &
                                    dSqPdPsiSq(1:nthe,j,k)*psival(j))
              vecd(1:nthe, j, k) = vecd(1:nthe, j, k) + jacobian(1:nthe,j,k)/fzet(k)**2 *dSqPdPsiSq(1:nthe,j,k)
           END DO
        ELSE            ! Picard method
           vecr(1:nthe, 1:npsi, k) = jacobian(1:nthe, 1:npsi, k) * dpdPsi(1:nthe, 1:npsi, k) / fzet(k)**2
        END IF Newton_method
     END DO
  ELSE       ! Anisotropic case
     DO k = 1, nzeta
        DO j = left, right
           DO i = 1, nthe
              IF (iOuterMethod == 1) THEN ! Picard method
                 xpr = (GradRhoGradZeta(i,j,k)**2 - GradRhoSq(i,j,k)*GradZetaSq(i,j,k))
                 xpt = (GradRhoGradZeta(i,j,k)*GradThetaGradZeta(i,j,k) &
                       - GradZetaSq(i,j,k)*GradRhoGradTheta(i,j,k))
                 c0  = -(f(j)*fzet(k)**2)/sigma(i,j,k)/bsq(i,j,k)
                 tr  = dPPerdRho(i,j,k) + 0.5*(1.-sigma(i,j,k))*dBsqdRho(i,j,k)
                 tt  = dPPerdTheta(i,j,k) + 0.5*(1.-sigma(i,j,k))*dBsqdTheta(i,j,k)
                 vecr(i,j,k) = jacobian(i,j,k)/fzet(k)**2*c0*(tr*xpr +tt*xpt)
              ELSE ! Newton method
!                 vecr(i,j,k) = jacobian(i,j,k) / fzet(k)**2 * (1./sigma(i,j,k) * dPperdPsi(i,j,k) &
!                      - 1./(sigma(i,j,k)*bsq(i,j,k)) * f(j) * fzet(k)**2 * (gradRhoGradZeta(i,j,k)* &
!                      gradThetaGradZeta(i,j,k) - gradRhoGradTheta(i,j,k)*gradZetaSq(i,j,k)) * &
!                      (dPperdTheta(i,j,k) + (1.-sigma(i,j,k)) * 0.5_dp * dBsqdTheta(i,j,k)) + &
!                      (1.-sigma(i,j,k)) / sigma(i,j,k) * 0.5_dp * dBsqdPsi(i,j,k))
!                 vecr(i,j,k) = vecr(i,j,k) - jacobian(i,j,k)/(fzet(k)**2) * dBBdPsi(i,j,k)*psival(j)
!                 vecd(i,j,k) = vecd(i,j,k) + jacobian(i,j,k)/(fzet(k)**2) * dBBdPsi(i,j,k)
              END IF
           END DO
        END DO
     END DO

  END IF Isotropic_case

  RETURN

END SUBROUTINE newj

END MODULE ModScbEquation
