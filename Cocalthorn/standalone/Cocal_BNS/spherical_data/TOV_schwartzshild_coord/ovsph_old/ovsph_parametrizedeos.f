      PROGRAM TOV
c
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
c     G = c = Msol = 1 unit is used in TOV integration.
c     assume pre = 10^37 dyn/cm^2 at rho = 10^16 gr/cm^3 to specify K.
c
c --- Desctiption for input parameters.
c
c  -- Input in 'ovpar_parameos.dat'
c
c     hini  : not used.
c     nstep : see Runge Kutta subroutine.
c     ndiv  : see Runge Kutta subroutine.
c     radini: initial radius for the surface of a neutron star.
c            (Iteration is made until "radi" to become a surface exactly.)
c     itype : itype = 0  iteration for compactness
c             itype = 1  single structure
c     chope : compactness to find (for itype = 0)
c
c
c  -- Input in "oveosparm.dat
c
c     nphase : total number of piecewise regions.
c     rhoini : central value of the density [g/cm^3]. 
c             (initial for the integration.)
c     rhocgs : value of the density in a piecewise region.  
c     abi    : value of the adiabatic index (Gamma) in a piecewise region.  
c
c --- Description for the variables.
c
c  -- Output in "ovphy.dat" : Mass, radius and central enthalpy
c
c       WRITE(1,630) XE,YN(1),YN(2),YN(3),YN(4),
c     &                 YN(5),YN(6),hini, compa, adm
c
c     XE    : Schwartzscild radial coordinate at the surface
c     YN(1) : Gravitational mass
c     YN(2) : Pressure
c     YN(3) : Proper rest mass 
c     YN(4) : Proper mass energy
c     YN(5) : Conformal factor (proper boundary condition is NOT applied)
c     YN(6) : Qantity relates to ADM mass energy.  
c     hini  : the relativistic enthalpy at the center.
c     adm   : ADM mass energy
c     compa : Compactness = YN(1)/XE = adm/XE
c               = Gravitational mass/
c                 Radius of star in Schwartzscild radial coordinate 
c
c  -- Output in "ovlas.dat" : profile of thermodynamic variables.
c
c     WRITE(8,600) XE,hh,rho0,pre,ene,YN(5)
c
c     XE    : Schwartzscild radial coordinate 
c     hh    : the relativistic enthalpy.
c     rho0  : the rest mass density.
c     pre   : the pressure.
c     ene   : the energy density.
c     YN(5) : Conformal factor (proper boundary condition is NOT applied)
c
c  -- Ouput in "oveos_phase.dat" : data at the interface of each piecewise 
c                                  region, rescaled in G=c=Msol=1 unit.
c
c      write(30,*)'#', nphase, qini
c      write(30,40) ii, abc(ii), abi(ii), rhoi(ii), qi(ii), hi(ii), 
c     &              rhocgs(ii)
c
c      qini   : p/rho at the center rescaled in G=c=Msol=1 unit.
c      abc    : adiabatic constant K_i rescaled in G=c=Msol=1 unit.
c      abi    : adiabatic index Gamma_i.
c      rhoi   : rest mass density rescaled in G=c=Msol=1 unit.
c      qi     : p/rho rescaled in G=c=Msol=1 unit.
c      hi     : the relativistic enthalpy rescaled in G=c=Msol=1 unit.
c      rhocgs : rest mass density in [g/cm^3].
c
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NEQ = 6)
      EXTERNAL TOVEQS
      common / misc / pinx, pi
      DIMENSION Y0(NEQ), YN(NEQ), WORK(NEQ,2)
      common / eospar / abc(0:10), abi(0:10), 
     &                 rhoi(0:10), qi(0:10), hi(0:10), nphase
c
      pi = 3.14159265358979d+0
      open (2,file='ovpar_parameos.dat',status='old')
      open (1,file='ovphy.dat',status='unknown')
      open (8,file='ovlas.dat',status='unknown')
      open (9,file='ovaux.dat',status='unknown')
c
      read(2,41) nstep, ndiv, radiini
      read(2,41) itype, idum, chope
 40   format(1p,2e15.7)
 41   format(2i5, 1p,1e15.7)
      close(2)
c
      ermx = 1.0d-10
      itmx = 100
c
      call eosparam(qini)
      call eoslookup(qini,qi,abc,abi,nphase,abct,abin,iph)
      call q2hprho(qini,hi(iph-1),qi(iph-1),abin,abct,hini,pre,rho0)
c
      hh = hini
c
c --- Iteration for compactness.
c
      itcom = 0
cc      iqc = 0
      dqdq = hini*0.1d0
c
 2222 continue
      itcom = itcom + 1
c
c --- Iteration for radius.
c
      radi = radiini
      itrad = 0
      dr = 0.2d0*radi
      error = dr/radi
c
 200  continue
      itrad = itrad + 1
c
      X0 = 0.0D+0
      Y0(1) = 0.0D+0
      Y0(2) = hini
      Y0(3) = 0.0D+0
      Y0(4) = 0.0D+0
      Y0(5) = 2.0D+0
      Y0(6) = 0.0D+0
      YN(1) = 0.0D+0
      YN(2) = hini
      YN(3) = 0.0D+0
      YN(4) = 0.0D+0
      YN(5) = 2.0D+0
      YN(6) = 0.0D+0
      XN = radi
      M = nstep
      H = (XN - X0) / M
      N = ndiv
c
      if (error.le.ermx.and.itype.eq.1) then
      call eoslookup(hini,hi,abc,abi,nphase,abct,abin,iph)
      call h2qprho(hini,hi(iph-1),qi(iph-1),abin,abct,q,pre,rho0,ene)
      WRITE(8,600) X0,hini,rho0,pre,ene,YN(5)
      end if
c
      DO 10 I = 1,M
       ii = i
       XE = X0 + H
       CALL RK(NEQ,TOVEQS,X0,XE,N,Y0,YN,WORK)
c
      hh = YN(2)
      call eoslookup(hh,hi,abc,abi,nphase,abct,abin,iph)
      call h2qprho(hh,hi(iph-1),qi(iph-1),abin,abct,q,pre,rho0,ene)
c
       if (error.le.ermx) then
       if (itype.eq.0.and.ii.eq.M) go to 2000
       if (itype.eq.1) then
       WRITE(8,600) XE,hh,rho0,pre,ene,YN(5)
       if (ii.eq.M) then
       YR = XE/YN(1)*(1.0d0-(1.0d0-2.0d0*YN(1)/XE)**0.5)
       adm =  YN(5)*YN(6)/YR
       compa = YN(1)/XE
       WRITE(6,630) XE,YN(1),YN(2),YN(3),YN(4), 
     &                 YN(5),YN(6),hini, compa, adm
       WRITE(1,630) XE,YN(1),YN(2),YN(3),YN(4),
     &                 YN(5),YN(6),hini, compa, adm
       stop ' end execution '
       end if
       end if
       end if
c
       small  = 1.0d0
       hh = YN(2)
       if (hh.le.small) go to 1000
c
  10  CONTINUE
c
c ---- Look for a precise radius.
c
 1000 continue
c
c --  using hh
c
       if (itrad.eq.1.and.ii.le.M.and.hh.le.1.0d0) then
       write(6,*) ' bad initial radius '
       radi = 0.95d0*radi
       go to 30
       end if
       if (hh.le.1.0d0) then
       radi = radib
       error = dr/radi
       dr = 0.2d0*dr
       WRITE(9,*)' back ',itrad, ii, radi, error
       end if
       if (ii.eq.M.and.hh.gt.1.0d0) then
       radib = radi
       radi = radi + dr
       WRITE(9,*)' hunt ',itrad, ii, radi, error, hh
       end if
c
 30    continue
c
 610   FORMAT(2i5,1p,4e12.4)
 600   FORMAT(1P ,7(E15.7))
c
      if (itrad.le.itmx) go to 200
      if (itrad.gt.itmx) stop ' iter '
c
c ---- Look for the radius END.
c
c ---- Searching for a certain value of compactness.
c
 2000 continue
c
      YR = XE/YN(1)*(1.0d0-(1.0d0-2.0d0*YN(1)/XE)**0.5)
      adm =  YN(5)*YN(6)/YR
      compa = YN(1)/XE
c
      if (itcom.eq.1) then
      if (compa.gt.chope) dqdq = - dqdq
      compab = compa
      hinib= hini
      hini = hini + dqdq
      dqdq = dabs(dqdq)
      erer = 1.0d0
      write(6,620) itcom, erer, hinib, compa
      go to 2222
      end if
c
      if (chope.ge.compa.and.compa.ge.compab) then
      compab = compa
      hinib= hini
      hini = hini + dqdq
      else if (chope.le.compa.and.compa.le.compab) then
      compab = compa
      hinib= hini
      hini = hini - dqdq
      else
      compab = compa
      hinibx = hini
      hini = 0.5d0*(hini+hinib)
      hinib = hinibx
      dqdq = 0.5d0*dabs(dqdq)
      end if
c
      erer = dabs((compa-chope)/chope)
      write(6,620) itcom, erer, hinib, compa
      if (erer.lt.ermx) then
      WRITE(6 ,630) XE,YN(1),YN(2),YN(3),YN(4), 
     &                 YN(5),YN(6),hinib, compa, adm
      WRITE(1 ,630) XE,YN(1),YN(2),YN(3),YN(4),
     &                 YN(5),YN(6),hinib, compa, adm
      stop ' end of execution '
      end if
c
      go to 2222
c
 620   FORMAT(1i5,1p,e12.4,3e16.8)
 630   FORMAT(1p,5e14.5)
      STOP 
      END
C
c ______________________________________________________________________
c ______________________________________________________________________
      include 'GR_pmeos.f'
c ______________________________________________________________________
c ______________________________________________________________________
c
      SUBROUTINE TOVEQS(X,Y,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NEQ = 6)
      common / misc / pinx, pi
      common / eospar / abc(0:10), abi(0:10), 
     &                 rhoi(0:10), qi(0:10), hi(0:10), nphase
      DIMENSION Y(NEQ), F(NEQ)
c
c --- TOV equations are in this routine.
c
      rr  = X
      rma = Y(1)
      hh  = Y(2)
      bma = Y(3)
      tma = Y(4)
      psi = Y(5)
      adm = Y(6)
      if (hh.le.1.0d0) then
      pre = 0.0d0
      rho0= 0.0d0
      ene = 0.0d0
      end if
      if (hh.gt.1.0d0) then
      call eoslookup(hh,hi,abc,abi,nphase,abct,abin,iph)
      call h2qprho(hh,hi(iph-1),qi(iph-1),abin,abct,q,pre,rho0,ene)
      end if
c
c
      if (rr.eq.0.0d0) then
      F(1) = 0.0d0
      F(2) = 0.0d0
      F(3) = 0.0d0
      F(4) = 0.0d0
      F(5) = 0.0d0
      F(6) = 0.0d0
      end if
      if (rr.ne.0.0d0) then
      F(1) = 4.0d0*pi*rr**2*ene
cc      F(2) = - (hh + 0.93d0)*(rma + 4.0d0*pi*pre*rr**3)
cc      F(2) = - (hh + 1.0d0)*(rma + 4.0d0*pi*pre*rr**3)
      F(2) = - hh*(rma + 4.0d0*pi*pre*rr**3)
     &        /(rr**2 - 2.0d0*rma*rr)
      F(3) = 4.0d0*pi*rr**2*rho0/(1.0d0-2.0d0*rma/rr)**0.5d0
      F(4) = 4.0d0*pi*rr**2*ene/(1.0d0-2.0d0*rma/rr)**0.5d0
      F(5) = 0.5d0*psi/rr*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rr)**0.5d0)
      F(6) = 4.0d0*pi*rr**2*ene/(1.0d0-2.0d0*rma/rr)**0.5d0/psi
      end if
c 
      RETURN
      END
c
c ______________________________________________________________________
c ______________________________________________________________________
c
      SUBROUTINE RK(NEQ,FUNC,X0,XE,N,Y0,YN,WORK)
**********************************************************************
*     SUBROUTINE RK NUMERICALLY INTEGRATES A SYSTEM OF NEQ           *
*     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
*             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
*     BY THE CLASSICAL RUNGE-KUTTA FORMULA.                          *
*                                                                    *
*     PARAMETERS                                                     *
*  === INPUT ===                                                     *
*     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
*     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
*                F(I)=DY(I)/DX                                       *
*     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
*     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
*     (5) N: NUMBER OF DIVISIONS                                     *
*        THE INTERVAL (X0, XE) IS DIVIDED INTO N SUBINTERVALS        *
*        WITH THE LENGTH (XE-X0)/N AND IN EACH SUBINTERVAL           *
*        THE CLASSICAL RUNGE-KUTTA FORMULA IS USED.                  *
*     (6) Y0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                    *
*  === OUTPUT ===                                                    *
*     (7) YN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE             *
*  === OTHER ===                                                     *
*     (8) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,2)) TO BE         *
*                 USED INSIDE RK                                     *
*     COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989, V. 1                *
**********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION Y0(NEQ),YN(NEQ),WORK(NEQ,2)
      H = (XE - X0) / N
      DO 10 I = 1,N
       CALL RKSTEP(NEQ,FUNC,X0,H,Y0,YN,WORK(1,1),WORK(1,2))
       X0 = X0 + H
       DO 20 J = 1,NEQ
        Y0(J) = YN(J)
   20  CONTINUE
   10 CONTINUE
      X0 = XE
      RETURN
      END
C
      SUBROUTINE RKSTEP(NEQ,FUNC,X,H,Y0,YN,AK,W)
       IMPLICIT REAL * 8(A-H,O-Z)
       PARAMETER(A2 = 0.5D+0, A3 = A2)
       PARAMETER(B2 = 0.5D+0, B3 = B2)
       PARAMETER(C1 = 1.0D+0 / 6, C2 = 1.0D+0 / 3, C3 = C2, C4 = C1)
       DIMENSION Y0(NEQ),YN(NEQ),AK(NEQ),W(NEQ)
      CALL FUNC(X,Y0,AK)
      DO 10 I = 1,NEQ
       YN(I) = Y0(I) + H * C1 * AK(I)
   10 CONTINUE
      DO 20 I = 1,NEQ
       W(I) = Y0(I) + H * B2 * AK(I)
   20 CONTINUE
      CALL FUNC(X + A2 * H,W,AK)
      DO 30 I = 1,NEQ
       YN(I) = YN(I) + H * C2 * AK(I)
   30 CONTINUE
      DO 40 I = 1,NEQ
       W(I) = Y0(I) + H * B3 * AK(I)
   40 CONTINUE
      CALL FUNC(X + A3 * H,W,AK)
      DO 50 I = 1,NEQ
       YN(I) = YN(I) + H * C3 * AK(I)
   50 CONTINUE
      DO 60 I = 1,NEQ
       W(I) = Y0(I) + H * AK(I)
   60 CONTINUE
      CALL FUNC(X + H,W,AK)
      DO 70 I = 1,NEQ
       YN(I) = YN(I) + H * C4 * AK(I)
   70 CONTINUE
      RETURN
      END
c

