      PROGRAM TOV
c
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
c --- Desctiption for parameters.
c
c     qini  : initial for emden function. 
c     pinx  : polytropic index
c     nstep : see Runge Kutta subroutine.
c     ndiv  : see Runge Kutta subroutine.
c     radini: initial radius for hunting radius from inside.
c     itype : itype = 0  iteration for compactness
c             itype = 1  single structure
c     chope : compactness to find (for itype = 0)
c
c --- Description for variables.
c
c     XE    : Schwartzscild radial coordinate 
c     YN(1) : Effective gravitational mass
c     YN(2) : Emden function
c     YN(3) : Proper rest mass 
c     YN(4) : Proper mass energy
c     YN(5) : Conformal factor (proper boundary condition not applied)
c     YN(6) : ADM mass energy  (proper boundary condition not applied)
c     adm   : ADM mass energy
c     compa : Compaction = YN(1)/XE = 
c               = Effective gravitational mass/
c                 Radius of star in Schwartzscild radial coordinate 
c
c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(NEQ = 6)
       EXTERNAL TEST
       common / misc / pinx, pi
       DIMENSION Y0(NEQ), YN(NEQ), WORK(NEQ,2)
c
      pi = 3.14159265358979d+0
      open (2,file='ovpar.dat',status='old')
      open (1,file='ovphy.dat',status='unknown')
      open (8,file='ovlas.dat',status='unknown')
      open (9,file='ovaux.dat',status='unknown')
      read(2,41) nstep, ndiv, radiini
      read(2,41) itype, idum, chope
 40   format(1p,2e15.7)
 41   format(2i5, 1p,1e15.7)
      close(2)
c
c --  read_parameter in grid_parameter.f90
c --  
      open (2,file='rnspar.dat',status='old')
      read(2,'(4i5)') ndum
      read(2,'(4i5)') ndum
      read(2,'(2i5)') ndum
      read(2,'(1p,3e10.3)') dum
      read(2,'(/,1i5)') ndum
      read(2,'(1p,2e10.3)') dum
      read(2,'(1p,2e10.3)') dum
      read(2,'(2(3x,a2),3x,3a1)') ndum
      read(2,'(1p,2e10.3)') dum
      read(2,'(/,2i5)') ndum
      read(2,'(1p,2e14.6)') emdc_ini, pinx
      read(2,'(1p,2e14.6)') restmass_sph, gravmass_sph
      read(2,'(1p,2e14.6)') rMoverR_sph, schwarz_radi_sph
c
c --  set compactness and initial
c
      chope = rMoverR_sph
      qini = emdc_ini
c
      ermx = 1.0d-10
      itmx = 100
c
c --- Iteration for compactness.
c
      itcom = 0
      iqc = 0
      dqdq = qini*0.1d0
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
      pini = qini**(pinx+1.0d0)
c
      X0 = 0.0D+0
      Y0(1) = 0.0D+0
      Y0(2) = qini
      Y0(3) = 0.0D+0
      Y0(4) = 0.0D+0
      Y0(5) = 2.0D+0
      Y0(6) = 0.0D+0
      YN(1) = 0.0D+0
      YN(2) = qini
      YN(3) = 0.0D+0
      YN(4) = 0.0D+0
      YN(5) = 2.0D+0
      YN(6) = 0.0D+0
      XN = radi
      M = nstep
      H = (XN - X0) / M
      N = ndiv
c
       if (error.le.ermx.and.itype.eq.1)
     & WRITE( 8 ,600) XE,YN(2),YN(2)**pinx,YN(2)**(1.0+pinx),YN(5)
c
      DO 10 I = 1,M
       ii = i
       XE = X0 + H
       CALL RK(NEQ,TEST,X0,XE,N,Y0,YN,WORK)
       if (error.le.ermx) then
       if (itype.eq.0.and.ii.eq.M) go to 2000
       if (itype.eq.1) then
       WRITE( 8 ,600) XE,YN(2),YN(2)**pinx,YN(2)**(1.0+pinx),YN(5)
       if (ii.eq.M) then
       YR = XE/YN(1)*(1.0d0-(1.0d0-2.0d0*YN(1)/XE)**0.5)
       adm =  YN(5)*YN(6)/YR
       compa = YN(1)/XE
       WRITE(6 ,630) XE,YN(1),YN(2),YN(3),YN(4), 
     &                  YN(5),YN(6),qini, compa, adm
       WRITE(1 ,630) XE,YN(1),YN(2),YN(3),YN(4),
     &                  YN(5),YN(6),qini, compa, adm
c
       stop ' end execution '
       end if
       end if
       end if
       if (YN(2).le.0.0d0) go to 1000
  10  CONTINUE
c
c ---- HUNTING for precise radius.
c
 1000 continue
c
       if (itrad.eq.1.and.ii.le.M.and.YN(2).le.0.0d0) then
       write(6,*) ' bad initial radius '
       radi = 0.95d0*radi
       go to 30
       end if
       if (YN(2).le.0.0d0) then
       radi = radib
       error = dr/radi
       dr = 0.2d0*dr
       WRITE(9,*)' back ',itrad, ii, radi, error
       end if
       if (ii.eq.M.and.YN(2).gt.0.0d0) then
       radib = radi
       radi = radi + dr
       WRITE(9,*)' hunt ',itrad, ii, radi, error
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
c ---- HUNTING END.
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
      qinib= qini
      qini = qini + dqdq
      dqdq = dabs(dqdq)
      erer = 1.0d0
      write(6,620) itcom, erer, qinib, compa
      go to 2222
      end if
c
      if (chope.ge.compa.and.compa.ge.compab) then
      compab = compa
      qinib= qini
      qini = qini + dqdq
      else if (chope.le.compa.and.compa.le.compab) then
      compab = compa
      qinib= qini
      qini = qini - dqdq
      else
      compab = compa
      qinibx = qini
      qini = 0.5d0*(qini+qinib)
      qinib = qinibx
      dqdq = 0.5d0*dabs(dqdq)
      end if
c
      erer = dabs((compa-chope)/chope)
ccc      write(6,*)compa, chope, erer
      write(6,620) itcom, erer, qinib, compa
      if (erer.lt.ermx) then
      WRITE(6 ,630) XE,YN(1),YN(2),YN(3),YN(4), 
     &                 YN(5),YN(6),qinib, compa, adm
      WRITE(1 ,630) XE,YN(1),YN(2),YN(3),YN(4),
     &                 YN(5),YN(6),qinib, compa, adm
c
      open(14,file='rnspar_add.dat',status='unknown') 
c
      emdc_ini = qini
      restmass_sph = YN(3)
      gravmass_sph = YN(1)
      rMoverR_sph   = compa
      schwarz_radi_sph = XE
      write(14,'(1p,2e14.6,a19)') emdc_ini, pinx, 
     & '   : emdc_ini, pinx'
      write(14,'(1p,2e14.6,a31)') restmass_sph, gravmass_sph, 
     & '   : restmass_sph, gravmass_sph'
      write(14,'(1p,2e14.6,a47)') rMoverR_sph, schwarz_radi_sph, 
     & '   : MoverR_sph,   schwarz_radi_sph  (K=1 unit)'
      write(6,'(/,/,1p,2e14.6,a19)') emdc_ini, pinx, 
     & '   : emdc_ini, pinx'
      write(6,'(1p,2e14.6,a31)') restmass_sph, gravmass_sph, 
     & '   : restmass_sph, gravmass_sph'
      write(6,'(1p,2e14.6,a47)') rMoverR_sph, schwarz_radi_sph, 
     & '   : MoverR_sph,   schwarz_radi_sph  (K=1 unit)'
c
      close(14)
c
      stop ' end of execution '
      end if
c
      go to 2222
c
 620   FORMAT(1i5,1p,e12.4,3e16.8)
 630   FORMAT(1p,5e14.6)
      STOP 
      END
C
c --------------------------------------------------------
c --- equations are in this routine.
c
      SUBROUTINE TEST(X,Y,F)
       IMPLICIT REAL*8(A-H,O-Z)
       PARAMETER(NEQ = 6)
       common / misc / pinx, pi
       DIMENSION Y(NEQ), F(NEQ)
c
      rr  = X
      rma = Y(1)
      emd = Y(2)
      bma = Y(3)
      tma = Y(4)
      psi = Y(5)
      adm = Y(6)
      if (emd.le.0.0d0) then
      rho0= 0.0d0
      rho = 0.0d0
      pre = 0.0d0
      hhh = 1.0d0
      end if
      if (emd.gt.0.0d0) then
      rho0= emd**pinx
      rho = emd**pinx*(1.0d0 + pinx*emd)
      pre = emd**(1.0d0 + pinx)
      hhh = 1.0d0 + (1.0d0+pinx)*emd
      end if
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
      F(1) = 4.0d0*pi*rr**2*rho
      F(2) = - 1.0d0/(1.0d0+pinx)*hhh*(rma + 4.0d0*pi*pre*rr**3)/
     &               (rr**2 - 2.0d0*rma*rr)
      F(3) =  4.0d0*pi*rr**2*rho0/(1.0d0-2.0d0*rma/rr)**0.5d0
      F(4) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0
      F(5) =  0.5d0*psi/rr*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rr)**0.5d0)
      F(6) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0/psi
      end if
c 
      RETURN
      END
c
c ---------------------------------------------------------
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
