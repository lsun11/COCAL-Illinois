subroutine update_parameter_axisym_WL_peos_drot(convf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, nrf_deform, nrf, ntgxy
  use def_metric, only : psi, alph, bvyd, tfkij, tfkijkij
  use def_metric_hij, only : hyyd
  use def_matter, only : omef
  use def_matter_parameter
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phig
  use interface_minv
  implicit none
  real(long), intent(in) :: convf
  real(long)     ::    sgg(1:4), gg(1:4,1:4)
! "bin" is beta(in) and so on.
  real(long)     ::    bin, bmid, bout
! ain is log(\alpha^((1/R_i)^2)) (in). R_in is in fact the previous R.
  real(long)     ::    ain, amid, aout
! pain = (\phi^2/\alpha)^((1/R_i)^2) (in).
  real(long)     ::    pain, pamid, paout
! ovin is \omega(in). \omega = \beta + \Omega \phi.
  real(long)     ::    ovin, ovmid, ovout
! "ov2in" is f_{ab} \omega^a \omega^b.
  real(long)     ::    ov2in, ov2mid, ov2out
! "jin" is Integral of J(Omega).
  real(long)     ::    jin, jmid, jout, jint_in, jint_mid, jint_out
  real(long)     ::    alin, almid, alout
  real(long)     ::    psin, psmid, psout
  real(long)     ::    hyin, hymid, hyout
  real(long)     ::    gin, gmid, gout
  real(long)     ::    omein, omemid, omeout
!
  real(long)     ::    termin, termmid, termout
  real(long)     ::    dtermdoin, dtermdomid, dtermdoout
  real(long)     ::    dtermdrin, dtermdrmid, dtermdrout      
  real(long)     ::    djintdoin, djintdomid, djintdoout
  real(long)     ::    dodoc_in, dodoc_mid, dodoc_out
  real(long)     ::    vphiyin, vphiymid, vphiyout
  real(long)     ::    loghc
  real(long)     ::    radiold, berold, omeold
  real(long)     ::    ddradi, ddber, ddome
  real(long)     ::    facfac, numero, error, hc, pre, rho, ene
  real(long)     ::    ome_mx, ome_eq, Jome, Jome_int, dum, small = 1.0d-30
  integer        ::    nin, nmx
!##############################################
  real(long)     ::    ut, hh, emdtest
!##############################################
  nin = nrf_deform
!
  sgg = 0.0d0
  gg  = 0.0d0
!
  call peos_q2hprho(emdc, hc, pre, rho, ene)
  call search_emdmax_xaxis_grid(nmx)
!  nmx = 0
!testtest
  loghc = dlog(hc)
  vphiyin  = vec_phig(nin,0,0,2)
  vphiymid = vec_phig(nmx,ntgxy,0,2)
  vphiyout = vec_phig(nrf,ntgxy,0,2)
  omeold = ome
  berold = ber
  radiold = radi
  numero = 1
! ################################################################ 
  ain  = dlog(alph(nin,0,0)**(1.d0/radi**2))
  pain = (psi(nin,0,0)**2/alph(nin,0,0))**(1.d0/radi**2)
  bin  = bvyd(nin,0,0)
  hyin = hyyd(nin,0,0)
  gin  = 1.0d0 + hyin
  psin = psi(nin,0,0)
  alin = alph(nin,0,0)
  omein = omef(nrf,0,0)
!
  amid = dlog(alph(nmx,ntgxy,0)**(1.d0/radi**2))
  pamid= (psi(nmx,ntgxy,0)**2/alph(nmx,ntgxy,0))**(1.d0/radi**2)
  bmid = bvyd(nmx,ntgxy,0)
  hymid= hyyd(nmx,ntgxy,0)
  gmid = 1.0d0 + hymid
  psmid= psi(nmx,ntgxy,0)
  almid= alph(nmx,ntgxy,0)
  omemid=omef(nmx,ntgxy,0)
!
  aout = dlog(alph(nrf,ntgxy,0)**(1.d0/radi**2))
  paout= (psi(nrf,ntgxy,0)**2/alph(nrf,ntgxy,0))**(1.d0/radi**2)
  bout = bvyd(nrf,ntgxy,0)
  hyout= hyyd(nrf,ntgxy,0)
  gout = 1.0d0 + hyout
  psout= psi(nrf,ntgxy,0)
  alout= alph(nrf,ntgxy,0)
  omeout=omef(nrf,ntgxy,0)
!
  ome_mx = omemid
  ome_eq = omeout
!
!  write (6,*) "test values"
!  write (6,'(1p,3e12.4)') ome, ber, radi
! ! write (6,'(1p,3e12.4)') psi(nin,0,0),psi(0,0,0),psi(nrf,ntgxy,0)
! ! write (6,'(1p,3e12.4)') alph(nin,0,0),alph(0,0,0),alph(nrf,ntgxy,0)
! ! write (6,'(1p,3e12.4)') error
!  write (6,'(1p,3e12.4)') vphiyin, vphiymid, vphiyout
!  write (6,'(1p,3e12.4)') alin, almid, alout
!  write (6,'(1p,3e12.4)') psin, psmid, psout
!  write (6,'(1p,3e12.4)') bin, bmid, bout
!  write (6,'(1p,3e12.4)') omein, omemid, omeout
!  write (6,'(1p,3e12.4)') jin, jmid, jout
!  write (6,'(1p,3e12.4)') jint_in, jint_mid, jint_out
!  write (6,'(1p,3e12.4)') ome_mx, ome_eq
!! stop
! ################################################################ 
  do
    call calc_omega_drot(vphiyin,alin,psin,bin,hyin,omein,Jome,Jome_int)
    jin  = Jome ; jint_in  = Jome_int
    call calc_omega_drot(vphiymid,almid,psmid,bmid,hymid,omemid,Jome,Jome_int)
    jmid = Jome ; jint_mid = Jome_int
    call calc_omega_drot(vphiyout,alout,psout,bout,hyout,omeout,Jome,Jome_int)
    jout = Jome ; jint_out = Jome_int
!
    ovin  = bin  + omein *vphiyin
    ovmid = bmid + omemid*vphiymid
    ovout = bout + omeout*vphiyout
! This is in fact only the \omega^y component. 
! Since only \phi^y is not 0. (We are at x axis.) 
! This is also the reason why I have those bvxd down there.
!      ov2in  = ovin**2 +bvxd(nin,0,0)**2+bvzd(nin,0,0)**2
!      ov2mid = ovmid**2+bvxd(0,0,0)**2  +bvzd(0,0,0)**2
!      ov2out = ovout**2+bvxd(nrf,ntgxy,0)**2+bvzd(nrf,ntgxy,0)**2
    ov2in  = gin *ovin**2
    ov2mid = gmid*ovmid**2
    ov2out = gout*ovout**2
    termin  = 1.0d0 - pain**( 2.d0*radi**2)*ov2in
    termmid = 1.0d0 - pamid**(2.d0*radi**2)*ov2mid
    termout = 1.0d0 - paout**(2.d0*radi**2)*ov2out
! termin is 1-(\phi^4/\alpha^2)^((R_0/R_i)^2) \omega^2            
!
    dodoc_in  = A2DR*index_DR*termin & 
    &         /(pain **(2.d0*radi**2)*gin *vphiyin **2 &
    &         + A2DR*index_DR*termin &
    &         + 2.0d0*pain **(2.d0*radi**2)*gin*vphiyin*ome**2*ovin*jin&
    &         - ome*termin*jin )
    dodoc_mid = A2DR*index_DR*termmid*(ome/(omemid+small))**(index_DR-1.0d0) &
    &         /(pamid**(2.d0*radi**2)*gmid*vphiymid**2 &
    &         + A2DR*index_DR*termmid*(ome/(omemid+small))**index_DR &
    &         + 2.0d0*pamid**(2.d0*radi**2)*gmid*vphiymid*omemid**2*ovmid*jmid&
    &         - omemid*termmid*jmid)
    dodoc_out = A2DR*index_DR*termout*(ome/(omeout+small))**(index_DR-1.0d0) &
    &         /(paout**(2.d0*radi**2)*gout*vphiyout**2 &
    &         + A2DR*index_DR*termout*(ome/(omeout+small))**index_DR &
    &         + 2.0d0*paout**(2.d0*radi**2)*gout*vphiyout*omeout**2*ovout*jout&
    &         - omeout*termout*jout)
!
  if (index_DR.ne.2.0d0) then
    djintdoin  = 0.0d0  !+ jin *dodoc_in
    djintdomid = A2DR*index_DR*omemid*(ome/(omemid+small))**(index_DR-1.0d0) &
    &                             /(2.0d0-index_DR) &
    &          - A2DR*index_DR*ome/(2.0d0-index_DR) + jmid*dodoc_mid
    djintdoout = A2DR*index_DR*omeout*(ome/(omeout+small))**(index_DR-1.0d0) &
    &                             /(2.0d0-index_DR) &
    &          - A2DR*index_DR*ome/(2.0d0-index_DR) + jout*dodoc_out
  else
    djintdoin  = 0.0d0 ! + jin *dodoc_in
    djintdomid = - 2.0d0*A2DR*ome*dlog(ome/(omemid+small)) + jmid*dodoc_mid
    djintdoout = - 2.0d0*A2DR*ome*dlog(ome/(omeout+small)) + jout*dodoc_out
  end if
!
    dtermdoin = - pain**( 2.d0*radi**2)*2.0d0*gin *ovin *vphiyin*dodoc_in
    dtermdomid= - pamid**(2.d0*radi**2)*2.0d0*gmid*ovmid*vphiymid*dodoc_mid
    dtermdoout= - paout**(2.d0*radi**2)*2.0d0*gout*ovout*vphiyout*dodoc_out
    dtermdrin =-dlog(pain )*4.d0*radi*pain**( 2.d0*radi**2)*ov2in
    dtermdrmid=-dlog(pamid)*4.d0*radi*pamid**(2.d0*radi**2)*ov2mid
    dtermdrout=-dlog(paout)*4.d0*radi*paout**(2.d0*radi**2)*ov2out
!
    sgg(1) = -(radi**2*ain  + 0.5d0*dlog(termin)  + jint_in  - dlog(ber))
    sgg(2) = -(radi**2*aout + 0.5d0*dlog(termout) + jint_out - dlog(ber))
    sgg(3) = -(radi**2*amid + 0.5d0*dlog(termmid) + jint_mid - dlog(ber) &
             + loghc)
!
!  ain is ln(\alpha)
    gg(1,1) =   0.5d0*dtermdoin/termin + djintdoin
    gg(1,2) = - 1.0d0/ber
    gg(1,3) = 2.0d0*radi*ain + 0.5d0*dtermdrin/termin
!
    gg(2,1) =   0.5d0*dtermdoout/termout + djintdoout
    gg(2,2) = - 1.0d0/ber
    gg(2,3) = 2.0d0*radi*aout + 0.5d0*dtermdrout/termout
!
    gg(3,1) =   0.5d0*dtermdomid/termmid + djintdomid
    gg(3,2) = - 1.0d0/ber
    gg(3,3) = 2.0d0*radi*amid + 0.5d0*dtermdrmid/termmid
!
!    test
!if (numero.eq.1) then
! write (6,'(1p,3e12.4)') sgg(1), sgg(2), sgg(3)
! write (6,'(1p,3e12.4)') gg(1,1), gg(1,2), gg(1,3)
! write (6,'(1p,3e12.4)') gg(2,1), gg(2,2), gg(2,3)
! write (6,'(1p,3e12.4)') gg(3,1), gg(3,2), gg(3,3)
!end if
!stop
!!pause
!    test
    call minv(gg,sgg,3,4)
    ddome = sgg(1)
    ddber = sgg(2)
    ddradi = sgg(3)
    facfac = dmin1(dble(numero)/5.0d0,1.0d0)
    ome = ome + ddome*facfac
    ber = ber + ddber*facfac
    radi = radi + ddradi*facfac
    error = dmax1(dabs(ddome/ome),dabs(ddber/ber),dabs(ddradi/radi))
    numero = numero + 1
    if (numero>1000) then
      write(6,*)' numero = ', numero, '   error =',error
    end if
!      if (iter<10.and.numero>10) exit
    if (numero>1010) exit
!###    if (dabs(error)<1.0d-08) exit
    if (dabs(error)<1.0d-12) exit
  end do
!    test
! write (6,*) "test values"
! write (6,*) sgg(1), sgg(2), sgg(3)
! write (6,*) 'update paramter', ome
!    test
!      write(6,*)' numero = ', numero, '  error = ',error
!
  if (radi < 0.d0) write(6,*) ' ### radi minus ###'
  if (radi < 0.d0) radi = - radi
  if (ome < 0.d0) write(6,*) ' ### ome minus ###'
  if (ome < 0.d0) ome = - ome
  if (ber < 0.d0) write(6,*) ' ### ber minus ###'
  if (ber < 0.d0) ber = - ber
  ome = convf*ome + (1.d0-convf)*omeold
  ber = convf*ber + (1.d0-convf)*berold
  radi = convf*radi + (1.d0-convf)*radiold
!
!
! --  Improving alpha and psi.  
!
  alph = alph**((radi/radiold)**2)
  psi = psi**((radi/radiold)**2)
!##############################################
!  ut = 1/sqrt(alph(nrf,ntgxy,0)**2 & 
!     -psi(nrf,ntgxy,0)**4*ov2out)
!  hh = ber*ut
!  emdtest = 1.0d0/(pinx+1.0d0)*(hh-1.0d0)
!!write (6,*) "h  = ", hh
!
end subroutine update_parameter_axisym_WL_peos_drot
