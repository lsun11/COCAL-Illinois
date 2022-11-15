subroutine update_parameter_spherical_peos(convf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, nrf_deform, nrf, ntgxy
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd, tfkij, tfkijkij
  use def_matter_parameter
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only   : sinphig, cosphig
  use def_vector_phi, only : vec_phig
  use interface_minv
  implicit none
  real(long), intent(in) :: convf
  real(long)     ::    sgg(1:4), gg(1:4,1:4)
  real(long)     ::    bin, bmid, bout
! "bin" is beta(in) and so on.
  real(long)     ::    ain, amid, aout
! ain is log(\alpha^((1/R_i)^2)) (in). R_in is in fact the previous R.
  real(long)     ::    pain, pamid, paout
! pain = (\phi^2/\alpha)^((1/R_i)^2) (in).
  real(long)     ::    ovin, ovmid, ovout
! ovin is \omega(in). \omega = \beta + \Omega \phi.
  real(long)     ::    ov2in, ov2mid, ov2out
! "ov2in" is f_{ab} \omega^a \omega^b.
  real(long)     ::    termin, termmid, termout
  real(long)     ::    dtermdoin, dtermdomid, dtermdoout
  real(long)     ::    dtermdrin, dtermdrmid, dtermdrout      
  real(long)     ::    vphiyin, vphiymid, vphiyout
  real(long)     ::    loghc
  real(long)     ::    radiold, berold, omeold
  real(long)     ::    ddradi, ddber, ddome
  real(long)     ::    facfac, numero, error, hc, pre, rho, ene
  integer        ::    nin
!##############################################
  real(long)     ::    ut, hh, emdtest
!##############################################
  nin = nrf_deform
!
  sgg = 0.0d0
  gg  = 0.0d0
!
  call peos_q2hprho(emdc, hc, pre, rho, ene)
  loghc = dlog(hc)
  vphiyin  = vec_phig(nin,0,0,2)
  vphiymid = vec_phig(0,0,0,2)
!out  vphiyout = vec_phig(nrf,ntgxy,0,2)
  omeold = ome
  berold = ber
  radiold = radi
  numero = 1
! ################################################################ 
  ain  = dlog(alph(nin,0,0)**(1.d0/radi**2))
  pain = (psi(nin,0,0)**2/alph(nin,0,0))**(1.d0/radi**2)
  bin  = bvyd(nin,0,0)
!
  amid = dlog(alph(0,0,0)**(1.d0/radi**2))
  pamid= (psi(0,0,0)**2/alph(0,0,0))**(1.d0/radi**2)
  bmid = bvyd(0,0,0)
!
!out  aout = dlog(alph(nrf,ntgxy,0)**(1.d0/radi**2))
!out  paout= (psi(nrf,ntgxy,0)**2/alph(nrf,ntgxy,0))**(1.d0/radi**2)
!out  bout = bvyd(nrf,ntgxy,0)
! write (6,*) "test values"
! write (6,'(1p,3e12.4)') ome, ber, radi
! write (6,'(1p,3e12.4)') psi(nin,0,0),psi(0,0,0)
! write (6,'(1p,3e12.4)') alph(nin,0,0),alph(0,0,0)
! write (6,'(1p,3e12.4)') error
! write (6,'(1p,3e12.4)') pain, pamid
!! ################################################################ 
 ome = 0.0d0
  do
    ovin  = bin  + ome*vphiyin
    ovmid = bmid + ome*vphiymid
!out    ovout = bout + ome*vphiyout
! This is in fact only the \omega^y component. 
! Since only \phi^y is not 0. (We are at x axis.) 
! This is also the reason why I have those bvxd down there.
!      ov2in  = ovin**2 +bvxd(nin,0,0)**2+bvzd(nin,0,0)**2
!      ov2mid = ovmid**2+bvxd(0,0,0)**2  +bvzd(0,0,0)**2
!      ov2out = ovout**2+bvxd(nrf,ntgxy,0)**2+bvzd(nrf,ntgxy,0)**2
    ov2in  = ovin**2 
    ov2mid = ovmid**2
!out    ov2out = ovout**2
    termin  = 1.0d0 - pain**( 2.d0*radi**2)*ov2in
    termmid = 1.0d0 - pamid**(2.d0*radi**2)*ov2mid
!out write(6,*)termin, termmid
!outstop
!out    termout = 1.0d0 - paout**(2.d0*radi**2)*ov2out
! termin is 1-(\phi^4/\alpha^2)^((R_0/R_i)^2) \omega^2            
!out    dtermdoin = - pain**(2.d0*radi**2)*2.0d0*ovin*vphiyin
!out    dtermdomid= - pamid**(2.d0*radi**2)*2.0d0*ovmid*vphiymid
!out    dtermdoout= - paout**(2.d0*radi**2)*2.0d0*ovout*vphiyout
    dtermdrin =-dlog(pain)*4.d0*radi*pain**(2.d0*radi**2)*ov2in
    dtermdrmid=-dlog(pamid)*4.d0*radi*pamid**(2.d0*radi**2)*ov2mid
!outwrite(6,*)    dtermdrin ,     dtermdrmid
!outstop
!out    dtermdrout=-dlog(paout)*4.d0*radi*paout**(2.d0*radi**2)*ov2out
!
!out    sgg(1) = -(radi**2*ain  + 0.5d0*dlog(termin)  - dlog(ber))
!out    sgg(2) = -(radi**2*aout + 0.5d0*dlog(termout) - dlog(ber))
!out    sgg(3) = -(radi**2*amid + 0.5d0*dlog(termmid) - dlog(ber) &
!out             + loghc)
!    sgg(1) = -(radi**2*ain  + 0.5d0*dlog(termin)  - dlog(ber))
!    sgg(2) = -(radi**2*amid + 0.5d0*dlog(termmid) - dlog(ber) &
!    &        + loghc)
    sgg(1) = -(radi**2*ain  - dlog(ber))
    sgg(2) = -(radi**2*amid - dlog(ber) + loghc)
!
!  ain is ln(\alpha)
!out    gg(1,1) =   0.5d0*dtermdoin/termin
!out    gg(1,2) = - 1.0d0/ber
!out    gg(1,3) = 2.0d0*radi*ain + 0.5d0*dtermdrin/termin
!    gg(1,1) = - 1.0d0/ber
!    gg(1,2) = 2.0d0*radi*ain + 0.5d0*dtermdrin/termin
    gg(1,1) = - 1.0d0/ber
    gg(1,2) = 2.0d0*radi*ain
!
!out    gg(2,1) =   0.5d0*dtermdoout/termout
!out    gg(2,2) = - 1.0d0/ber
!out    gg(2,3) = 2.0d0*radi*aout + 0.5d0*dtermdrout/termout
!
!out    gg(3,1) =   0.5d0*dtermdomid/termmid
!out    gg(3,2) = - 1.0d0/ber
!out    gg(3,3) = 2.0d0*radi*amid + 0.5d0*dtermdrmid/termmid
!    gg(2,1) = - 1.0d0/ber
!    gg(2,2) = 2.0d0*radi*amid + 0.5d0*dtermdrmid/termmid
    gg(2,1) = - 1.0d0/ber
    gg(2,2) = 2.0d0*radi*amid
!
!    test
!if (numero.le.10) then
! write (6,'(1p,3e12.4)') sgg(1), sgg(2)
! write (6,'(1p,3e12.4)') gg(1,1), gg(1,2)
! write (6,'(1p,3e12.4)') gg(2,1), gg(2,2)
!end if
!if (numero.eq.10) stop
!!stop
!!pause
!    test
!out    call minv(gg,sgg,3,4)
!out    ddome = sgg(1)
!out    ddber = sgg(2)
!out    ddradi = sgg(3)
    call minv(gg,sgg,2,4)
    ddber = sgg(1)
    ddradi = sgg(2)
    facfac = dmin1(dble(numero)/5.0d0,1.0d0)
!out    ome = ome + ddome*facfac
    ome = 0.0d0
    ber = ber + ddber*facfac
    radi = radi + ddradi*facfac
!out    error = dmax1(dabs(ddome/ome),dabs(ddber/ber),dabs(ddradi/radi))
    error = dmax1(dabs(ddber/ber),dabs(ddradi/radi))
    numero = numero + 1
    if (numero>1000) then
      write(6,*)' numero = ', numero, '   error =',error
    end if
!      if (iter<10.and.numero>10) exit
    if (numero>1010) exit
    if (dabs(error)<1.0d-08) exit
  end do
!    test
! write (6,*) "test values"
! write (6,*) sgg(1), sgg(2), sgg(3)
! write (6,*) ome, ber, radi
!    test
!      write(6,*)' numero = ', numero, '  error = ',error
!
  if (radi < 0.d0) write(6,*) ' ### radi minus ###'
  if (radi < 0.d0) radi = - radi
  if (ome < 0.d0) write(6,*) ' ### ome minus ###'
  if (ome < 0.d0) ome = - ome
  if (ber < 0.d0) write(6,*) ' ### ber minus ###'
  if (ber < 0.d0) ber = - ber
!out  ome = convf*ome + (1.d0-convf)*omeold
  ber = convf*ber + (1.d0-convf)*berold
  radi = convf*radi + (1.d0-convf)*radiold
!
!
! --  Improving alpha and psi.  
!
  alph = alph**((radi/radiold)**2)
  psi = psi**((radi/radiold)**2)
!##############################################
  ut = 1/sqrt(alph(nrf,ntgxy,0)**2 & 
     -psi(nrf,ntgxy,0)**4*ov2out)
  hh = ber*ut
  emdtest = 1.0d0/(pinx+1.0d0)*(hh-1.0d0)
!write (6,*) "h  = ", hh
!
end subroutine update_parameter_spherical_peos
