subroutine update_parameter_triaxial_WL_peos(convf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, nrf_deform, nrf, &
  &                          ntgpolp, ntgeq, ntgxy, npgxzp, npgyzp
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd, tfkij, tfkijkij
  use def_metric_hij, only : hxxd, hyyd
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
  real(long)     ::    vphixin, vphiymid, vphiyout
  real(long)     ::    loghc
  real(long)     ::    radiold, berold, omeold
  real(long)     ::    ddradi, ddber, ddome
  real(long)     ::    facfac, numero, error, hc, pre, rho, ene
  real(long)     ::    gin, gmid, gout
  integer        ::    nin
!
! ################################################################       
  real(long)     ::    ut, hh, emdtest
! ################################################################       
  nin = nrf_deform
!
  sgg = 0.0d0
  gg  = 0.0d0
!
  call peos_q2hprho(emdc, hc, pre, rho, ene)
  loghc = dlog(hc)
  vphixin =  vec_phig(nin,ntgeq,npgyzp,1)
  vphiymid = vec_phig(0,0,0,2)
  vphiyout = vec_phig(nrf,ntgeq,0,2)
  omeold = ome
  berold = ber
  radiold = radi
  numero = 1
! ################################################################ 
  ain = dlog(alph(nin,ntgeq,npgyzp)**(1.d0/radi**2))
  pain= (psi(nin,ntgeq,npgyzp)**2/alph(nin,ntgeq,npgyzp))**(1.d0/radi**2)
  bin = bvxd(nin,ntgeq,npgyzp)
  gin = 1.0d0+hxxd(nin,ntgeq,npgyzp)
!
  amid = dlog(alph(0,0,0)**(1.d0/radi**2))
  pamid= (psi(0,0,0)**2/alph(0,0,0))**(1.d0/radi**2)
  bmid = bvyd(0,0,0)
  gmid = 1.0d0+hyyd(0,0,0)
!
  aout = dlog(alph(nrf,ntgeq,0)**(1.d0/radi**2))
  paout= (psi(nrf,ntgeq,0)**2/alph(nrf,ntgeq,0))**(1.d0/radi**2)
  bout = bvyd(nrf,ntgeq,0)
  gout = 1.0d0+hyyd(nrf,ntgeq,0)
! ################################################################
  do
    ovin  = bin  + ome*vphixin
    ovmid = bmid + ome*vphiymid
    ovout = bout + ome*vphiyout
!
    ov2in  = gin*ovin**2
    ov2mid = gmid*ovmid**2
    ov2out = gout*ovout**2
    termin  = 1.0d0 - pain**( 2.d0*radi**2)*ov2in
    termmid = 1.0d0 - pamid**(2.d0*radi**2)*ov2mid
    termout = 1.0d0 - paout**(2.d0*radi**2)*ov2out
! termin is 1-(\phi^4/\alpha^2)^((R_0/R_i)^2) \omega^2            
    dtermdoin = - pain**(2.d0*radi**2)*2.0d0*gin*ovin*vphixin
    dtermdomid= - pamid**(2.d0*radi**2)*2.0d0*gmid*ovmid*vphiymid
    dtermdoout= - paout**(2.d0*radi**2)*2.0d0*gout*ovout*vphiyout
    dtermdrin =-dlog(pain)*4.d0*radi*pain**(2.d0*radi**2)*ov2in
    dtermdrmid=-dlog(pamid)*4.d0*radi*pamid**(2.d0*radi**2)*ov2mid
    dtermdrout=-dlog(paout)*4.d0*radi*paout**(2.d0*radi**2)*ov2out
!
    sgg(1) = -(radi**2*ain  + 0.5d0*dlog(termin)  - dlog(ber))
    sgg(2) = -(radi**2*aout + 0.5d0*dlog(termout) - dlog(ber))
    sgg(3) = -(radi**2*amid + 0.5d0*dlog(termmid) - dlog(ber) + loghc)
!
!  ain is ln(\alpha)
    gg(1,1) =   0.5d0*dtermdoin/termin
    gg(1,2) = - 1.0d0/ber
    gg(1,3) = 2.0d0*radi*ain + 0.5d0*dtermdrin/termin
!
    gg(2,1) =   0.5d0*dtermdoout/termout
    gg(2,2) = - 1.0d0/ber
    gg(2,3) = 2.0d0*radi*aout + 0.5d0*dtermdrout/termout
!
    gg(3,1) =   0.5d0*dtermdomid/termmid
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
    if (dabs(error)<1.0d-08) exit
  end do
!    test
!write (6,*) "values"
!write (6,*) omeold, berold, radiold
!write (6,*) ome, ber, radi
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
! --  Improving alpha and psi.  
!
  alph = alph**((radi/radiold)**2)
  psi = psi**((radi/radiold)**2)
!
end subroutine update_parameter_triaxial_WL_peos
