subroutine peos_paimproco_frg(ber,radi,convf,iter,fmax0,isw)
!
! --- Improve parameters
!
!  use CB_fR_mesh_grav, only : nrgsf
  use def_metric_1D, only : alph, psi
  use CB_fR_param_physp, only : surr
  use grid_parameter_1D, only : nrf, emdc, pinx, nrg
  implicit none
  real(8) :: berb, dlogamid, dlogasf, dloghmid, &
  &          dloghsf, f0, f1, hhmid, hhsf, radib, pre, rho, ene
  real(8), intent(inout) :: ber, radi, convf, fmax0
  integer, intent(in)    :: iter, isw
!
! --- for interporation.
!
!  nrgsf = nrf
!
  call peos_q2hprho(emdc, hhmid, pre, rho, ene)
!  hhmid  = 1.d0+(pinx+1.d0)*emdc
  hhsf   = 1.d0!+(pinx+1.d0)*surr
  dloghmid = dlog(hhmid)
  dloghsf  = dlog(hhsf)
!
  f0 = alph(0)
  f1 = alph(nrf)
  dlogamid = dlog(f0**(1.d0/radi**2))
  dlogasf = dlog(f1**(1.d0/radi**2))
!
  berb = ber
  radib = radi
!
  radi = dsqrt(dabs(-(dloghmid - dloghsf)/(dlogamid - dlogasf)))
  ber  = dexp(0.5d0*(radi**2*(dlogamid + dlogasf)&
     &                         + (dloghmid + dloghsf)))
!
!
  if (radi <= 0.d0) then
    write(6,*) ' ### radi minus ###'
    radi = - radi
  end if
  if (ber <= 0.d0) then
    write(6,*) ' ### ber minus ###'
    ber = - ber
  end if
  ber = convf*ber + (1.d0-convf)*berb
  radi = convf*radi + (1.d0-convf)*radib
!
  if (isw == 1) then
    write(6,"(a22,es12.4,',  value =',es12.4)")'  -- ber   --, error =', &
       &        2.0d0*abs((ber - berb)/(ber + berb)), ber
    write(6,"(a22,es12.4,',  value =',es12.4)")'  -- radi  --, error =', &
       &        2.0d0*dabs((radi - radib)/(radi + radib)), radi
  end if
  fmax0 = dmax1(2.0d0*dabs((ber - berb)/(ber + berb)),&
     &              2.0d0*dabs((radi - radib)/(radi + radib)),fmax0)
!
! --  Improving alpha and psi.  
!
  alph(0:nrg) = alph(0:nrg)**((radi/radib)**2)
  psi(0:nrg) = psi(0:nrg)**((radi/radib)**2)
!
end subroutine peos_paimproco_frg
