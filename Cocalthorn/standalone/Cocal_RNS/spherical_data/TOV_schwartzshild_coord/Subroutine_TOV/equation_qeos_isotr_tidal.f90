! --- tov equations are in this routine.
!
SUBROUTINE equation(x, y, f)
  use phys_constant, only : pi
  use def_qeos_parameter, only : rhosurf_gcm1
  implicit none
!
  integer, parameter :: ne = 12
  integer, parameter :: ell = 2
  real(8), intent(inout) :: y(ne), f(ne), x
  real(8) :: rr, rma, hh, bma, tma, psi, adm, c_s, dpde, &
          &  rho0, pre, ene, q, bcy, alp, rsh, rfac, bca, dpdrho0, dht, ht
  real(8), external :: quark_rho2p,  quark_rho2ene, quark_rho2h
!
  rr  = x
  rma = y(1)
!  hh  = y(2)
  rho0= y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
  bcy = y(7)
  alp = y(8)
  rsh = y(9)
  bca = y(10)
  ht  = y(11)
  dht = y(12)

  q = quark_rho2p(rho0)/rho0
  if (rho0 <= rhosurf_gcm1) then
    hh  = 1.0d0 !!!
    pre = 0.0d0
    rho0= 0.0d0
    ene = 0.0d0
  else
     call quark_rho2phenedpdrho(rho0,pre,hh,ene,dpdrho0)
     call quark_sound_speed(c_s,rho0)
     dpde = c_s*c_s
  end if
!
  if (rr == 1.0d-14) then
    rsh =1.0d-14
    f(1:10) = 0.0d0
    f(1) = 4.0d0*pi*ene*1.0d-28
    f(9)= 1.0d0
!    f(9)= 1.4254285581322195d0
    f(11) = 2.0d-14
  else
    rfac = rsh*(1.0d0-2.0d0*rma/rsh)**0.5d0/rr
    f(1) = 4.0d0*pi*rsh**2*ene*rfac
!c      f(2) = - (hh + 0.93d0)*(rma + 4.0d0*pi*pre*rr**3)
!c      f(2) = - (hh + 1.0d0)*(rma + 4.0d0*pi*pre*rr**3)
    
    f(2) = - (pre+ene)*(rma + 4.0d0*pi*pre*rsh**3)/((rsh**2 -2.0d0*rma*rsh)*dpdrho0)*rfac


    f(3) = 4.0d0*pi*rsh**2*rho0/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(4) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(5) = 0.5d0*psi/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(6) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0/psi*rfac

    f(7) = 0.5d0/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(8) = alp*(rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
    f(9) = rfac
    f(10)= (rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
    f(11)= dht*rfac
    f(12)= (-(2.0d0/rsh+(1.0d0/(1.0d0-2.0d0*rma/rsh))*(2.0d0*rma/(rsh**2)+4.0d0*pi*rsh*(pre-ene)))*dht -  &
    &        ((1.0d0/(1.0d0-2.0d0*rma/rsh))*(-6.0d0/rsh**2+4.0d0*pi*(ene+pre)/dpde +    &
!    &         ((1.0d0/(1.0d0-2.0d0*rma/rsh))*(-6.0d0/rsh**2+4.0d0*pi*dpde +    &
    &        4.0d0*pi*(5.0d0*ene+9.0d0*pre))  -      &
    &        ((2.0d0*(rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh))**2))*ht)*rfac

  end if
!
  return
END SUBROUTINE equation
