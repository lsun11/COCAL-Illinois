! --- tov equations are in this routine.
!
SUBROUTINE equation_isotr(x, y, f, rfac0)
  use phys_constant, only : pi
  implicit none
!
  integer, parameter :: ne = 10
  real(8), intent(in) :: rfac0
  real(8), intent(inout) :: y(ne), f(ne), x
  real(8) :: rr, rma, hh, bma, tma, psi, adm, &
          &  rho0, pre, ene, q, bcy, alp, rsh, rfac, bca
!
  rr  = x
  rma = y(1)
  hh  = y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
  bcy = y(7)
  alp = y(8)
  rsh = y(9)
  bca = y(10)

!  if (q <= 0.0d0) then
  if (hh <= 1.0d0) then
    q = 0.0d0 !!!
    pre = 0.0d0
    rho0= 0.0d0
    ene = 0.0d0
  else
    call peos_h2qprho(hh, q, pre, rho0, ene)
  end if
!
  if (rr == 0.0d0) then
    f(1:10) = 0.0d0
    f(9)=rfac0
!    f(9)=1.4255d0
  else
    rfac = rsh*(1.0d0-2.0d0*rma/rsh)**0.5d0/rr
    f(1) = 4.0d0*pi*rsh**2*ene*rfac
!c      f(2) = - (hh + 0.93d0)*(rma + 4.0d0*pi*pre*rr**3)
!c      f(2) = - (hh + 1.0d0)*(rma + 4.0d0*pi*pre*rr**3)
    f(2) = - hh*(rma + 4.0d0*pi*pre*rsh**3) &
      &    /(rsh**2 - 2.0d0*rma*rsh)*rfac
    f(3) = 4.0d0*pi*rsh**2*rho0/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(4) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(5) = 0.5d0*psi/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(6) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0/psi*rfac

    f(7) = 0.5d0/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(8) = alp*(rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
    f(9) = rfac
!    f(10)= 1.0d0/rsh/(1.0d0-2.0d0*rma/rsh)**0.5
    f(10)= (rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
  end if
!
  return
END SUBROUTINE equation_isotr
