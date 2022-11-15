! --- tov equations are in this routine.
!
SUBROUTINE equation(x, y, f)
  use phys_constant, only : pi
  implicit none
!
  integer, parameter :: neq = 6
  real(8), intent(inout) :: y(neq), f(neq), x
  real(8) :: rr, rma, hh, bma, tma, psi, adm, &
          &  rho0, pre, ene, q
!
  rr  = x
  rma = y(1)
  hh  = y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
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
    f(1:6) = 0.0d0
  else
    f(1) = 4.0d0*pi*rr**2*ene
!c      f(2) = - (hh + 0.93d0)*(rma + 4.0d0*pi*pre*rr**3)
!c      f(2) = - (hh + 1.0d0)*(rma + 4.0d0*pi*pre*rr**3)
    f(2) = - hh*(rma + 4.0d0*pi*pre*rr**3) &
      &    /(rr**2 - 2.0d0*rma*rr)
    f(3) = 4.0d0*pi*rr**2*rho0/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(4) = 4.0d0*pi*rr**2*ene/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(5) = 0.5d0*psi/rr*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rr)**0.5d0)
    f(6) = 4.0d0*pi*rr**2*ene/(1.0d0-2.0d0*rma/rr)**0.5d0/psi
  end if
!
  return
END SUBROUTINE equation
