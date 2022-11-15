! --- tov equations are in this routine.
!
SUBROUTINE equation(x, y, f)
  use phys_constant, only : pi
  use def_qeos_parameter, only : enesurf_gcm1, aq
  implicit none
!
  integer, parameter :: ne = 10
  real(8), intent(inout) :: y(ne), f(ne), x
  real(8) :: rr, rma, hh, bma, tma, psi, adm, &
          &  rho0, pre, ene, q, bcy, alp, rsh, rfac, bca
!
  rr  = x
  rma = y(1)
  ene = y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
  bcy = y(7)
  alp = y(8)
  rsh = y(9)
  bca = y(10)

!  if (q <= 0.0d0) then
  if (ene <= enesurf_gcm1) then
    pre = 0.0d0
    rho0= 0.0d0
    ene = 0.0d0
  else
    call mit_ene2prho(ene, pre, rho0)
  end if
!
  if (rr == 0.0d0) then
    f(1:10) = 0.0d0
    f(9)= 1.0d0
!    f(9)= 1.4254285581322195d0
  else
    rfac = rsh*(1.0d0-2.0d0*rma/rsh)**0.5d0/rr
    f(1) = 4.0d0*pi*rsh**2*ene*rfac
!c      f(2) = - (hh + 0.93d0)*(rma + 4.0d0*pi*pre*rr**3)
!c      f(2) = - (hh + 1.0d0)*(rma + 4.0d0*pi*pre*rr**3)
     f(2) = - ((aq+1.0d0)*ene-aq*enesurf_gcm1)*(rma + 4.0d0*pi*aq*(ene-enesurf_gcm1)*rsh**3)/(rsh**2 - &
      &       2.0d0*rma*rsh)/aq*rfac

    f(3) = 4.0d0*pi*rsh**2*rho0/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(4) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0*rfac
    f(5) = 0.5d0*psi/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(6) = 4.0d0*pi*rsh**2*ene/(1.0d0-2.0d0*rma/rsh)**0.5d0/psi*rfac

    f(7) = 0.5d0/rsh*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rsh)**0.5d0)*rfac
    f(8) = alp*(rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
    f(9) = rfac
    f(10)= (rma+4.0d0*pi*rsh**3*pre)/(rsh**2-2.0d0*rma*rsh)*rfac
  end if
!
  return
END SUBROUTINE equation
