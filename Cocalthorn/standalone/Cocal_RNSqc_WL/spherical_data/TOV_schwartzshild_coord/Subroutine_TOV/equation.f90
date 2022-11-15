! --------------------------------------------------------
! --- equations are in this routine.
!
subroutine equation(x, y, f)
  use def_matter, only : pinx
  use phys_constant, only : pi
  implicit none
!
  integer, parameter     :: neq = 6
  real(8), intent(inout) :: y(neq), f(neq), x
  real(8) :: rr, rma, emd, bma, tma, psi, adm, &
          &  rho0, rho, pre, hhh
!
  rr  = x
  rma = y(1)
  emd = y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
  if (emd <= 0.0d0) then
    rho0= 0.0d0
    rho = 0.0d0
    pre = 0.0d0
    hhh = 1.0d0
  else
    rho0= emd**pinx
    rho = emd**pinx*(1.0d0 + pinx*emd)
    pre = emd**(1.0d0 + pinx)
    hhh = 1.0d0 + (1.0d0+pinx)*emd
  end if
!
  if (rr == 0.0d0) then
    f(1:6) = 0.0d0
  else
    f(1) = 4.0d0*pi*rr**2*rho
    f(2) = - 1.0d0/(1.0d0+pinx)*hhh*(rma + 4.0d0*pi*pre*rr**3)/ &
      &    (rr**2 - 2.0d0*rma*rr)
    f(3) =  4.0d0*pi*rr**2*rho0/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(4) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(5) =  0.5d0*psi/rr*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rr)**0.5d0)
    f(6) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0/psi
  end if
!
  return
end subroutine equation
