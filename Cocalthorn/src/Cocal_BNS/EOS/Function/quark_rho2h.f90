function quark_rho2h(rho)
  use def_qeos_parameter
  implicit none
  real(8) :: rho, quark_rho2h
  integer :: ii
  quark_rho2h=0.0d0
  do ii=1,2*nphase
     quark_rho2h=abch(ii)*rho**abih(ii)+quark_rho2h
  end do
     quark_rho2h=quark_rho2h+1.0d0+eneconst_gcm1
  return
END function quark_rho2h
