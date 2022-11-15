function quark_rho2ene(rho)
  use def_qeos_parameter
  implicit none
  real(8) :: rho, quark_rho2ene
  integer :: ii
  quark_rho2ene=0.0d0
  do ii=1,nphase
     quark_rho2ene=abcene(ii)*rho**abiene(ii)+quark_rho2ene
  end do
     quark_rho2ene=quark_rho2ene+rho+rho*eneconst_gcm1
  return
END function quark_rho2ene
