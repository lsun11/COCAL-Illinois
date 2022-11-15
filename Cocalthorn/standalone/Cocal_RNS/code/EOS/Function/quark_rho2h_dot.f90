function quark_rho2h_dot(rho)
  use def_qeos_parameter
  implicit none
  real(8) :: rho, quark_rho2h_dot
  integer :: ii
  quark_rho2h_dot=0.0d0
  do ii=1,2*nphase
     quark_rho2h_dot=abchdot(ii)*rho**abihdot(ii)+quark_rho2h_dot
  end do
  return
END function quark_rho2h_dot
