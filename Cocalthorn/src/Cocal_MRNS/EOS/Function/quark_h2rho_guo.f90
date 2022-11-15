function quark_h2rho_guo(hh)
  use def_qeos_parameter
  implicit none
  real(8) :: hh, quark_h2rho_guo
  integer :: ii
     quark_h2rho_guo=((hh-1.0d0-eneconst_gcm1)/(abch(1)+abch(1+nphase)))**&
&(1.0d0/abih(1))
  return
END function quark_h2rho_guo
