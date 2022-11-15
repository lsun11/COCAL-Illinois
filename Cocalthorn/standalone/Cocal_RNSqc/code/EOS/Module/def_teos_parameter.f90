module def_peos_parameter
  use phys_constant         !nnteos
  implicit none
  real(8) :: enei(0:nnteos), prei(0:nnteos), rhoi(0:nnteos), qi(0:nnteos), hi(0:nnteos)
  real(8) :: rhocgs(0:nnteos), precgs(0:nnteos), enecgs(0:nnteos)
  real(8) :: rhoini_cgs, rhoini_gcm1, emdini_gcm1  !used in TOV solver
  real(8) :: kappa_crust, gamma_crust 
  integer :: nphase
end module def_peos_parameter
