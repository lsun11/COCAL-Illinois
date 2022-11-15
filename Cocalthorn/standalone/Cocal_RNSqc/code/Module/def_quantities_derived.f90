module def_quantities_derived
  use phys_constant, only : long
  implicit none
! -- Other derived quantities
  real(long) :: coord_radius_normalized_x
  real(long) :: coord_radius_normalized_y
  real(long) :: coord_radius_normalized_z
  real(long) :: coord_axis_ratio_yx, proper_axis_ratio_yx
  real(long) :: coord_axis_ratio_zx, proper_axis_ratio_zx
  real(long) :: coord_eccentricity_zx, proper_eccentricity_zx
  real(long) :: omega, omega_M, omega2_M2
  real(long) :: E_bind_over_M, J_over_M2, J_over_Madm2
  real(long) :: Madm_MK_virial, M0_minus_M0sph_normalized
  real(long) :: rho_c_normalized, pre_c_normalized
  real(long) :: epsilon_c_normalized, p_over_rho_c_normalized
!
end module def_quantities_derived
