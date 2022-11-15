subroutine copy_from_mpatch_all_BNS_CF(impt)
!
  implicit none
  integer :: impt
!
  call copy_grid_parameter_from_mpt(impt)
  call copy_grid_parameter_binary_excision_from_mpt(impt)
  call copy_def_peos_parameter_from_mpt(impt)
  call copy_def_matter_parameter_from_mpt(impt)
!  call copy_def_bh_parameter_from_mpt(impt)
  call copy_grid_points_binary_excision_from_mpt(impt)
  call copy_coordinate_grav_r_from_mpt(impt)
  call copy_coordinate_grav_theta_from_mpt(impt)
  call copy_coordinate_grav_phi_from_mpt(impt)
  call copy_coordinate_grav_extended_from_mpt(impt)
  call copy_def_binary_parameter_from_mpt(impt)
  call copy_def_quantities_BNS_from_mpt(impt)
!  call copy_def_quantities_ns_from_mpt(impt)
  call copy_def_quantities_derived_from_mpt(impt)
  call copy_def_vector_x_from_mpt(impt)
  call copy_def_vector_phi_from_mpt(impt)
!  call copy_def_vector_bh_from_mpt(impt)
  call copy_legendre_fn_grav_from_mpt(impt)
  call copy_trigonometry_grav_theta_from_mpt(impt)
  call copy_trigonometry_grav_phi_from_mpt(impt)
  call copy_weight_midpoint_grav_from_mpt(impt)
  call copy_weight_midpoint_fluid_from_mpt(impt)
  call copy_weight_midpoint_binary_excision_from_mpt(impt)
!
end subroutine copy_from_mpatch_all_BNS_CF
