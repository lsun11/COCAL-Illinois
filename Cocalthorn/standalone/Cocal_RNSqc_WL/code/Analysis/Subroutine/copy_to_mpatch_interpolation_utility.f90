subroutine copy_to_mpatch_interpolation_utility(impt)
!
implicit none
  integer :: impt
!
  call copy_grid_parameter_to_mpt(impt)
  call copy_grid_parameter_binary_excision_to_mpt(impt)
  call copy_def_peos_parameter_to_mpt(impt)
  call copy_def_matter_parameter_to_mpt(impt)
  call copy_grid_points_binary_excision_to_mpt(impt)
  call copy_coordinate_grav_r_to_mpt(impt)
  call copy_coordinate_grav_theta_to_mpt(impt)
  call copy_coordinate_grav_phi_to_mpt(impt)
  call copy_coordinate_grav_extended_to_mpt(impt)
  call copy_def_binary_parameter_to_mpt(impt)
  call copy_def_quantities_BNS_to_mpt(impt)
  call copy_def_quantities_derived_to_mpt(impt)
  call copy_def_vector_x_to_mpt(impt)
  call copy_def_vector_phi_to_mpt(impt)
  call copy_trigonometry_grav_theta_to_mpt(impt)
  call copy_trigonometry_grav_phi_to_mpt(impt)
  call copy_weight_midpoint_grav_to_mpt(impt)
  call copy_weight_midpoint_fluid_to_mpt(impt)
  call copy_weight_midpoint_binary_excision_to_mpt(impt)
!
end subroutine copy_to_mpatch_interpolation_utility
 
