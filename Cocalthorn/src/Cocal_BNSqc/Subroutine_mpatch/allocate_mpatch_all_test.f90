subroutine allocate_mpatch_all_test
!
   implicit none
!
  call allocate_def_quantities_derived_mpt
  call allocate_def_quantities_mpt
  call allocate_def_vector_phi_mpt
  call allocate_def_vector_x_mpt
  call allocate_grid_points_binary_excision_mpt
  call allocate_legendre_fn_grav_mpt
!  call allocate_metric_and_matter_mpt
  call allocate_radial_green_fn_grav_mpt
  call allocate_trigonometry_grav_phi_mpt
  call allocate_weight_midpoint_binary_excision_mpt
  call allocate_weight_midpoint_fluid_mpt
  call allocate_weight_midpoint_grav_mpt
!
end subroutine allocate_mpatch_all_test
 
