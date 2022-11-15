subroutine coordinate_update_binary
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  use def_vector_phi
  implicit none
! update  ex_nrg
  ex_nrg =   
  call calc_parameter_binary_excision
  call calc_grid_points_binary_excision
  call calc_weight_midpoint_binary_excision
!
  call calc_vector_phi_grav
  call calc_vector_phi_matter
end subroutine coordinate_update_binary
