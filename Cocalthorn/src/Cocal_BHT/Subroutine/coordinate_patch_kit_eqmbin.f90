subroutine coordinate_patch_kit_eqmbin
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  implicit none
! call subroutines. the order is important.
  call read_parameter_binary_excision
  call binary_excision_points
  call weight_calc_midpoint_binary_excision
end subroutine coordinate_patch_kit_eqmbin

