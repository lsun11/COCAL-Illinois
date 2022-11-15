subroutine coordinate_patch_kit_grav_1D
  use grid_parameter_1D
  use coordinate_grav_r_1D
  use weight_grav_1D
  use radial_perm_fn_grav_1D
  implicit none
! call subroutines. the order is important.
  call read_parameter_1D
  call grid_r_1D
  call calc_hfsn
  call weight_calc_grav_1D
end subroutine coordinate_patch_kit_grav_1D
