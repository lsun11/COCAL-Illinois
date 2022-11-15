subroutine coordinate_patch_cartesian
  use grid_parameter_cartesian
  use coordinate_grav_xyz
  implicit none
! call subroutines. the order is important.
  call read_parameter_cartesian
  call grid_xyz
end subroutine coordinate_patch_cartesian

