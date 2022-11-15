subroutine coordinate_patch_kit_grav
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use weight_midpoint_grav
  use weight_midpoint_fluid
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use radial_green_fn_grav
  use legendre_fn_grav
  implicit none
! call subroutines. the order is important.
  call read_parameter
  call grid_r
  call allocate_hgfn
  call calc_hgfn
  call grid_theta
  call trig_grav_theta
  call allocate_legendre
  call legendre
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call allocate_weight_midpoint_grav
  call weight_calc_midpoint_grav
!  call weight_calc_midpoint_grav_th4th
  call allocate_weight_midpoint_fluid
  call weight_calc_midpoint_fluid
  call grid_extended
end subroutine coordinate_patch_kit_grav

