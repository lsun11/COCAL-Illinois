subroutine coordinate_patch_kit_grav_mpt
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use weight_midpoint_grav
  use weight_midpoint_fluid
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use legendre_fn_grav
  implicit none
! call subroutines. the order is important.
  call grid_r
!  call grid_r_bhex('eBH')
  call grid_theta
  call trig_grav_theta
  call legendre
  call grid_phi
  call trig_grav_phi
  call weight_calc_midpoint_grav
  call weight_calc_midpoint_grav_th4th
  call weight_calc_midpoint_fluid
  call grid_extended
end subroutine coordinate_patch_kit_grav_mpt

