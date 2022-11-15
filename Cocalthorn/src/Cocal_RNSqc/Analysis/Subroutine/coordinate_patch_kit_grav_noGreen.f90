subroutine coordinate_patch_kit_grav_noGreen
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use weight_midpoint_grav
  implicit none
! call subroutines. the order is important.
  call read_parameter
  call grid_r
!  call grid_r_bhex('eBH')
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call grid_extended
!
end subroutine coordinate_patch_kit_grav_noGreen
