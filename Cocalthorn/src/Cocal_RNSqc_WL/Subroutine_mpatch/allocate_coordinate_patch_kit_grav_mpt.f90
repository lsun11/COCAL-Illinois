subroutine allocate_coordinate_patch_kit_grav_mpt
  use weight_midpoint_grav
  use weight_midpoint_fluid
  use trigonometry_grav_phi
  use legendre_fn_grav
  use weight_midpoint_fluid_sphcoord
  implicit none
! call subroutines. the order is important.
  call allocate_legendre
  call allocate_trig_grav_mphi
  call allocate_weight_midpoint_grav
  call allocate_weight_midpoint_fluid
  call allocate_weight_midpoint_fluid_sphcoord
end subroutine allocate_coordinate_patch_kit_grav_mpt
