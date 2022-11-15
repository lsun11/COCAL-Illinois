!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine allocate_weight_midpoint_fluid_sphcoord_mpt
  use weight_midpoint_fluid_sphcoord_mpt
  use grid_parameter, only : nrf, ntg, npg
  use phys_constant, only : nmpt
  use make_array_3d
  use make_array_4d
  implicit none
! weight binary excision
  call alloc_array4d(hwrtpg_fc_, 1, nrf, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array3d(hwtpg_fc_, 1, ntg, 1, npg, 1, nmpt)
!
end subroutine allocate_weight_midpoint_fluid_sphcoord_mpt
