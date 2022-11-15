!  weight for numerical integration using mid-point rule
!______________________________________________
subroutine allocate_weight_midpoint_binary_excision_mpt
  use weight_midpoint_binary_excision_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use make_array_3d
  use make_array_4d
  implicit none
! weight binary excision
  call alloc_array4d(hwrtpg_ex_, 1, nrg, 1, ntg, 1, npg, 1, nmpt)
  call alloc_array3d(hwtpg_ex_, 1, ntg, 1, npg, 1, nmpt)
!
end subroutine allocate_weight_midpoint_binary_excision_mpt
