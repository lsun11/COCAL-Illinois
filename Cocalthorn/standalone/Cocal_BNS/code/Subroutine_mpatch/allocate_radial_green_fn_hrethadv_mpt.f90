subroutine allocate_radial_green_fn_hrethadv_mpt
  use radial_green_fn_hrethadv_mpt
  use grid_parameter, only : nrg, nlg
  use phys_constant, only : nmpt
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array5d(bsjy_, 1, nrg, 0, nlg, 0, nlg, 0, nrg, 1, nmpt)
  call alloc_array4d(sbsjy_, 0, nlg, 0, nlg, 0, nrg, 1, nmpt)
  call alloc_array4d(sbsjyp_, 0, nlg, 0, nlg, 0, nrg, 1, nmpt)
!
end subroutine allocate_radial_green_fn_hrethadv_mpt
