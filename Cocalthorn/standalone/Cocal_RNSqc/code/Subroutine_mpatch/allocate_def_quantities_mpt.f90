subroutine allocate_def_quantities_mpt
  use phys_constant, only : nmpt
  use def_quantities_mpt
  use make_array_2d
  implicit none
!
  call alloc_array2d(def_quantities_real_ , 1, 200, 1, nmpt)
!
end subroutine allocate_def_quantities_mpt