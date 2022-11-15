subroutine allocate_def_quantities_bh_mpt
  use phys_constant, only : nmpt
  use def_quantities_bh_mpt
  use make_array_2d
  implicit none
! -- Other bh quantities
  call alloc_array2d(def_quantities_bh_ , 1, 10, 1, nmpt)
!
end subroutine allocate_def_quantities_bh_mpt
