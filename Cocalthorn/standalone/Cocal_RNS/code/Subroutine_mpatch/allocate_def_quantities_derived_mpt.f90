subroutine allocate_def_quantities_derived_mpt
  use phys_constant, only : nmpt
  use def_quantities_derived_mpt
  use make_array_2d
  implicit none
! -- Other derived quantities
  call alloc_array2d(def_quantities_derived_real_ , 1, 50, 1, nmpt)
!
end subroutine allocate_def_quantities_derived_mpt
