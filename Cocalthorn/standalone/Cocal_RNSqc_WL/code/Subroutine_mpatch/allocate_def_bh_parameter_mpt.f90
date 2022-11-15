subroutine allocate_def_bh_parameter_mpt
  use phys_constant, only : nmpt
  use def_bh_parameter_mpt
  use make_array_2d
  use make_char2_array_2d
  implicit none
!
  call alloc_array2d(def_bh_param_real_ , 1, 50, 1, nmpt)
  call alloc_char2_array2d(def_bh_param_char2_ , 1, 5, 1, nmpt)
!
end subroutine allocate_def_bh_parameter_mpt
