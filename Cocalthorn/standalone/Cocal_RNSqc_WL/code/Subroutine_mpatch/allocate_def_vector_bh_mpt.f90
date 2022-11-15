subroutine allocate_def_vector_bh_mpt
  use def_vector_bh_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : ntg, npg
  use make_array_4d
  implicit none
  call alloc_array4d( vec_bh_cm_xg_, 0,ntg, 0,npg, 1,3, 1,nmpt)
  call alloc_array4d(hvec_bh_cm_xg_, 1,ntg, 1,npg, 1,3, 1,nmpt)
  call alloc_array4d( vec_bh_cm_phig_, 0,ntg, 0,npg, 1,3, 1,nmpt)
  call alloc_array4d(hvec_bh_cm_phig_, 1,ntg, 1,npg, 1,3, 1,nmpt)
  call alloc_array4d( vec_bh_cbh_xg_, 0,ntg, 0,npg, 1,3, 1,nmpt)
  call alloc_array4d(hvec_bh_cbh_xg_, 1,ntg, 1,npg, 1,3, 1,nmpt)
  call alloc_array4d( vec_bh_cbh_phig_, 0,ntg, 0,npg, 1,3, 1,nmpt)
  call alloc_array4d(hvec_bh_cbh_phig_, 1,ntg, 1,npg, 1,3, 1,nmpt)
  call alloc_array4d( vec_bh_cbh_spin_, 0,ntg, 0,npg, 1,3, 1,nmpt)
  call alloc_array4d(hvec_bh_cbh_spin_, 1,ntg, 1,npg, 1,3, 1,nmpt)
end subroutine allocate_def_vector_bh_mpt
