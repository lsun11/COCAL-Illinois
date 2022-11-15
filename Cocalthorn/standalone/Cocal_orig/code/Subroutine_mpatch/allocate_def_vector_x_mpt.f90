subroutine allocate_def_vector_x_mpt
  use def_vector_x_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use make_array_5d
  implicit none
  call alloc_array5d(vec_xg_, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, nmpt)
  call alloc_array5d(vec_xf_, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, nmpt)
  call alloc_array5d(hvec_xg_, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, nmpt)
  call alloc_array5d(hvec_xf_, 1, nrf, 1, ntf, 1, npf, 1, 3, 1, nmpt)
end subroutine allocate_def_vector_x_mpt