subroutine allocate_def_vector_phi_mpt
  use def_vector_phi_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use make_array_4d
  use make_array_5d
  implicit none
  call alloc_array5d(vec_phig_, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, nmpt)
  call alloc_array5d(vec_phif_, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, nmpt)
  call alloc_array5d(hvec_phig_, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, nmpt)
  call alloc_array5d(hvec_phif_, 1, nrf, 1, ntf, 1, npf, 1, 3, 1, nmpt)
  call alloc_array4d(hvec_phif_surface_, 1, ntf, 1, npf, 1, 3, 1, nmpt)
end subroutine allocate_def_vector_phi_mpt