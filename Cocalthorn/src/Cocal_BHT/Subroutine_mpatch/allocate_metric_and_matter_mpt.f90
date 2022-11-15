subroutine allocate_metric_and_matter_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use def_vector_x
  use def_vector_phi
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(rs_, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emd_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(utf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(psi_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(alph_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvxd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvyd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvzd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
!
  call allocate_def_vector_x_mpt
  call allocate_def_vector_phi_mpt
!
end subroutine allocate_metric_and_matter_mpt
