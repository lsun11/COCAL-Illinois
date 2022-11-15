subroutine allocate_metric_and_matter_BHNS_test
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_metric
  use def_matter
  use def_vector_x
  use def_vector_phi
  use def_vector_bh
  use def_vector_irg
  use make_array_3d
  use make_array_2d
  implicit none
!
  call alloc_array2d(rs, 0, ntf, 0, npf)
  call alloc_array3d(emd, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(emdg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(utf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(psi, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
!
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_vector_bh
  call allocate_vector_irg
!
end subroutine allocate_metric_and_matter_BHNS_test
