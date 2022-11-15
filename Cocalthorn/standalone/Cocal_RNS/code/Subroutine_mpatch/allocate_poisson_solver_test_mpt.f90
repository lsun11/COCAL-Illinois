subroutine allocate_poisson_solver_test_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : ntf, nrf, npf, nrg, ntg, npg
  use def_metric_mpt
  use def_matter_mpt
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(rs_, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emd_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emdg_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(psi_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(alph_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
!
end subroutine allocate_poisson_solver_test_mpt
