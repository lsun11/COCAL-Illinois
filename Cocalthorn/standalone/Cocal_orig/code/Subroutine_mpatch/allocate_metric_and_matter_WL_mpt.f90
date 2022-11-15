subroutine allocate_metric_and_matter_WL
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_matter
  use def_gamma_crist
  use def_gamma_crist_grid
  use def_metric_hij
  use def_metric_rotshift
  use def_ricci_tensor
  use def_metric_excurve_grid
  use def_shift_derivatives
  use def_shift_derivatives_grid
  use def_Lie_derivatives
  use def_Lie_derivatives_grid
  use def_matter_parameter
  use def_cristoffel
  use def_cristoffel_grid
  use def_matter
  use def_vector_x
  use def_vector_phi
  use make_array_2d
  use make_array_4d
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array4d(rs, 0, ntf, 0, npf)
  call alloc_array4d(em_d, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(utf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(lambda_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(psi_, 0, nrg, 0, ntg, 0, npg,1,nmpt, 1, nmpt)
  call alloc_array4d(alph_, 0, nrg, 0, ntg, 0, npg,1,nmpt, 1, nmpt)
  call alloc_array4d(bvxd_, 0, nrg, 0, ntg, 0, npg,1,nmpt, 1, nmpt)
  call alloc_array4d(bvyd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvzd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
!
  call alloc_array4d(alps2_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvxu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvyu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(bvzu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxxd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxyd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxzd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hyyd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hyzd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hzzd_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxxu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxyu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hxzu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hyyu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hyzu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(hzzu_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)

!
  call allocate_vector_x
  call allocate_vector_phi
!
end subroutine allocate_metric_and_matter_WL
