subroutine allocate_metric_and_matter_qeos
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_metric_excurve_grid
  use def_matter
  use def_vector_x
  use def_vector_phi
  use make_array_2d
  use make_array_3d
  use make_array_5d
  implicit none
!
  call alloc_array2d(rs, 0, ntf, 0, npf)
  call alloc_array2d(ergoin, 0, ntg, 0, npg)
  call alloc_array2d(ergoout, 0, ntg, 0, npg)
  call alloc_array3d(rhog, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rhof, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(utg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(utf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(omeg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(omef, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg_int, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef_int, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(psi, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(tfkijkij, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array5d(tfkij, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 3)
  call alloc_array5d(tfkij_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
  call alloc_array3d(tfkijkij_grid, 0, nrg, 0, ntg, 0, npg)
!
  call allocate_vector_x
  call allocate_vector_phi
!
end subroutine allocate_metric_and_matter_qeos
