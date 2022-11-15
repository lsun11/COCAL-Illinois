subroutine allocate_BBH_CF_AH
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_metric_excurve_grid
  use def_vector_x
  use def_vector_phi
  use def_vector_bh
  use def_vector_irg
  use make_array_3d
  use make_array_5d
  implicit none
!
  call alloc_array3d(psi , 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzu, 0, nrg, 0, ntg, 0, npg)
!
  call alloc_array3d(tfkijkij,1,nrg,1,ntg,1,npg)
  call alloc_array5d(tfkij,1,nrg,1,ntg,1,npg,1,3,1,3)
  call alloc_array3d(tfkijkij_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array5d(tfkij_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
  call alloc_array3d(trk,1,nrg,1,ntg,1,npg)
  call alloc_array3d(trk_grid, 0, nrg, 0, ntg, 0, npg)
!  
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_vector_bh
  call allocate_vector_irg
!
end subroutine allocate_BBH_CF_AH
