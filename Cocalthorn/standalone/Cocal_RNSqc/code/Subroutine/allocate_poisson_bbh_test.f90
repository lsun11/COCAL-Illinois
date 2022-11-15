subroutine allocate_poisson_bbh_test
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_vector_x
  use def_vector_phi
  use def_vector_bh
  use def_vector_irg
  use make_array_2d
  use make_array_3d
  implicit none
!
  call alloc_array3d(psi, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_vector_bh
  call allocate_vector_irg
! 
  psi(0:nrg,0:ntg,0:npg) =1.0d0
  alph(0:nrg,0:ntg,0:npg)=1.0d0
  alps(0:nrg,0:ntg,0:npg)=1.0d0
!
end subroutine allocate_poisson_bbh_test
