subroutine allocate_metric_BH_WL
  use grid_parameter, only : nrg, ntg, npg
  use def_vector_x
  use def_vector_phi
  use def_vector_bh
  use def_matter, only : omeg
  use make_array_3d
  implicit none
!
  call allocate_metric_WL
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_vector_bh
!
  call alloc_array3d(omeg, 0, nrg, 0, ntg, 0, npg)
!
end subroutine allocate_metric_BH_WL
