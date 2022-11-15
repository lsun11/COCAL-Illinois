subroutine allocate_BHT_metric_and_matter_WL
  use phys_constant, only : long
  use grid_parameter
  use def_matter
  use def_vector_x
  use def_vector_phi
  use make_array_2d
  use make_array_3d
  implicit none
!
  call allocate_BHT_metric_WL
!  call allocate_metric_on_SFC_WL
  call allocate_BHT_matter_emdrho
  call allocate_BHT_matter_4velocity
  call allocate_BHT_matter_rotlaw
!  call alloc_array3d(hhf, 0, nrg, 0, ntg, 0, npg)
!  call alloc_array3d(lambda, 0, nrg, 0, ntg, 0, npg)
!
  call allocate_BHT_vector_x
  call allocate_BHT_vector_phi
!
end subroutine allocate_BHT_metric_and_matter_WL
