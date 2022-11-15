subroutine allocate_metric_and_matter_WL
  use phys_constant, only : long
  use grid_parameter
  use def_matter
  use def_vector_x
  use def_vector_phi
  use make_array_2d
  use make_array_3d
  implicit none
!
  call allocate_metric_WL
  call allocate_metric_on_SFC_WL
  call allocate_matter_emdrsrho
  call allocate_matter_4velocity
  call allocate_matter_rotlaw
  call alloc_array3d(hhf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(lambda, 0, nrf, 0, ntf, 0, npf)
!
  call allocate_vector_x
  call allocate_vector_phi
!
end subroutine allocate_metric_and_matter_WL
