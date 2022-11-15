subroutine allocate_BHT_matter_WL
  use phys_constant, only : long
  use grid_parameter
  use def_matter
  use def_vector_x
  use def_vector_phi
  implicit none
!
  call allocate_BHT_matter_emdrho
  call allocate_BHT_matter_4velocity
  call allocate_BHT_matter_rotlaw
!
  call allocate_BHT_vector_x
  call allocate_BHT_vector_phi
!
end subroutine allocate_BHT_matter_WL
