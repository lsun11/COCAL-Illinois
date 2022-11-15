subroutine allocate_BHT_matter_emdrho
  use phys_constant, only : long
  use grid_parameter
  use def_matter, only : emdg, rhog
  use make_array_2d
  use make_array_3d
  implicit none
!
  call alloc_array3d(emdg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rhog, 0, nrg, 0, ntg, 0, npg)
!
end subroutine allocate_BHT_matter_emdrho
