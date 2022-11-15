subroutine allocate_horizon_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : ntg, npg
  use def_horizon_mpt, only : ahz_
  use make_array_3d
  implicit none
!
  call alloc_array3d(ahz_, 0,ntg, 0,npg, 1,nmpt)
!
end subroutine allocate_horizon_mpt
