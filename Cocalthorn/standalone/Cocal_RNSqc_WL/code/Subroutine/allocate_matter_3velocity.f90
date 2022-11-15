subroutine allocate_matter_3velocity
  use phys_constant, only : long
  use grid_parameter
  use def_matter_velocity, only  : vxu, vyu, vzu
  use make_array_3d
  implicit none
!
  call alloc_array3d(vxu, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(vyu, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(vzu, 0,nrf, 0,ntf, 0,npf)
!
end subroutine allocate_matter_3velocity
