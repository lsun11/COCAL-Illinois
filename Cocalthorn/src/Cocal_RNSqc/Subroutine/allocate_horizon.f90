subroutine allocate_horizon
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  use def_horizon, only : ahz
  use make_array_2d
  implicit none
!
  call alloc_array2d(ahz,0,ntg,0,npg)
!
end subroutine allocate_horizon
