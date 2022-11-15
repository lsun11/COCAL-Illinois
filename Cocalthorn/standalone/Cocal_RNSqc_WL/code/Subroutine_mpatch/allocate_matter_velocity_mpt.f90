subroutine allocate_matter_velocity_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf
  use def_matter_velocity
  use def_matter_velocity_mpt
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array4d(vxu_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(vyu_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(vzu_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!
end subroutine allocate_matter_velocity_mpt
