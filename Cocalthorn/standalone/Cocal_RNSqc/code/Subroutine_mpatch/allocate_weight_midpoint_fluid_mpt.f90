subroutine allocate_weight_midpoint_fluid_mpt
  use weight_midpoint_fluid_mpt
  use grid_parameter, only : nrf, ntf, npf
  use phys_constant, only : nmpt
  use make_array_4d 
  implicit none
!
  call alloc_array4d(hwrtpf_, 1, nrf, 1, ntf, 1, npf, 1, nmpt)
  call alloc_array4d(tzwrtpf_, 0, nrf, 1, ntf, 1, npf, 1, nmpt)
  call alloc_array4d(siwrtpf_, 0, nrf, 1, ntf, 1, npf, 1, nmpt)
  call alloc_array4d(rtsiwrtpf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!
end subroutine allocate_weight_midpoint_fluid_mpt