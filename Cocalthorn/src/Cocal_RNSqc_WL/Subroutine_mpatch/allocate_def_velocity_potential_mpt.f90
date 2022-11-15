subroutine allocate_def_velocity_potential_mpt
  use grid_parameter, only : nrf, ntf, npf, nlg
  use phys_constant, only : nmpt
  use def_velocity_potential_mpt
  use make_array_4d
  implicit none
!
  call alloc_array4d(vep, vep_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!
end subroutine allocate_def_velocity_potential_mpt
