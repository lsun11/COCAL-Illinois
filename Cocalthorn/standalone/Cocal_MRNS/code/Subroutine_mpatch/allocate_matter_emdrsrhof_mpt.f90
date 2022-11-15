subroutine allocate_matter_emdrsrhof_mpt
  use grid_parameter , only : nrf, ntf, npf
  use matter_emdrsrhof_mpt
  use phys_constant, only : nmpt
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(rs_, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emd_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(rhof_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!
end subroutine allocate_matter_emdrsrhof_mpt
