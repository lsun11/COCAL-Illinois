subroutine allocate_matter_emdrsrhof
  use grid_parameter
  use def_matter, only : emd, rs, rhof
  use make_array_2d
  use make_array_3d
  implicit none
!
  call alloc_array2d(rs, 0, ntf, 0, npf)
  call alloc_array3d(emd, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(rhof, 0, nrf, 0, ntf, 0, npf)
!
end subroutine allocate_matter_emdrsrhof
