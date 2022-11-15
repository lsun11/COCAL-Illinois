subroutine allocate_matter_emdrsrho
  use phys_constant, only : long
  use grid_parameter
  use def_matter, only : rs, emd, emdg, rhof, rhog
  use make_array_2d
  use make_array_3d
  implicit none
!
  call alloc_array2d(rs, 0, ntf, 0, npf)
  call alloc_array3d(emd, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(rhof, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(emdg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rhog, 0, nrg, 0, ntg, 0, npg)
!
end subroutine allocate_matter_emdrsrho
