subroutine allocate_matter_rotlaw
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_matter,  only : omeg, omef, jomeg, jomef, jomeg_int, jomef_int
  use make_array_3d
  implicit none
!
  call alloc_array3d(omeg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(omef, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg_int, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef_int, 0, nrf, 0, ntf, 0, npf)
!
end subroutine allocate_matter_rotlaw
