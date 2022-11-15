subroutine allocate_matter
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_matter
  use def_velocity_rot
  use make_array_2d
  use make_array_3d
  implicit none
!
  call alloc_array2d(rs,        0, ntf, 0, npf)
  call alloc_array3d(emdg,      0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(emd,       0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(utg,       0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(utf,       0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(omeg,      0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(omef,      0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg,     0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg_int, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef_int, 0, nrf, 0, ntf, 0, npf)
!
  call alloc_array3d(vep,       0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepxf,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepyf,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepzf,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(wxspf,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(wyspf,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(wzspf,     0, nrf, 0, ntf, 0, npf)
!
end subroutine allocate_matter
