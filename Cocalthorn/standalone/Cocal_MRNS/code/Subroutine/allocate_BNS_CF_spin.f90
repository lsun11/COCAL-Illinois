subroutine allocate_BNS_CF_spin
  use phys_constant, only : long
  use grid_parameter
  use def_matter
  use def_velocity_potential
  use def_velocity_rot
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(vep  , 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepxf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepyf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepzf, 0, nrf, 0, ntf, 0, npf)

  call alloc_array3d(vepxg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vepyg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vepzg, 0, nrg, 0, ntg, 0, npg)

  call alloc_array2d(alm, 0, nlg, 0, nlg)
!
  call alloc_array3d(wxspf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(wyspf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(wzspf, 0, nrf, 0, ntf, 0, npf)

  call alloc_array3d(wxspg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(wyspg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(wzspg, 0, nrg, 0, ntg, 0, npg)

  wxspf=0.0d0;  wyspf=0.0d0;  wzspf=0.0d0
  wxspg=0.0d0;  wyspg=0.0d0;  wzspg=0.0d0
!
end subroutine allocate_BNS_CF_spin
