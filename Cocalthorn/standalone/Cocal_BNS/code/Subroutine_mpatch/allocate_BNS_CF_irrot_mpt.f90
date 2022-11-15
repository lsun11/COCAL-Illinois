subroutine allocate_BNS_CF_irrot_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg, nlg
  use def_metric_mpt
  use def_matter_mpt
  use def_velocity_potential_mpt
!  use def_velocity_rot_mpt
  use make_array_3d
  use make_array_4d
  use make_array_5d
  implicit none
!
!  call alloc_array5d(vrot_, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, nmpt)

  call alloc_array4d(vep_  , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(vepxf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(vepyf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(vepzf_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)

  call alloc_array4d(vepxg_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(vepyg_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(vepzg_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)

!  call alloc_array3d(alm_, 0, nlg, 0, nlg, 1, nmpt)
!
end subroutine allocate_BNS_CF_irrot_mpt
