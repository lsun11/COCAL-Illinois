subroutine allocate_matter_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_metric
  use def_metric_mpt
  use def_matter
  use def_matter_mpt
  use def_velocity_rot_mpt
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(rs_       , 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emd_      , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(utf_      , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(omef_     , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(jomef_    , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(jomef_int_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)

  call alloc_array4d(emdg_     , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(utg_      , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(omeg_     , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(jomeg_    , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(jomeg_int_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
!
  call alloc_array4d(vep_      , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(wxspf_    , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(wyspf_    , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(wzspf_    , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!
end subroutine allocate_matter_mpt
