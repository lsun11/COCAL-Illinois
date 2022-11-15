subroutine allocate_BNS_CF_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_metric_mpt
  use def_matter_mpt
  use def_metric_excurve_grid_mpt
  use make_array_3d
  use make_array_4d
  use make_array_6d
  implicit none
!
  call alloc_array3d(rs_       , 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(emdg_     , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(emd_      , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(utg_      , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(utf_      , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(omeg_     , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(omef_     , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(jomeg_    , 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(jomef_    , 0, nrf, 0, ntf, 0, npf, 1, nmpt)
  call alloc_array4d(jomeg_int_, 0, nrg, 0, ntg, 0, npg, 1, nmpt)
  call alloc_array4d(jomef_int_, 0, nrf, 0, ntf, 0, npf, 1, nmpt)
!    
  call alloc_array4d(psi_ , 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alph_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alps_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvxd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvyd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvzd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!
  call alloc_array4d(trk_, 1,nrg, 1,ntg, 1,npg, 1,nmpt)
!!!  call alloc_array4d(tfkijkij_, 1,nrg, 1,ntg, 1,npg, 1,nmpt)
!!!  call alloc_array6d(tfkij_,    1,nrg, 1,ntg, 1,npg, 1,3, 1,3, 1,nmpt)
  call alloc_array4d(trk_grid_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!!!  call alloc_array4d(tfkijkij_grid_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!!!  call alloc_array6d(tfkij_grid_, 0,nrg, 0,ntg, 0,npg, 1,3, 1,3, 1,nmpt)
!
end subroutine allocate_BNS_CF_mpt
