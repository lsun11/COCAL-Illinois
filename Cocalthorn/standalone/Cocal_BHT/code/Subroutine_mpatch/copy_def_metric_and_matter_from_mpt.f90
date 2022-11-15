subroutine copy_def_metric_and_matter_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_metric_excurve_grid
  use def_metric_excurve_grid_mpt
  use def_matter
  use def_matter_mpt
  use copy_array_4dto3d_mpt
  use copy_array_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, psi_ , psi , 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, alph_, alph, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, alps_, alps, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvxd_, bvxd, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvyd_, bvyd, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvzd_, bvzd, 0, nrg, 0, ntg, 0, npg)
!
  call copy_array4dto3d_mpt(impt, trk_grid_, trk_grid, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, trk_,      trk,      1, nrg, 1, ntg, 1, npg)
!
  call copy_array3dto2d_mpt(impt, rs_       , rs       , 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, emdg_     , emdg     , 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, emd_      , emd      , 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, utg_      , utg      , 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, utf_      , utf      , 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, omeg_     , omeg     , 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, omef_     , omef     , 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, jomeg_    , jomeg    , 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, jomef_    , jomef    , 0, nrf, 0, ntf, 0, npf)
  call copy_array4dto3d_mpt(impt, jomeg_int_, jomeg_int, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, jomef_int_, jomef_int, 0, nrf, 0, ntf, 0, npf)

end subroutine copy_def_metric_and_matter_from_mpt
