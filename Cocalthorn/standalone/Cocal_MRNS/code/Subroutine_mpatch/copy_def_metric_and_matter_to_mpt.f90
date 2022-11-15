subroutine copy_def_metric_and_matter_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric
  use def_metric_mpt
  use def_metric_excurve_grid
  use def_metric_excurve_grid_mpt
  use def_matter
  use def_matter_mpt
  use copy_array_3dto4d_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, psi , psi_ , 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, alph, alph_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, alps, alps_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvxd, bvxd_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvyd, bvyd_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvzd, bvzd_, 0, nrg, 0, ntg, 0, npg)
!
  call copy_array3dto4d_mpt(impt, trk_grid, trk_grid_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, trk,      trk_,      1, nrg, 1, ntg, 1, npg)
!
  call copy_array2dto3d_mpt(impt, rs       , rs_       , 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, emdg     , emdg_     , 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, emd      , emd_      , 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, utg      , utg_      , 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, utf      , utf_      , 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, omeg     , omeg_     , 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, omef     , omef_     , 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, jomeg    , jomeg_    , 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, jomef    , jomef_    , 0, nrf, 0, ntf, 0, npf)
  call copy_array3dto4d_mpt(impt, jomeg_int, jomeg_int_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, jomef_int, jomef_int_, 0, nrf, 0, ntf, 0, npf)

end subroutine copy_def_metric_and_matter_to_mpt
