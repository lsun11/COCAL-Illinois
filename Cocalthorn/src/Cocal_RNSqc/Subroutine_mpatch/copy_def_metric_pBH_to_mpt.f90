subroutine copy_def_metric_pBH_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH
  use def_metric_pBH_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt,     wme,     wme_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, log_wme, log_wme_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, log_N  , log_N_  , 0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_def_metric_pBH_to_mpt
