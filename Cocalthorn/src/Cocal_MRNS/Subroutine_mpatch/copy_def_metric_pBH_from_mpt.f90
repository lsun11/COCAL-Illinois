subroutine copy_def_metric_pBH_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH
  use def_metric_pBH_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt,     wme_,     wme, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, log_wme_, log_wme, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, log_N_,   log_N,   0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_def_metric_pBH_from_mpt
