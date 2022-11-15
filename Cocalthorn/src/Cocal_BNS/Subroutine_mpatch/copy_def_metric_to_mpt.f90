subroutine copy_def_metric_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use def_metric_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, psi, psi_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, alph, alph_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, alps, alps_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvxd, bvxd_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvyd, bvyd_, 0, nrg, 0, ntg, 0, npg)
  call copy_array3dto4d_mpt(impt, bvzd, bvzd_, 0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_def_metric_to_mpt
