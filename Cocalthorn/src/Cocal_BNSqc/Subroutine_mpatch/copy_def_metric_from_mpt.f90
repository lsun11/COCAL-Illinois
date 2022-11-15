subroutine copy_def_metric_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use def_metric_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, psi_, psi, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, alph_, alph, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, alps_, alps, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvxd_, bvxd, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvyd_, bvyd, 0, nrg, 0, ntg, 0, npg)
  call copy_array4dto3d_mpt(impt, bvzd_, bvzd, 0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_def_metric_from_mpt
