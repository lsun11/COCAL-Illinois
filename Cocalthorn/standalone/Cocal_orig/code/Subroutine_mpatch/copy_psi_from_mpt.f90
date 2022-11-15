subroutine copy_psi_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi
  use def_metric_mpt, only : psi_
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, psi_, psi, 0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_psi_from_mpt
