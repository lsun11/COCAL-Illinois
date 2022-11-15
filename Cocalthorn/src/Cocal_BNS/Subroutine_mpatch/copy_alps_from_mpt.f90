subroutine copy_alps_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : alps
  use def_metric_mpt, only : alps_
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, alps_, alps, 0, nrg, 0, ntg, 0, npg)
!
end subroutine copy_alps_from_mpt
