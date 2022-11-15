subroutine copy_def_horizon_to_mpt(impt)
  use grid_parameter, only : ntg, npg
  use def_horizon
  use def_horizon_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array2dto3d_mpt(impt, ahz, ahz_, 0, ntg, 0, npg)
!
end subroutine copy_def_horizon_to_mpt
