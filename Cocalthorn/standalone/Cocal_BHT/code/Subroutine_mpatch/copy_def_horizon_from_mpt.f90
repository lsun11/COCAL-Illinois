subroutine copy_def_horizon_from_mpt(impt)
  use grid_parameter, only : ntg, npg
  use def_horizon
  use def_horizon_mpt
  use copy_array_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto2d_mpt(impt, ahz_, ahz, 0, ntg, 0, npg)
!
end subroutine copy_def_horizon_from_mpt
