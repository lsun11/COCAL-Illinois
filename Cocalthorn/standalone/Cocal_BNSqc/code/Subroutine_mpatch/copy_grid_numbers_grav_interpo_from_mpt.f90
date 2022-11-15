subroutine copy_grid_numbers_grav_interpo_from_mpt(impt)
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use grid_parameter_mpt, only : grid_param_int_
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1; nrg_itp = grid_param_int_(i,impt)
  i=i+1; ntg_itp = grid_param_int_(i,impt)
  i=i+1; npg_itp = grid_param_int_(i,impt)
!
end subroutine copy_grid_numbers_grav_interpo_from_mpt
