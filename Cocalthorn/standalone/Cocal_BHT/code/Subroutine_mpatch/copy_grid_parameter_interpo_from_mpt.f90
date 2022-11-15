subroutine copy_grid_parameter_interpo_from_mpt(impt)
  use grid_parameter_interpo
  use grid_parameter_mpt, only : grid_param_int_, grid_param_real_
  implicit none
  integer :: i, impt
!
  i=0
  i=i+1; nrg_itp = grid_param_int_(i,impt)
  i=i+1; ntg_itp = grid_param_int_(i,impt)
  i=i+1; npg_itp = grid_param_int_(i,impt)
!
  i=0
  i=i+1; rgin_itp  = grid_param_real_(i,impt)
  i=i+1; rgmid_itp = grid_param_real_(i,impt)
  i=i+1; rgout_itp = grid_param_real_(i,impt)
  i=i+1; ratio_itp = grid_param_real_(i,impt)
!
end subroutine copy_grid_parameter_interpo_from_mpt
