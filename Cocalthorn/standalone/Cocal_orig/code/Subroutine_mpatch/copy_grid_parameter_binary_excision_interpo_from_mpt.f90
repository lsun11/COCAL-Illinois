subroutine copy_grid_parameter_binary_excision_interpo_from_mpt(impt)
  use grid_parameter_binary_excision_mpt
  use grid_parameter_binary_excision_interpo
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; ex_nrg_itp  = grid_param_bin_ex_int_(i,impt)
  i=i+1; ex_ndis_itp = grid_param_bin_ex_int_(i,impt)

  i=0
  i=i+1; ex_radius_itp = grid_param_bin_ex_real_(i,impt)
  i=i+1; ex_rgin_itp   = grid_param_bin_ex_real_(i,impt)
  i=i+1; ex_rgmid_itp  = grid_param_bin_ex_real_(i,impt)
  i=i+1; ex_rgout_itp  = grid_param_bin_ex_real_(i,impt)
!  
end subroutine copy_grid_parameter_binary_excision_interpo_from_mpt
