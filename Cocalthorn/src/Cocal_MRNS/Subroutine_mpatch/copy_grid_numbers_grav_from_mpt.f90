subroutine copy_grid_numbers_grav_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_parameter_mpt, only : grid_param_int_
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; nrg = grid_param_int_(i,impt)
  i=i+1; ntg = grid_param_int_(i,impt)
  i=i+1; npg = grid_param_int_(i,impt)
!
end subroutine copy_grid_numbers_grav_from_mpt
