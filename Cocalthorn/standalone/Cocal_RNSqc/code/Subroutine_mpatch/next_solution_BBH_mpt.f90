subroutine next_solution_BBH_mpt
  use phys_constant,  only : long
  use grid_parameter, only : rgin
  use def_bh_parameter, only : rgin_deform
  implicit none
  integer :: impt
!
  do impt = 1, 2
    call copy_grid_parameter_from_mpt(impt)
    call copy_def_bh_parameter_from_mpt(impt)
    rgin = rgin + rgin_deform
    call copy_grid_parameter_to_mpt(impt)
    call copy_def_bh_parameter_to_mpt(impt)
  end do
!
end subroutine next_solution_BBH_mpt
