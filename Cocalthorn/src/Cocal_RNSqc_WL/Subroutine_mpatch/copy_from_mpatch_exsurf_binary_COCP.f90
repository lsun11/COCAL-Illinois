subroutine copy_from_mpatch_exsurf_binary_COCP(impt)
!
  implicit none
  integer :: impt
!
  call copy_grid_parameter_from_mpt(impt)
  call copy_grid_parameter_binary_excision_from_mpt(impt)
  call copy_grid_points_binary_excision_from_mpt(impt)
  call copy_coordinate_grav_r_from_mpt(impt)
  call copy_coordinate_grav_theta_from_mpt(impt)
  call copy_coordinate_grav_phi_from_mpt(impt)
  call copy_coordinate_grav_extended_from_mpt(impt)
  call copy_def_binary_parameter_from_mpt(impt)
!
end subroutine copy_from_mpatch_exsurf_binary_COCP
