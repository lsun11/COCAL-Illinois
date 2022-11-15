subroutine update_coordinates_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : rgin
  use def_bh_parameter, only : ome_bh
  use def_binary_parameter, only : dis
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_asymptotic_patch
  use grid_points_binary_in_asympto
  use weight_midpoint_binary_excision
  implicit none
  integer :: impt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_def_binary_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_bh_parameter_from_mpt(impt)
    call coordinate_patch_kit_grav_mpt
!!!    call calc_parameter_binary_excision
    if (impt.ne.nmpt) then
      call calc_grid_points_binary_excision
      call calc_grid_points_binary_in_asympto(impt,nmpt)
      call copy_grid_points_binary_in_asympto_to_mpt(impt)
    end if
    call calc_weight_midpoint_binary_excision
    call calc_vector_x_grav(0)
    call calc_vector_phi_grav(0)
    call calc_vector_bh(0)
!
    call copy_to_mpatch_all_BBH_CF(impt)
  end do
  call copy_from_mpatch_all_BBH_CF(nmpt)
! -- coordinates of asymptotic patch in central patches
  do impt = 1, 2
    call calc_grid_points_asymptotic_patch(impt,nmpt)
    call copy_grid_points_asymptotic_patch_to_mpt(impt)
  end do
!                     
end subroutine update_coordinates_mpt
