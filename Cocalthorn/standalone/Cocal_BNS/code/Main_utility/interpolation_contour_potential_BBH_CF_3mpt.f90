!______________________________________________
include '../Include_file/include_modulefiles_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_modulefiles_analysis_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_analysis_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_subroutines_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_subroutines_analysis_plot_BBH_CF_3mpt.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_BBH_CF_3mpt
!
  use grid_parameter_binary_excision
  use grid_points_asymptotic_patch
  use phys_constant
  use def_metric, only : psi
  use def_metric_cartesian, only : psica
  use interface_modules_cartesian
  use grid_points_binary_excision
  use trigonometry_grav_phi
  implicit none
  integer :: impt, impt_ex
!
! -- Read parameters
  call allocate_grid_parameter_mpt
  call allocate_grid_parameter_binary_excision_mpt
  call allocate_def_bh_parameter_mpt
  do impt = 1, nmpt
    call read_parameter_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call read_parameter_binary_excision_mpt(impt)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
    call read_parameter_bh_mpt(impt)
    call copy_def_bh_parameter_to_mpt(impt)
  end do
!
! -- Allocate arrays
  call set_allocate_size_mpt
!
  call allocate_grid_points_binary_excision
  call allocate_grid_points_asymptotic_patch
  call allocate_metric_CF
  call allocate_grid_points_binary_excision_mpt
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_metric_CF_mpt
  call allocate_trig_grav_mphi
  call coordinate_patch_cartesian
  call allocate_metric_and_matter_cartesian
  call allocate_mpatch_all_BBH_CF
  call allocate_coordinate_patch_kit_grav_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_bh_parameter_from_mpt(impt)
    call coordinate_patch_kit_grav_noGreen_mpt
    call calc_parameter_binary_excision
    call calc_grid_points_binary_excision
    call copy_to_mpatch_interpolation_utility(impt)
    call IO_input_converged_solution_3D_CF_BH_mpt(impt)
    call copy_def_metric_to_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call copy_def_bh_parameter_to_mpt(impt)
    call copy_def_binary_parameter_to_mpt(impt)
  end do
  if (nmpt.gt.2) then
    call copy_from_mpatch_interpolation_utility(nmpt)
    do impt = 1, nmpt
! -- coordinates of asymptotic patch in central patches
      call calc_grid_points_asymptotic_patch(impt,nmpt)
      call copy_grid_points_asymptotic_patch_to_mpt(impt)
    end do
  end if
!
  do impt = 1, nmpt
    call copy_from_mpatch_interpolation_utility(impt)
    call copy_def_metric_from_mpt(impt)
    call IO_output_plot_xyz_BBH_CF_mpt(impt)
  end do
!
END PROGRAM interpolation_contour_potential_BBH_CF_3mpt
