!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/grid_points_asymptotic_patch.f90'
include '../Module/make_char1_array_2d.f90'
include '../Module/make_char2_array_2d.f90'
include '../Module/make_int_array_3d.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module_mpatch/grid_points_asymptotic_patch_mpt.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Include_file/include_subroutines_peos.f90'
include '../Include_file/include_modulefiles_mpatch.f90'
!
include '../Analysis/Module/grid_parameter_cartesian.f90'
include '../Analysis/Module/coordinate_grav_xyz.f90'
include '../Analysis/Module/def_metric_cartesian.f90'
include '../Analysis/Module/def_matter_cartesian.f90'
include '../Analysis/Module/interface_modules_cartesian.f90'
!
include '../Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90'
include '../Analysis/Subroutine/coordinate_patch_kit_grav_noGreen_mpt.f90'
include '../Analysis/Subroutine/coordinate_patch_cartesian.f90'
include '../Analysis/Subroutine/allocate_metric_and_matter_cartesian.f90'
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_test.f90'
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_test_mpt.f90'
include '../Analysis/Subroutine/interpo_gr2cgr_4th.f90'
include '../Analysis/Subroutine/interpo_fl2cgr_4th.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian_mpt.f90'
include '../Analysis/Subroutine/interpolation_matter.f90'
include '../Analysis/Subroutine/IO_output_plot_xyz.f90'
include '../Analysis/Subroutine/IO_output_plot_xyz_mpt.f90'
include '../Analysis/Subroutine/IO_output_plot_averaged_error_mpt.f90'
include '../Analysis/Subroutine/copy_to_mpatch_interpolation_utility.f90'
include '../Analysis/Subroutine/copy_from_mpatch_interpolation_utility.f90'
!
include '../Subroutine/allocate_metric_and_matter_BHNS_test.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/interpo_lag4th_2Dsurf.f90'
include '../Subroutine/IO_input_potential_test_3D.f90'
include '../Subroutine_mpatch/IO_input_potential_test_3D_mpt.f90'
include '../Subroutine_mpatch/test_analytic_BHNS_solution_mpt.f90'
include '../Subroutine_mpatch/test_analytic_bns_solution_mpt.f90'
include '../Subroutine_mpatch/test_analytic_solution_bhex_mpt.f90'
include '../Subroutine_mpatch/allocate_grid_points_asymptotic_patch_mpt.f90'
include '../Subroutine_mpatch/allocate_metric_and_matter_BHNS_test_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_asymptotic_patch_to_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_asymptotic_patch_from_mpt.f90'
include '../Include_file/include_subroutines_mpatch.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_binary_test_mpt
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
  do impt = 1, nmpt
    call read_parameter_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call read_parameter_binary_excision_mpt(impt)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
  end do
!
! -- Allocate arrays
  call set_allocate_size_mpt
!
  call allocate_grid_points_binary_excision
  call allocate_grid_points_asymptotic_patch
  call allocate_metric_and_matter_BHNS_test
  call allocate_grid_points_binary_excision_mpt
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_metric_and_matter_BHNS_test_mpt
  call allocate_trig_grav_mphi
  call coordinate_patch_cartesian
  call allocate_metric_and_matter_cartesian
  call allocate_mpatch_all_test
  call allocate_coordinate_patch_kit_grav_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call coordinate_patch_kit_grav_noGreen_mpt
    call calc_parameter_binary_excision
    call calc_grid_points_binary_excision
    call copy_to_mpatch_interpolation_utility(impt)
    call IO_input_potential_test_3D_mpt(impt)
    call copy_poisson_solver_test_to_mpt(impt)
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
    call copy_poisson_solver_test_from_mpt(impt)
!    call test_analytic_BHNS_solution_mpt(impt)
    call test_analytic_solution_bhex_mpt(impt)
    call IO_output_plot_xyz_mpt(impt)
    call IO_output_plot_averaged_error_mpt(impt)
    if(impt.eq.1) impt_ex = 2
    if(impt.eq.2) impt_ex = 1
    call interpolation_fillup_cartesian_mpt(psi, psica, impt, impt_ex)
    call IO_output_cartesian_contour_potential_test_mpt(impt)
  end do
!
END PROGRAM interpolation_contour_potential_binary_test_mpt
