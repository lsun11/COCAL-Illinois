!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Include_file/include_subroutines_peos.f90'
!
include '../Module/def_matter_velocity.f90'
include '../Subroutine/interpo_lag4th_2Dsurf.f90'
!
include '../Analysis/Module/grid_parameter_cartesian.f90'
include '../Analysis/Module/coordinate_grav_xyz.f90'
include '../Analysis/Module/def_metric_cartesian.f90'
include '../Analysis/Module/def_matter_cartesian.f90'
include '../Analysis/Module/interface_modules_cartesian.f90'
!
include '../Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90'
include '../Analysis/Subroutine/coordinate_patch_cartesian.f90'
include '../Analysis/Subroutine/allocate_metric_and_matter_cartesian.f90'
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_test.f90'
include '../Analysis/Subroutine/IO_output_surface.f90'
include '../Analysis/Subroutine/interpo_gr2cgr_4th.f90'
include '../Analysis/Subroutine/interpo_fl2cgr_4th.f90'
!!!include '../Analysis/Subroutine/interpolation_metric.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian.f90'
include '../Analysis/Subroutine/interpolation_matter.f90'
include '../Analysis/Subroutine/calc_potential_minus_rinv_excision.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/IO_input_potential_test_3D.f90'
!include '../Analysis/Subroutine/interpolation_cartesian.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_wave_binary_test
!
  use def_metric, only : psi
  use def_metric_cartesian, only : psica
  use grid_parameter_binary_excision
  use interface_modules_cartesian
  implicit none
!
  call coordinate_patch_kit_grav_noGreen
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call coordinate_patch_cartesian
  call allocate_poisson_solver_test
  call allocate_metric_and_matter_cartesian
  call IO_input_potential_test_3D
  call calc_vector_x_grav(0)
!!!  call calc_potential_minus_rinv
  call calc_potential_minus_rinv_excision
!!!  call interpolation_metric(psi,psica)
  call interpolation_fillup_cartesian(psi,psica)
!
  call IO_output_cartesian_contour_potential_test
!
END PROGRAM interpolation_contour_wave_binary_test
