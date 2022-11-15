!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Include_file/include_subroutines_peos.f90'
!
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
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_2pot_test.f90'
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_BBH_CF.f90'
include '../Analysis/Subroutine/IO_output_cartesian_contour_potential_BBH_CF_v1.f90'
include '../Analysis/Subroutine/IO_output_cartesian_planes.f90'
include '../Analysis/Subroutine/interpo_gr2cgr_4th.f90'
include '../Analysis/Subroutine/interpo_fl2cgr_4th.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian_parity.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian_bh_parity.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian_bh.f90'
include '../Analysis/Subroutine/interpolation_fillup_cartesian_bh_all.f90'
include '../Analysis/Subroutine/interpolation_matter.f90'
include '../Analysis/Subroutine/IO_output_bhex_2pot.f90'
include '../Analysis/Subroutine/IO_output_BBH_CF.f90'
include '../Analysis/Subroutine/average_error.f90'
include '../Analysis/Subroutine/makeline.f90'
!include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/allocate_poisson_bbh_test.f90'
include '../Subroutine/IO_input_potential_test_3D.f90'
include '../Subroutine/test_analytic_solution_bhex_psialph.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_BBH_CF
!
  use def_metric
  use def_metric_cartesian  
  use interface_modules_cartesian
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_parameter
  implicit none
  integer :: irg, itg, ipg
!
  call coordinate_patch_kit_grav_noGreen
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call coordinate_patch_cartesian
  call allocate_poisson_bbh_test
  call allocate_metric_and_matter_cartesian
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call IO_input_potential_test_3D
  call interpolation_fillup_cartesian_bh_all
!  call interpolation_fillup_cartesian(psi,psica)
!  call interpolation_fillup_cartesian(alph,alphca)
!  call interpolation_fillup_cartesian_parity(bvxd,bvxdca,-1.0d0)
!  call interpolation_fillup_cartesian_parity(bvyd,bvydca,-1.0d0)
!  call interpolation_fillup_cartesian_parity(bvzd,bvzdca,+1.0d0)
!  call interpolation_fillup_cartesian_bh(psi,psica)
!  call interpolation_fillup_cartesian_bh(alph,alphca)
!  call interpolation_fillup_cartesian_bh_parity(bvxd,bvxdca,-1.0d0)
!  call interpolation_fillup_cartesian_bh_parity(bvyd,bvydca,-1.0d0)
!  call interpolation_fillup_cartesian_bh_parity(bvzd,bvzdca,+1.0d0)
  call makeline
!
  call IO_output_cartesian_contour_potential_BBH_CF
  call IO_output_cartesian_contour_potential_BBH_CF_v1
  call IO_output_cartesian_planes
!  call test_analytic_solution_bhex_psialph
  call IO_output_BBH_CF
!  call average_error
!
END PROGRAM interpolation_contour_potential_BBH_CF
