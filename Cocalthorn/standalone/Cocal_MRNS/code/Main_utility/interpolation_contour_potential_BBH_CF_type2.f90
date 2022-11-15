!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/weight_midpoint_binary_excision.f90'
!
include '../Module/def_metric_excurve_grid.f90'
include '../Module/def_bh_parameter.f90'
include '../Module/def_vector_bh.f90'
include '../Module/def_vector_irg.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Function/lagint_2nd.f90'
!
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_surface_int.f90'
include '../Module_interface/interface_sourceterm_HaC_CF.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_nosym.f90'
include '../Module_interface/interface_grgrad_4th_gridpoint_bhex.f90'
include '../Module_interface/interface_surf_source_adm_mass.f90'
include '../Module_interface/interface_surf_source_komar_mass.f90'
include '../Module_interface/interface_source_ang_mom_inf.f90'
include '../Module_interface/interface_source_ang_mom_thr.f90'
include '../Module_interface/interface_source_ang_mom_smarr.f90'
include '../Module_interface/interface_surf_int_grav_rg.f90'
include '../Module_interface/interface_vol_int_grav_bhex.f90'
include '../Analysis/Module/grid_parameter_cartesian.f90'
include '../Analysis/Module/coordinate_grav_xyz.f90'
include '../Analysis/Module/def_metric_cartesian.f90'
include '../Analysis/Module/def_matter_cartesian.f90'
include '../Analysis/Module/interface_modules_cartesian.f90'
!
include '../Subroutine/interpo_lag4th_2Dsurf.f90'
!
include '../Include_file/include_subroutines_peos.f90'
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
include '../Analysis/Subroutine/read_omega.f90'

!include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/allocate_poisson_bbh_test.f90'
include '../Subroutine/allocate_BBH_CF.f90'
include '../Subroutine/allocate_BBH_CF_AH.f90'
include '../Subroutine/IO_input_potential_test_3D.f90'
include '../Subroutine/test_analytic_solution_bhex_psialph.f90'

include '../Subroutine/interpo_linear1p_type0_2Dsurf.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_surface_int.f90'
include '../Subroutine/sourceterm_HaC_CF.f90'
include '../Subroutine/grdr_gridpoint_type0_nosym.f90'
include '../Subroutine/grgrad_4th_gridpoint_bhex.f90'
include '../Subroutine/read_parameter_bh.f90'
include '../Subroutine/calc_vector_bh.f90'
include '../Subroutine/calc_vector_irg.f90'
include '../Subroutine/excurve_CF_gridpoint.f90'
include '../Subroutine/excurve_CF_gridpoint_bhex.f90'

include '../Subroutine/calc_physical_quantities_BBH_CF.f90'
include '../Subroutine/calc_mass_BBH_CF.f90'
include '../Subroutine/calc_mass_BBH_CF_adm_vol.f90'
include '../Subroutine/surf_source_adm_mass.f90'
include '../Subroutine/surf_source_komar_mass.f90'
include '../Subroutine/source_ang_mom_inf.f90'
include '../Subroutine/source_ang_mom_thr.f90'
include '../Subroutine/source_ang_mom_smarr.f90'
include '../Subroutine/printout_physq_console_BBH.f90'
include '../Subroutine/calc_mass_ir.f90'
!include '../Subroutine/calc_circular_orbit.f90'
include '../Subroutine/calc_ang_mom_BBH_CF_inf.f90'
include '../Subroutine/calc_ang_mom_BBH_CF_thr.f90'
include '../Subroutine/calc_ang_mom_BBH_CF_smarr.f90'
include '../Subroutine/save_solution.f90'
include '../Subroutine/write_omega_last.f90'
include '../Subroutine/test_excurve.f90'
include '../Subroutine/modify_r0_excurve.f90'
include '../Subroutine/calc_app_hor_area_BBH_CF.f90'
include '../Subroutine/surf_int_grav_rg.f90'
include '../Subroutine/vol_int_grav_bhex.f90'
!include '../Subroutine/.f90'
!include '../Subroutine/.f90'
!include '../Subroutine/.f90'
!include '../Subroutine/.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_BBH_CF_type1
!
  use def_metric
  use def_metric_cartesian 
  use interface_modules_cartesian
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_parameter
  use weight_midpoint_grav
  use weight_midpoint_binary_excision
  implicit none
  integer :: irg, itg, ipg
!
  call coordinate_patch_kit_grav_noGreen
  call allocate_weight_midpoint_grav
  call weight_calc_midpoint_grav
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call coordinate_patch_cartesian
!  call allocate_poisson_bbh_test
!  call allocate_BBH_CF
  call allocate_BBH_CF_AH
  call allocate_metric_and_matter_cartesian
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call read_parameter_bh
  call allocate_weight_midpoint_binary_excision
  call calc_weight_midpoint_binary_excision
  call IO_input_initial_3D
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
!!  call makeline
!
  call IO_output_cartesian_contour_potential_BBH_CF
!!  call IO_output_cartesian_contour_potential_BBH_CF_v1
!!  call IO_output_cartesian_planes
!  call test_analytic_solution_bhex_psialph
! writes plot_x, plot_y, plot_z 
!!  call IO_output_BBH_CF
!  call average_error
!
! The following are used for testing circular orbit parameters
  call excurve
  call excurve_CF_gridpoint_bhex
  call calc_mass_BBH_CF
  call calc_mass_BBH_CF_adm_vol
  call calc_ang_mom_BBH_CF_inf
!  call modify_r0_excurve
  call calc_ang_mom_BBH_CF_thr
  call calc_ang_mom_BBH_CF_smarr
  call calc_app_hor_area_BBH_CF
!  call test_excurve 
  call printout_physq_console_BBH
!
END PROGRAM interpolation_contour_potential_BBH_CF_type1
