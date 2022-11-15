!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_grav_bhex_nb.f90'
include '../Module/radial_green_fn_grav_bhex_dd.f90'
include '../Module/radial_green_fn_grav_bhex_nd.f90'
include '../Module/radial_green_fn_grav_bhex_nh.f90'
include '../Module/radial_green_fn_grav_bhex_dh.f90'
include '../Module/copy_array_3d.f90'
include '../Module/gnufor2.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric_type0.f90'
include '../Module_interface/interface_error_metric_type1.f90'
include '../Module_interface/interface_error_metric_type2.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_nosym.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_3rd_nosym.f90'
include '../Module_interface/interface_grgrad_4th_gridpoint_bhex.f90'
include '../Module_interface/interface_grgrad_midpoint_r3rd_type0.f90'
include '../Module_interface/interface_grgrad_midpoint_r4th_type0.f90'
include '../Module_interface/interface_sourceterm_poisson_solver_test.f90'
include '../Module_interface/interface_sourceterm_surface_int_homosol.f90'
include '../Module_interface/interface_poisson_solver_binary_bhex_homosol.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_surface_int.f90'
include '../Module_interface/interface_poisson_solver_binary_bhex.f90'
include '../Module_interface/interface_poisson_solver_binary_surf_int.f90'
include '../Module_interface/interface_poisson_solver_binary_vol_int.f90'
include '../Module_interface/interface_poisson_solver_bhex_surf_int.f90'
include '../Module_interface/interface_poisson_solver_bhex_surf_int_all.f90'
include '../Module_interface/interface_copy_to_hgfn_and_gfnsf.f90'
include '../Module_interface/interface_interpo_gr2gr_4th.f90'
include '../Module_interface/interface_update_parameter.f90'
include '../Module_interface/interface_interpolation_fillup_binary.f90'
!
include '../Module_interface/interface_bh_boundary_d_psi.f90'
include '../Module_interface/interface_bh_boundary_d_alps.f90'
include '../Module_interface/interface_bh_boundary_d_bvxd.f90'
include '../Module_interface/interface_bh_boundary_d_bvyd.f90'
include '../Module_interface/interface_bh_boundary_d_bvzd.f90'
include '../Module_interface/interface_bh_boundary_d_Bfun.f90'
include '../Module_interface/interface_bh_boundary_n_Bfun.f90'
include '../Module_interface/interface_bh_boundary_d_potx.f90'
include '../Module_interface/interface_bh_boundary_d_poty.f90'
include '../Module_interface/interface_bh_boundary_d_potz.f90'
include '../Module_interface/interface_outer_boundary_d_psi.f90'
include '../Module_interface/interface_outer_boundary_d_alps.f90'
include '../Module_interface/interface_outer_boundary_d_bvxd.f90'
include '../Module_interface/interface_outer_boundary_d_bvyd.f90'
include '../Module_interface/interface_outer_boundary_d_bvzd.f90'
include '../Module_interface/interface_outer_boundary_d_Bfun.f90'
include '../Module_interface/interface_outer_boundary_n_Bfun.f90'
include '../Module_interface/interface_outer_boundary_d_potx.f90'
include '../Module_interface/interface_outer_boundary_d_poty.f90'
include '../Module_interface/interface_outer_boundary_d_potz.f90'
include '../Module_interface/interface_sourceterm_HaC_CF.f90'
include '../Module_interface/interface_sourceterm_trG_CF.f90'
include '../Module_interface/interface_sourceterm_MoC_CF.f90'
include '../Module_interface/interface_sourceterm_MoC_CF_type1_bhex.f90'
include '../Module_interface/interface_sourceterm_Bfun.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary_parity.f90'
include '../Module_interface/interface_interpolation_fillup_binary_parity.f90'
include '../Module_interface/interface_compute_shift_v2.f90'
include '../Module_interface/interface_compute_dBfun.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
include '../Subroutine/error_metric_type0.f90'
include '../Subroutine/error_metric_type1.f90'
include '../Subroutine/error_metric_type2.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/printout_error_all_metric.f90'
include '../Subroutine/coordinate_patch_kit_bhex.f90'
include '../Subroutine/copy_hgfn_nb_to_hgfn.f90'
include '../Subroutine/copy_hgfn_dd_to_hgfn.f90'
include '../Subroutine/copy_to_hgfn_and_gfnsf.f90'
include '../Subroutine/sourceterm_poisson_solver_test.f90'
include '../Subroutine/sourceterm_surface_int_homosol.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_surface_int.f90'
include '../Subroutine/poisson_solver_binary_bhex.f90'
include '../Subroutine/poisson_solver_binary_bhex_homosol.f90'
include '../Subroutine/poisson_solver_binary_vol_int.f90'
include '../Subroutine/poisson_solver_binary_surf_int.f90'
include '../Subroutine/poisson_solver_bhex_surf_int.f90'
include '../Subroutine/poisson_solver_bhex_surf_int_all.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/grdr_gridpoint_type0_nosym.f90'
include '../Subroutine/grdr_gridpoint_type0_3rd_nosym.f90'
include '../Subroutine/grgrad_4th_gridpoint_bhex.f90'
include '../Subroutine/grgrad_midpoint_r3rd_type0.f90'
include '../Subroutine/grgrad_midpoint_r4th_type0.f90'
include '../Subroutine/interpo_gr2gr_4th.f90'
include '../Subroutine/interpolation_fillup_binary.f90'
!
include '../Subroutine/allocate_BBH_CF.f90'
include '../Subroutine/iteration_BBH_CF_type1.f90'
include '../Subroutine/bh_boundary_d_psi.f90'
include '../Subroutine/bh_boundary_d_alps.f90'
include '../Subroutine/bh_boundary_d_bvxd.f90'
include '../Subroutine/bh_boundary_d_bvyd.f90'
include '../Subroutine/bh_boundary_d_bvzd.f90'
include '../Subroutine/bh_boundary_d_Bfun.f90'
include '../Subroutine/bh_boundary_n_Bfun.f90'
include '../Subroutine/bh_boundary_d_potx.f90'
include '../Subroutine/bh_boundary_d_poty.f90'
include '../Subroutine/bh_boundary_d_potz.f90'
include '../Subroutine/outer_boundary_d_psi.f90'
include '../Subroutine/outer_boundary_d_alps.f90'
include '../Subroutine/outer_boundary_d_bvxd.f90'
include '../Subroutine/outer_boundary_d_bvyd.f90'
include '../Subroutine/outer_boundary_d_bvzd.f90'
include '../Subroutine/outer_boundary_d_Bfun.f90'
include '../Subroutine/outer_boundary_n_Bfun.f90'
include '../Subroutine/outer_boundary_d_potx.f90'
include '../Subroutine/outer_boundary_d_poty.f90'
include '../Subroutine/outer_boundary_d_potz.f90'
include '../Subroutine/sourceterm_HaC_CF.f90'
include '../Subroutine/sourceterm_trG_CF.f90'
include '../Subroutine/sourceterm_MoC_CF.f90'
include '../Subroutine/sourceterm_MoC_CF_type1_bhex.f90'
include '../Subroutine/sourceterm_Bfun.f90'
include '../Subroutine/interpolation_fillup_binary_parity.f90'
include '../Subroutine/printout_error_metric_combined.f90'
include '../Subroutine/printout_error_metric_combined_B.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary_parity.f90'
include '../Subroutine/compute_shift_v2.f90'
include '../Subroutine/compute_dBfun.f90'
include '../Subroutine/reset_bh_boundary_all.f90'
include '../Subroutine/reset_metric_CF.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_BBH_CF_type1
!
  use grid_parameter, only : outdata_type, iter_max
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
  use radial_green_fn_grav_bhex_nd
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  use interface_copy_to_hgfn_and_gfnsf
  implicit none
  integer :: iseq, iter_count, total_iteration
!
  call coordinate_patch_kit_bhex
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_dd
  call calc_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
  call calc_hgfn_bhex_nd
! -- No boundary Green's fn
  call allocate_hgfn_bhex_nb
  call calc_hgfn_bhex_nb
!  call copy_to_hgfn_and_gfnsf(hgfn_nb,gfnsf_nb)
! --
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call allocate_weight_midpoint_binary_excision
  call calc_weight_midpoint_binary_excision
  call allocate_BBH_CF
!  call test_analytic_solution_bhex_psialph
! 
  call iteration_BBH_CF_type1(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!  call printout_NS_shape_seq(iseq)
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!!  call printout_debug
!  call printout_NS_shape
!
END PROGRAM Main_BBH_CF_type1
