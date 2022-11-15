!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/grid_points_asymptotic_patch.f90'
include '../Module/grid_points_binary_in_asympto.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_grav_bhex_nb.f90'
include '../Module/radial_green_fn_grav_bhex_di.f90'
include '../Module/radial_green_fn_grav_bhex_dd.f90'
include '../Module/radial_green_fn_grav_bhex_dh.f90'
include '../Module/radial_green_fn_grav_bhex_nd.f90'
include '../Module/radial_green_fn_grav_bhex_nh.f90'
include '../Module/copy_array_3d.f90'
include '../Module/make_char2_array_2d.f90'
include '../Module/make_char1_array_2d.f90'
include '../Module/make_int_array_3d.f90'
include '../Module/grid_parameter_interpo.f90'
include '../Module/grid_parameter_binary_excision_interpo.f90'
include '../Module/coordinate_grav_extended_interpo.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric.f90'
include '../Module_interface/interface_error_metric_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_nosym.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_3rd_nosym.f90'
include '../Module_interface/interface_grgrad_midpoint_r3rd_type0.f90'
include '../Module_interface/interface_grgrad_midpoint_r4th_type0.f90'
include '../Module_interface/interface_sourceterm_poisson_solver_test.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_surface_int.f90'
include '../Module_interface/interface_sourceterm_surface_int_homosol.f90'
include '../Module_interface/interface_sourceterm_outsurf_interpo_from_asympto_mpt.f90' 
include '../Module_interface/interface_sourceterm_insurf_asympto_interpo_from_mpt.f90'
include '../Module_interface/interface_sourceterm_volume_int_bbh_2pot_test.f90'
include '../Module_interface/interface_sourceterm_2pot_test.f90'
include '../Module_interface/interface_poisson_solver_binary.f90'
include '../Module_interface/interface_poisson_solver_binary_bhex_homosol.f90'
include '../Module_interface/interface_poisson_solver_binary_bhex.f90'
include '../Module_interface/interface_poisson_solver_binary_surf_int.f90'
include '../Module_interface/interface_poisson_solver_binary_vol_int.f90'
include '../Module_interface/interface_poisson_solver_bhex_surf_int.f90'
include '../Module_interface/interface_poisson_solver_bhex_surf_int_all.f90'
include '../Module_interface/interface_poisson_solver_asymptotic_patch_homosol.f90'
include '../Module_interface/interface_copy_to_hgfn_and_gfnsf.f90'
include '../Module_interface/interface_bh_boundary_test.f90'
include '../Module_interface/interface_bh_boundary_test_mpt.f90'
include '../Module_interface/interface_bh_boundary_nh_psi_test.f90'
include '../Module_interface/interface_bh_boundary_nh_alph_test.f90'
include '../Module_interface/interface_bh_boundary_dh_psi_test.f90'
include '../Module_interface/interface_bh_boundary_dh_alph_test.f90'
include '../Module_interface/interface_bh_boundary_psi_test_mpt.f90'
include '../Module_interface/interface_bh_boundary_alph_test_mpt.f90'
include '../Module_interface/interface_bh_boundary_psi_test2_mpt.f90'
include '../Module_interface/interface_bh_boundary_alph_test2_mpt.f90'
include '../Module_interface/interface_interpo_gr2gr_4th_mpt.f90'
include '../Module_interface/interface_interpolation_fillup_binary_mpt.f90'
include '../Module_interface/interface_interpo_patch_to_active_patch_mpt.f90'
include '../Include_file/include_modulefiles_mpatch.f90'
include '../Module_mpatch/grid_points_asymptotic_patch_mpt.f90'
include '../Module_mpatch/grid_points_binary_in_asympto_mpt.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
!include '../Subroutine/iteration_poisson_bbh_test.f90'
include '../Subroutine/error_metric.f90'
include '../Subroutine/error_metric_type0.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/coordinate_patch_kit_bhex.f90'
include '../Subroutine/copy_hgfn_nb_to_hgfn.f90'
include '../Subroutine/copy_hgfn_di_to_hgfn.f90'
include '../Subroutine/copy_hgfn_dd_to_hgfn.f90'
include '../Subroutine/copy_hgfn_dh_to_hgfn.f90'
include '../Subroutine/copy_to_hgfn_and_gfnsf.f90'
include '../Subroutine/sourceterm_poisson_solver_test.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_surface_int.f90'
include '../Subroutine/sourceterm_surface_int_homosol.f90'
include '../Subroutine/sourceterm_volume_int_bbh_2pot_test.f90'
include '../Subroutine/sourceterm_2pot_test.f90'
include '../Subroutine/poisson_solver_binary.f90'
include '../Subroutine/poisson_solver_binary_bhex.f90'
include '../Subroutine/poisson_solver_binary_bhex_homosol.f90'
include '../Subroutine/poisson_solver_binary_vol_int.f90'
include '../Subroutine/poisson_solver_binary_surf_int.f90'
include '../Subroutine/poisson_solver_bhex_surf_int.f90'
include '../Subroutine/poisson_solver_bhex_surf_int_all.f90'
include '../Subroutine/poisson_solver_asymptotic_patch_homosol.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/grdr_gridpoint_type0_nosym.f90'
include '../Subroutine/grdr_gridpoint_type0_3rd_nosym.f90'
include '../Subroutine/grgrad_midpoint_r3rd_type0.f90'
include '../Subroutine/grgrad_midpoint_r4th_type0.f90'
include '../Subroutine/reset_bh_boundary.f90'
include '../Subroutine/bh_boundary_test.f90'
include '../Subroutine/bh_boundary_nh_psi_test.f90'
include '../Subroutine/bh_boundary_nh_alph_test.f90'
include '../Subroutine/bh_boundary_dh_psi_test.f90'
include '../Subroutine/bh_boundary_dh_alph_test.f90'
include '../Subroutine/test_source.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/test_analytic_solution.f90'
include '../Subroutine/allocate_metric_and_matter_BHNS_test.f90'
!
include '../Subroutine_mpatch/test_source_mpt.f90'
include '../Subroutine_mpatch/IO_output_poisson_test_3D_mpt.f90'
!include '../Subroutine_mpatch/iteration_poisson_bbh_2pot_test_mpt.f90'
include '../Subroutine_mpatch/iteration_poisson_bbh_2pot_test_3mpt.f90'
include '../Subroutine_mpatch/allocate_metric_and_matter_BHNS_test_mpt.f90'
include '../Subroutine_mpatch/allocate_grid_points_asymptotic_patch_mpt.f90'
include '../Subroutine_mpatch/allocate_grid_points_binary_in_asympto_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_asymptotic_patch_to_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_asymptotic_patch_from_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_binary_in_asympto_to_mpt.f90'
include '../Subroutine_mpatch/copy_grid_points_binary_in_asympto_from_mpt.f90'
include '../Subroutine_mpatch/copy_metric_and_matter_BHNS_test_from_mpt.f90'
include '../Subroutine_mpatch/copy_metric_and_matter_BHNS_test_to_mpt.f90'
include '../Subroutine_mpatch/bh_boundary_test_mpt.f90'
include '../Subroutine_mpatch/bh_boundary_psi_test_mpt.f90'
include '../Subroutine_mpatch/bh_boundary_alph_test_mpt.f90'
include '../Subroutine_mpatch/bh_boundary_psi_test2_mpt.f90'
include '../Subroutine_mpatch/bh_boundary_alph_test2_mpt.f90'
include '../Subroutine_mpatch/test_analytic_solution_bhex_mpt.f90'
include '../Subroutine_mpatch/test2_analytic_solution_bhex_mpt.f90'
include '../Subroutine_mpatch/interpolation_fillup_binary_mpt.f90'
include '../Subroutine_mpatch/interpo_gr2gr_4th_mpt.f90'
include '../Subroutine_mpatch/interpo_patch_to_active_patch_mpt.f90'
include '../Subroutine_mpatch/copy_coordinate_grav_extended_interpo_from_mpt.f90'
include '../Subroutine_mpatch/copy_grid_numbers_grav_from_mpt.f90'
include '../Subroutine_mpatch/copy_grid_numbers_grav_interpo_from_mpt.f90'
include '../Subroutine_mpatch/copy_grid_parameter_interpo_from_mpt.f90'
include '../Subroutine_mpatch/copy_grid_parameter_binary_excision_interpo_from_mpt.f90'
!include '../Module_mpatch/grid_parameter_mpt.f90'
include '../Include_file/include_subroutines_mpatch.f90'
include '../Subroutine_mpatch/sourceterm_outsurf_interpo_from_asympto_mpt.f90' 
include '../Subroutine_mpatch/sourceterm_insurf_asympto_interpo_from_mpt.f90'
!include '../Module/grid_points_binary_excision.f90'
!include '../Subroutine_mpatch/read_parameter_mpt.f90'
!include '../Subroutine_mpatch/read_parameter_binary_excision_mpt.f90'
!include '../Subroutine_mpatch/allocate_grid_parameter_mpt.f90'
!include '../Subroutine_mpatch/copy_grid_parameter_to_mpt.f90'
!include '../Subroutine_mpatch/copy_grid_parameter_from_mpt.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_poisson_binary_2pot_test2_3mpt
!
  use phys_constant, only : nmpt
  use grid_parameter, only : indata_type, outdata_type, iter_max
!  use grid_parameter_binary_excision, only : ex_nrg, ex_ndis
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_asymptotic_patch
  use grid_points_binary_in_asympto
  use weight_midpoint_binary_excision
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
  use radial_green_fn_grav_bhex_nd
  use radial_green_fn_grav_bhex_dh
  use radial_green_fn_grav_bhex_nh
  use interface_copy_to_hgfn_and_gfnsf
  use def_vector_x
  use def_vector_phi
  implicit none
  integer :: impt, itg
  integer :: iseq, iter_count, total_iteration
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
  call allocate_coordinate_patch_kit_grav_mpt
  call allocate_grid_points_binary_excision
  call allocate_grid_points_asymptotic_patch
  call allocate_grid_points_binary_in_asympto
  call allocate_weight_midpoint_binary_excision
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_nb
  call allocate_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
  call allocate_metric_and_matter_BHNS_test
!
  call allocate_mpatch_all_test
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_grid_points_binary_in_asympto_mpt
  call allocate_metric_and_matter_BHNS_test_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call coordinate_patch_kit_grav_mpt
    call calc_parameter_binary_excision
    if (impt.ne.nmpt) then 
      call calc_grid_points_binary_excision
      call calc_grid_points_binary_in_asympto(impt,nmpt)
      call copy_grid_points_binary_in_asympto_to_mpt(impt)
    end if
    call calc_weight_midpoint_binary_excision
    call calc_weight_midpoint_binary_excision_hybrid
    call calc_vector_x_grav(0)
    call calc_vector_x_matter(0)
    call calc_vector_phi_grav(0)
    call calc_vector_phi_matter(0)
!!! unnecessary for BBH    call test_source_mpt(impt)
    call copy_to_mpatch_all_test(impt)
    call copy_metric_and_matter_BHNS_test_to_mpt(impt)
  end do
  call copy_from_mpatch_all_test(nmpt)
! -- coordinates of asymptotic patch in central patches
  do impt = 1, 2
    call calc_grid_points_asymptotic_patch(impt,nmpt)
    call copy_grid_points_asymptotic_patch_to_mpt(impt)
  end do
! -- Analytic solution for 2 potential problem on all patches.
  do impt = 1, nmpt
    call copy_from_mpatch_all_test(impt)
    call copy_metric_and_matter_BHNS_test_from_mpt(impt)
    call test_analytic_solution_bhex_mpt(impt)
    call copy_to_mpatch_all_test(impt)
    call copy_metric_and_matter_BHNS_test_to_mpt(impt)
  end do
!
  call iteration_poisson_bbh_2pot_test_3mpt(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!  call printout_NS_shape_seq(iseq)
  do impt = 1, nmpt
    call copy_from_mpatch_all_test(impt)
    call copy_metric_and_matter_BHNS_test_from_mpt(impt)
    if (outdata_type.eq.'3D') call IO_output_poisson_test_3D_mpt(impt)
  end do
!
END PROGRAM Main_poisson_binary_2pot_test2_3mpt
