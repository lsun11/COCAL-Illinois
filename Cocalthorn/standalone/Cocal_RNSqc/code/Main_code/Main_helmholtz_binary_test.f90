!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/grid_points_asymptotic_patch.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_helmholtz.f90'
include '../Module/radial_green_fn_hrethadv.f90'
include '../Module/radial_green_fn_hrethadv_homosol.f90'
include '../Module/copy_array_3d.f90'
include '../Module/copy_array_4d.f90'
include '../Module/make_char2_array_2d.f90'
include '../Module/make_char1_array_2d.f90'
include '../Module/make_int_array_3d.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric.f90'
include '../Module_interface/interface_error_metric_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_nosym.f90'
include '../Module_interface/interface_sourceterm_helmholtz_solver_test.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_outsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_surface_int.f90'
include '../Module_interface/interface_sourceterm_surface_int_homosol.f90'
include '../Module_interface/interface_sourceterm_insurf_asymptotic_patch.f90'
include '../Module_interface/interface_helmholtz_solver_vol_int.f90'
include '../Module_interface/interface_helmholtz_solver_surf_int.f90'
include '../Module_interface/interface_helmholtz_solver_binary.f90'
include '../Module_interface/interface_helmholtz_solver_binary_vol_int.f90'
include '../Module_interface/interface_helmholtz_solver_binary_surf_int.f90'
include '../Module_interface/interface_helmholtz_solver_outer_surf_int.f90'
include '../Module_interface/interface_helmholtz_solver_asymptotic_patch_homosol.f90'
include '../Module_interface/interface_copy_to_bsjy_and_sbsjy.f90'
include '../Module_interface/interface_interpo_binary_to_asymptotic_patch.f90'
!
include '../Include_file/include_modulefiles_mpatch.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
include '../Subroutine/iteration_helmholtz_solver_binary_test.f90'
include '../Subroutine/error_metric.f90'
include '../Subroutine/error_metric_type0.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/sourceterm_helmholtz_solver_test.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_outsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_surface_int.f90'
include '../Subroutine/sourceterm_surface_int_homosol.f90'
include '../Subroutine/helmholtz_solver_vol_int.f90'
include '../Subroutine/helmholtz_solver_surf_int.f90'
include '../Subroutine/helmholtz_solver_binary.f90'
include '../Subroutine/helmholtz_solver_binary_vol_int.f90'
include '../Subroutine/helmholtz_solver_binary_surf_int.f90'
include '../Subroutine/helmholtz_solver_outer_surf_int.f90'
include '../Subroutine/helmholtz_solver_asymptotic_patch_homosol.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/grdr_gridpoint_type0_nosym.f90'
include '../Subroutine/test_source_helical_binary.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/allocate_metric_and_matter_BHNS_test.f90'
include '../Subroutine/test_analytic_solution.f90'
include '../Subroutine/calc_radial_green_fn_hrethadv.f90'
include '../Subroutine/calc_radial_green_fn_hrethadv_homosol.f90'
include '../Subroutine/interpo_binary_to_asymptotic_patch.f90'
include '../Subroutine/copy_to_bsjy_and_sbsjy.f90'
include '../Subroutine/sphbess.f90'
include '../Subroutine/sphbess_and_dx.f90'
include '../Subroutine/sphbess_dx.f90'
include '../Subroutine/bessjy.f90'
include '../Subroutine/beschb.f90'
include '../Include_file/include_subroutines_mpatch.f90'
include '../Subroutine_mpatch/IO_output_poisson_test_3D_mpt.f90'
include '../Subroutine_mpatch/allocate_metric_and_matter_BHNS_test_mpt.f90'
include '../Subroutine_mpatch/copy_metric_and_matter_BHNS_test_to_mpt.f90'
include '../Subroutine_mpatch/copy_metric_and_matter_BHNS_test_from_mpt.f90'
include '../Subroutine_mpatch/copy_to_mpatch_helmholtz_test.f90'
include '../Subroutine_mpatch/copy_from_mpatch_helmholtz_test.f90'
include '../Subroutine_mpatch/sourceterm_insurf_asymptotic_patch.f90'
include '../Function/chebev.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_helmholtz_binary_test
!
  use phys_constant, only : nmpt
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_asymptotic_patch
  use weight_midpoint_grav
  use weight_midpoint_fluid
  use weight_midpoint_binary_excision
  use def_vector_x
  use def_vector_phi
  use radial_green_fn_helmholtz
  use radial_green_fn_hrethadv
  use radial_green_fn_hrethadv_homosol
  implicit none
  integer :: iseq, iter_count, total_iteration, impt, impt_bin
!
! -- Read parameters
  call allocate_grid_parameter_mpt
  call allocate_grid_parameter_binary_excision_mpt
  call read_parameter_mpt(1)
  call copy_grid_parameter_to_mpt(1)
  call read_parameter_binary_excision_mpt(1)
  call copy_grid_parameter_binary_excision_to_mpt(1)
!  
  call read_parameter_mpt(2)
  call copy_grid_parameter_to_mpt(2)
  call read_parameter_binary_excision_mpt(2)
  call copy_grid_parameter_binary_excision_to_mpt(2)
!  
! -- Allocate arrays
  call set_allocate_size_mpt
!
  call allocate_coordinate_patch_kit_grav_mpt
  call allocate_grid_points_binary_excision
  call allocate_weight_midpoint_grav
  call allocate_weight_midpoint_fluid
  call allocate_weight_midpoint_binary_excision
  call allocate_vector_x
  call allocate_vector_phi
  call allocate_radial_green_fn_helmholtz
  call allocate_radial_green_fn_hrethadv
  call allocate_radial_green_fn_hrethadv_homosol
  call allocate_metric_and_matter_BHNS_test
  call allocate_metric_and_matter_BHNS_test_mpt
!
  call allocate_mpatch_all_test
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call coordinate_patch_kit_grav_mpt
    if (impt.eq.1) then
      call calc_parameter_binary_excision
      call calc_grid_points_binary_excision
      call calc_weight_midpoint_binary_excision
      call test_source_helical_binary
    end if
!
    if (impt.eq.2) then
      impt_bin = 1
      call allocate_grid_points_asymptotic_patch
      call calc_parameter_binary_excision
      call calc_grid_points_asymptotic_patch(impt_bin,impt)
    end if
    call calc_vector_x_grav(0)
    call calc_vector_x_matter(0)
    call calc_vector_phi_grav(0)
    call calc_vector_phi_matter(0)
    call copy_to_mpatch_helmholtz_test(impt)
    call copy_metric_and_matter_BHNS_test_to_mpt(impt)
  end do
!
  call iteration_helmholtz_solver_binary_test(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!
  do impt = 1, nmpt
    call copy_from_mpatch_helmholtz_test(impt)
    call copy_metric_and_matter_BHNS_test_from_mpt(impt)
    if (outdata_type.eq.'3D') call IO_output_poisson_test_3D_mpt(impt)
  end do
!
!  call printout_NS_shape_seq(iseq)
!!  call printout_debug
!  call printout_NS_shape
!
END PROGRAM Main_helmholtz_binary_test
