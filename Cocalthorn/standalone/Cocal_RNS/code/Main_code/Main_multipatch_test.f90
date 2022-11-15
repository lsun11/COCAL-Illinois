!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_grav_bhex_di.f90'
include '../Module/radial_green_fn_grav_bhex_nb.f90'
include '../Module/copy_array_3d.f90'
include '../Module/make_char2_array_2d.f90'
include '../Module/make_char1_array_2d.f90'
include '../Module/make_int_array_3d.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_grdr_gridpoint_type0_nosym.f90'
include '../Module_interface/interface_sourceterm_poisson_solver_test.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_surface_int.f90'
include '../Module_interface/interface_poisson_solver_binary_bhex.f90'
include '../Module_interface/interface_poisson_solver_binary_surf_int.f90'
include '../Module_interface/interface_poisson_solver_binary_vol_int.f90'
include '../Module_interface/interface_poisson_solver_bhex_surf_int.f90'
include '../Module_interface/interface_copy_to_hgfn_and_gfnsf.f90'
include '../Include_file/include_modulefiles_mpatch.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
include '../Subroutine/iteration_poisson_bbh_test.f90'
include '../Subroutine/error_metric.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/coordinate_patch_kit_bhex.f90'
include '../Subroutine/copy_hgfn_di_to_hgfn.f90'
include '../Subroutine/copy_to_hgfn_and_gfnsf.f90'
include '../Subroutine/sourceterm_poisson_solver_test.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_surface_int.f90'
include '../Subroutine/poisson_solver_binary_bhex.f90'
include '../Subroutine/poisson_solver_binary_vol_int.f90'
include '../Subroutine/poisson_solver_binary_surf_int.f90'
include '../Subroutine/poisson_solver_bhex_surf_int.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/grdr_gridpoint_type0_nosym.f90'
include '../Subroutine/reset_bh_boundary.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/test_analytic_solution.f90'
!include '../Module_mpatch/grid_parameter_mpt.f90'
include '../Include_file/include_subroutines_mpatch.f90'
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
PROGRAM Main_multipatch_test
!
  use phys_constant, only : nmpt
  use grid_parameter, only :  nrg, mass_eps, sw_art_deform
!  use grid_parameter_binary_excision, only : ex_nrg, ex_ndis
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use interface_copy_to_hgfn_and_gfnsf
!-- 
!
  implicit none
  integer :: impt
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
  call allocate_weight_midpoint_binary_excision
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_nb
  call allocate_poisson_solver_test
!
  call allocate_mpatch_all_test
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call coordinate_patch_kit_grav_mpt
    call calc_grid_points_binary_excision
    call calc_weight_midpoint_binary_excision
    call calc_hgfn_bhex_nb
    call copy_to_hgfn_and_gfnsf(hgfn_nb,gfnsf_nb)
    call calc_vector_x_grav(1)
    call calc_vector_x_matter(1)
    call calc_vector_phi_grav(1)
    call calc_vector_phi_matter(1)
    call copy_to_mpatch_all_test(impt)
  end do
!
  call copy_from_mpatch_all_test(1)
  call copy_from_mpatch_all_test(2)
!  call copy_to_mpatch_all_test
!  call copy_from_mpatch_all_test
!
!  call copy_grid_parameter_from_mpt(1)
!  write(6,*) nrg, mass_eps, sw_art_deform
!  call copy_grid_parameter_from_mpt(2)
!  write(6,*) nrg, mass_eps, sw_art_deform 
  !call coordinate_patch_kit_grav
  !call read_parameter_binary_excision
  !call calc_parameter_binary_excision
  !call allocate_grid_points_binary_excision
  !call calc_grid_points_binary_excision
!
END PROGRAM Main_multipatch_test
