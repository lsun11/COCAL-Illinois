!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_helmholtz.f90'
include '../Module/radial_green_fn_hrethadv.f90'
include '../Module/copy_array_3d.f90'
include '../Module/copy_array_4d.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_sourceterm_helmholtz_solver_test.f90'
include '../Module_interface/interface_sourceterm_exsurf_eqm_binary.f90'
include '../Module_interface/interface_sourceterm_outsurf_eqm_binary.f90'
include '../Module_interface/interface_helmholtz_solver_binary.f90'
include '../Module_interface/interface_helmholtz_solver_binary_vol_int.f90'
include '../Module_interface/interface_helmholtz_solver_binary_surf_int.f90'
include '../Module_interface/interface_helmholtz_solver_outer_surf_int.f90'
include '../Module_interface/interface_copy_to_bsjy_and_sbsjy.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
include '../Subroutine/iteration_helmholtz_solver_binary_test.f90'
include '../Subroutine/error_metric.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/sourceterm_helmholtz_solver_test.f90'
include '../Subroutine/sourceterm_exsurf_eqm_binary.f90'
include '../Subroutine/sourceterm_outsurf_eqm_binary.f90'
include '../Subroutine/helmholtz_solver_binary.f90'
include '../Subroutine/helmholtz_solver_binary_vol_int.f90'
include '../Subroutine/helmholtz_solver_binary_surf_int.f90'
include '../Subroutine/helmholtz_solver_outer_surf_int.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/test_source_helical_binary.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/test_analytic_solution.f90'
include '../Subroutine/calc_radial_green_fn_hrethadv.f90'
include '../Subroutine/copy_to_bsjy_and_sbsjy.f90'
include '../Subroutine/sphbess.f90'
include '../Subroutine/sphbess_and_dx.f90'
include '../Subroutine/sphbess_dx.f90'
include '../Subroutine/bessjy.f90'
include '../Subroutine/beschb.f90'
include '../Function/chebev.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_helmholtz_binary_test
!
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  use radial_green_fn_helmholtz
  use radial_green_fn_hrethadv
  implicit none
  integer :: iseq, iter_count, total_iteration
!
  call coordinate_patch_kit_grav
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call allocate_weight_midpoint_binary_excision
  call calc_weight_midpoint_binary_excision
  call allocate_poisson_solver_test
  call allocate_radial_green_fn_helmholtz
  call allocate_radial_green_fn_hrethadv
  call test_source_helical_binary
! call test_analytic_solution
! 
  call iteration_helmholtz_solver_binary_test(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!  call printout_NS_shape_seq(iseq)
  if (outdata_type.eq.'3D') call IO_output_poisson_test_3D
!!  call printout_debug
!  call printout_NS_shape
!
END PROGRAM Main_helmholtz_binary_test
