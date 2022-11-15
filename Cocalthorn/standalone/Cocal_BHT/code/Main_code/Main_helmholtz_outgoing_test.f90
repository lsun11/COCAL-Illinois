!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Module/grid_parameter_binary_excision.f90'
include '../Module/grid_points_binary_excision.f90'
include '../Module/weight_midpoint_binary_excision.f90'
include '../Module/radial_green_fn_hrethadv.f90'
include '../Module/radial_green_fn_hret_mi_hadv.f90'
include '../Module/radial_green_fn_helmholtz.f90'
include '../Module/copy_array_3d.f90'
include '../Module/copy_array_4d.f90'
!
include '../Include_file/include_interface_modulefiles_peos.f90'
include '../Module_interface/interface_error_metric.f90'
include '../Module_interface/interface_grdr_gridpoint_type0.f90'
include '../Module_interface/interface_sourceterm_helmholtz_solver_test.f90'
include '../Module_interface/interface_helmholtz_solver_outgoing.f90'
include '../Module_interface/interface_helmholtz_solver_vol_int.f90'
include '../Module_interface/interface_helmholtz_solver_vol_int_hret_mi_hadv.f90'
include '../Module_interface/interface_copy_to_bsjy_and_sbsjy.f90'
!
include '../Include_file/include_subroutines_peos.f90'
include '../Subroutine/IO_output_poisson_test_3D.f90'
include '../Subroutine/iteration_helmholtz_solver_outgoing_test.f90'
include '../Subroutine/error_metric.f90'
include '../Subroutine/printout_error_metric.f90'
include '../Subroutine/sourceterm_helmholtz_solver_test.f90'
include '../Subroutine/helmholtz_solver_outgoing.f90'
include '../Subroutine/helmholtz_solver_vol_int.f90'
include '../Subroutine/helmholtz_solver_vol_int_hret_mi_hadv.f90'
include '../Subroutine/grdr_gridpoint_type0.f90'
include '../Subroutine/test_source_helical.f90'
include '../Subroutine/allocate_poisson_solver_test.f90'
include '../Subroutine/test_analytic_solution.f90'
include '../Subroutine/calc_radial_green_fn_hrethadv.f90'
include '../Subroutine/calc_radial_green_fn_hret_mi_hadv.f90'
include '../Subroutine/copy_to_bsjy_and_sbsjy.f90'
include '../Subroutine/sphbess.f90'
include '../Subroutine/sphbess_and_dx.f90'
include '../Subroutine/bessjy.f90'
include '../Subroutine/beschb.f90'
include '../Function/chebev.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_helmholtz_outgoing_test
!
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max
  use radial_green_fn_helmholtz
  use radial_green_fn_hrethadv
  use radial_green_fn_hret_mi_hadv
  implicit none
  integer :: iseq, iter_count, total_iteration
!
  call coordinate_patch_kit_grav
  call allocate_poisson_solver_test
  call allocate_radial_green_fn_helmholtz
  call allocate_radial_green_fn_hrethadv
  call allocate_radial_green_fn_hret_mi_hadv
  call test_source_helical
! 
  call iteration_helmholtz_solver_outgoing_test(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
  if (outdata_type.eq.'3D') call IO_output_poisson_test_3D
!
END PROGRAM Main_helmholtz_outgoing_test
