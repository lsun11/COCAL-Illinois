!______________________________________________
include '../Include_file/include_modulefiles_BBH_CF.f90'
include '../Include_file/include_interface_modulefiles_BBH_CF.f90'
include '../Include_file/include_subroutines_BBH_CF.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_BBH_CF
!
  use grid_parameter, only : indata_type, outdata_type, iter_max
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
  call read_parameter_bh
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call allocate_weight_midpoint_binary_excision
  call calc_weight_midpoint_binary_excision
  call allocate_BBH_CF
  if (indata_type.eq.'IN') call initial_metric_CF
  if (indata_type.eq.'3D') call IO_input_initial_3D
!  call test_analytic_solution_bhex_psialph
! 
  call iteration_BBH_CF(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!  call printout_NS_shape_seq(iseq)
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!!  call printout_debug
!  call printout_NS_shape
!
END PROGRAM Main_BBH_CF
