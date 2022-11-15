!______________________________________________
include '../Include_file/include_modulefiles_BBH_CF_circ.f90'
include '../Include_file/include_interface_modulefiles_BBH_CF_circ.f90'
include '../Include_file/include_subroutines_BBH_CF_circ.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_BBH_CF_circ
!
  use grid_parameter, only : indata_type, outdata_type, iter_max
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
  use radial_green_fn_grav_bhex_nd
  use radial_green_fn_grav_bhex_rd
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use weight_midpoint_binary_excision
  use interface_copy_to_hgfn_and_gfnsf
  implicit none
  integer :: iseq, iter_count, total_iteration, flag_omega, count_adj
!
  call coordinate_patch_kit_bhex
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_dd
  call calc_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
  call calc_hgfn_bhex_nd
  call allocate_hgfn_bhex_rd
  call calc_hgfn_bhex_rd
! -- No boundary Green's fn
  call allocate_hgfn_bhex_nb
  call calc_hgfn_bhex_nb
!  call copy_to_hgfn_and_gfnsf(hgfn_nb,gfnsf_nb)
! --
  call read_parameter_bh
  call read_parameter_binary_excision
  call calc_parameter_binary_excision
  call IO_printout_grid_data
  call allocate_grid_points_binary_excision
  call calc_grid_points_binary_excision
  call allocate_weight_midpoint_binary_excision
  call calc_weight_midpoint_binary_excision
!  call calc_weight_midpoint_binary_excision_hybrid
  call allocate_BBH_CF_AH
  if (indata_type.eq.'IN') call initial_metric_CF
  if (indata_type.eq.'3D') call IO_input_initial_3D
!  call test_analytic_solution_bhex_psialph
!
  total_iteration = 0

  call calc_circular_orbit(iter_count)
!
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!
END PROGRAM Main_BBH_CF_circ
