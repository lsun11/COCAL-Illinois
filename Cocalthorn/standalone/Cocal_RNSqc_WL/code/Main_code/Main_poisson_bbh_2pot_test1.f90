!______________________________________________
include '../Include_file/include_modulefiles_poisson_bbh_2pot_test1.f90'
include '../Include_file/include_interface_modulefiles_poisson_bbh_2pot_test1.f90'
include '../Include_file/include_subroutines_poisson_bbh_2pot_test1.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_poisson_bbh_2pot_test1
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
  use def_bh_parameter, only: bh_sptype
  implicit none
  integer :: iseq, iter_count, total_iteration
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
  call allocate_poisson_bbh_test
  if (indata_type.eq.'IN') call initial_metric_CF
  if (indata_type.eq.'3D') call IO_input_initial_3D
!  call test_analytic_solution_bhex_psialph
!
  total_iteration = 0

  call iteration_poisson_bbh_2pot_test1(total_iteration)

  if (outdata_type.eq.'3D') call IO_output_solution_3D
!
END PROGRAM Main_poisson_bbh_2pot_test1
