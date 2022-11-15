!______________________________________________
include '../Include_file/include_modulefiles_pBH.f90'
include '../Include_file/include_interface_modulefiles_pBH.f90'
include '../Include_file/include_subroutines_pBH.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_pBH_CF
!
  use grid_parameter, only : indata_type, outdata_type, iter_max
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
  use radial_green_fn_grav_bhex_nd
  use radial_green_fn_grav_bhex_sd
  use interface_copy_to_hgfn_and_gfnsf
  implicit none
  integer :: iseq, iter_count, total_iteration
!
  call read_parameter_pbh
  call coordinate_patch_kit_bhex
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_dd
  call calc_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
  call calc_hgfn_bhex_nd
  call allocate_hgfn_bhex_sd
  call calc_hgfn_bhex_sd
! -- No boundary Green's fn
  call allocate_hgfn_bhex_nb
  call calc_hgfn_bhex_nb
  call copy_to_hgfn_and_gfnsf(hgfn_nb,gfnsf_nb)
! --
  call allocate_BBH_CF
  call allocate_pBH_CF
  if (indata_type.eq.'IN') call initial_metric_CF_pBH
  if (indata_type.eq.'3D') call IO_input_initial_3D
! 
  call iteration_pBH_CF(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!
  call allocate_horizon
  call initial_AHfinder
  call excurve_TrpBH
  call excurve_TrpBH_gridpoint
  call copy_Aij_pBH_to_tfkij
  call iteration_AHfinder(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge AH finder **'
  end if
  call copy_Aij_pBH_to_tfkij
  call calc_AHarea_AHfinder
  call IO_output_AHfinder
  call IO_output_AHfinder_gnuplot
!
  iseq = 1
  call calc_physical_quantities_BH
  call printout_physq_BH(iseq)
  call IO_output_plot_xyz_pBH_CF
!
END PROGRAM Main_pBH_CF
