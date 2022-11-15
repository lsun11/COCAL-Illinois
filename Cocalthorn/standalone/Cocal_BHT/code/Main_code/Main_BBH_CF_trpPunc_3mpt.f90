!______________________________________________
include '../Include_file/include_modulefiles_BBH_CF_trpPunc_3mpt.f90'
include '../Include_file/include_interface_modulefiles_BBH_CF_trpPunc_3mpt.f90'
include '../Include_file/include_subroutines_BBH_CF_trpPunc_3mpt.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_BBH_CF_trpPunc_3mpt
!
  use phys_constant, only : nmpt
  use grid_parameter, only : indata_type, outdata_type, iter_max, &
  &                          num_sol_seq
  use def_bh_parameter, only : bh_soltype, mass_pBH
!  use grid_parameter_binary_excision, only : ex_nrg, ex_ndis
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_asymptotic_patch
  use grid_points_binary_in_asympto
  use weight_midpoint_binary_excision
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
!  use radial_green_fn_grav_bhex_di
  use radial_green_fn_grav_bhex_nd
  use radial_green_fn_grav_bhex_dh
  use radial_green_fn_grav_bhex_nh
  use radial_green_fn_grav_bhex_sd
  use interface_copy_to_hgfn_and_gfnsf
  implicit none
  integer :: impt, itg
  integer :: iseq, iter_count, total_iteration
!
! -- Read parameters
  call allocate_grid_parameter_mpt
  call allocate_grid_parameter_binary_excision_mpt
  call allocate_def_bh_parameter_mpt
  do impt = 1, nmpt
    call read_parameter_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call read_parameter_binary_excision_mpt(impt)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
    call read_parameter_pbh_mpt(impt)
    call copy_def_bh_parameter_to_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
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
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_nb
  call allocate_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
  call allocate_hgfn_bhex_sd
!  call allocate_hgfn_bhex_di
  call allocate_BBH_CF_AH
  call allocate_pBH_CF
!
  call allocate_mpatch_all_BBH_CF
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_grid_points_binary_in_asympto_mpt
  call allocate_BBH_CF_AH_mpt
  call allocate_pBH_CF_mpt
!
! -- Solution sequence
do iseq = 1, num_sol_seq
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_bh_parameter_from_mpt(impt)
    call coordinate_patch_kit_grav_mpt
    call calc_parameter_binary_excision
    if (impt.ne.nmpt) then 
      call calc_grid_points_binary_excision
      call calc_grid_points_binary_in_asympto(impt,nmpt)
      call copy_grid_points_binary_in_asympto_to_mpt(impt)
    end if
    call calc_weight_midpoint_binary_excision
!    call calc_weight_midpoint_binary_excision_hybrid
    call calc_vector_x_grav(2)
    call calc_vector_phi_grav(2)
    call calc_vector_bh(2)
!
    call copy_to_mpatch_all_BBH_CF(impt)
    call copy_def_metric_to_mpt(impt)
    call copy_def_metric_pBH_to_mpt(impt)
  end do
  call copy_from_mpatch_all_BBH_CF(nmpt)
! -- coordinates of asymptotic patch in central patches
  do impt = 1, 2
    call calc_grid_points_asymptotic_patch(impt,nmpt)
    call copy_grid_points_asymptotic_patch_to_mpt(impt)
  end do
if(iseq.eq.1) then
  do impt = 1, nmpt
    call copy_from_mpatch_all_BBH_CF(impt)
    if (indata_type.eq.'IN') call initial_metric_CF_pBH_mpt(impt)
    if (indata_type.eq.'3D') call IO_input_initial_3D_CF_BH_mpt(impt)
    call compute_alps2wmeN_mpt(impt)
    call copy_def_metric_to_mpt(impt)
    call copy_def_metric_pBH_to_mpt(impt)
    call initial_AHfinder
    call copy_def_horizon_to_mpt(impt)
  end do
end if
!
  call iteration_BBH_CF_trpPunc_3mpt(iter_count)
  write(6,*)'== Solution sequence # == ', iseq
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!
  call calc_physical_quantities_BBH_trpPunc_CF_mpt
  do impt = 1, nmpt
    call copy_from_mpatch_all_BBH_CF(impt)
    call printout_physq_BBH_trpPunc_mpt(iseq,impt)
  end do
!
  if (bh_soltype.eq.'SQ') call next_solution_BBH_mpt
  if (bh_soltype.ne.'SQ') exit
!
end do
!
  do impt = 1, nmpt
    call copy_from_mpatch_all_BBH_CF(impt)
    call copy_def_metric_from_mpt(impt)
    call copy_def_metric_pBH_from_mpt(impt)
    call copy_def_horizon_from_mpt(impt)
    if (outdata_type.eq.'3D') call IO_output_solution_3D_CF_BH_mpt(impt)
    call IO_output_plot_xyz_pBH_CF_mpt(impt)
    call IO_output_AHfinder_gnuplot_mpt(impt)
  end do
!
END PROGRAM Main_BBH_CF_trpPunc_3mpt
