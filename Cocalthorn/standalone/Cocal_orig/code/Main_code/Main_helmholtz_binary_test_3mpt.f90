!______________________________________________
include '../Include_file/include_modulefiles_helmholtz_binary_test_3mpt.f90'
include '../Include_file/include_interface_modulefiles_helmholtz_binary_test_3mpt.f90'
include '../Include_file/include_subroutines_helmholtz_binary_test_3mpt.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_helmholtz_binary_test_3mpt
!
  use phys_constant, only : nmpt
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_binary_excision_mpt
  use grid_points_asymptotic_patch
  use grid_points_binary_in_asympto
  use weight_midpoint_grav
  use weight_midpoint_fluid
  use weight_midpoint_binary_excision
  use def_vector_x
  use def_vector_phi
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
  use radial_green_fn_grav_bhex_sd
  use radial_green_fn_helmholtz
  use radial_green_fn_hrethadv
  use radial_green_fn_hrethadv_homosol
  use radial_green_fn_hret_mi_hadv
  use radial_green_fn_hret_mi_hadv_homosol
  implicit none
  integer :: iseq, iter_count, total_iteration, impt, impt_bin
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
    call read_parameter_bh_mpt(impt)
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
  call allocate_hgfn_bhex_sd
  call allocate_radial_green_fn_helmholtz
  call allocate_radial_green_fn_hrethadv
  call allocate_radial_green_fn_hrethadv_homosol
  call allocate_radial_green_fn_hret_mi_hadv
  call allocate_radial_green_fn_hret_mi_hadv_homosol
  call allocate_metric_and_matter_BHNS_test
!!
  call allocate_mpatch_all_BBH_CF
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_grid_points_binary_in_asympto_mpt
  call allocate_metric_and_matter_BHNS_test_mpt
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
    call test_source_helical_binary_mpt(impt)
    call calc_weight_midpoint_binary_excision
    call calc_vector_x_grav(0)
    call calc_vector_x_matter(0)
    call calc_vector_phi_grav(0)
    call calc_vector_phi_matter(0)
    call calc_vector_bh(0)
!
    call copy_to_mpatch_all_BBH_CF(impt)
    call copy_metric_and_matter_BHNS_test_to_mpt(impt)
  end do
!
  call copy_from_mpatch_all_BBH_CF(nmpt)
! -- coordinates of asymptotic patch in central patches
  do impt = 1, 2
    call calc_grid_points_asymptotic_patch(impt,nmpt)
    call copy_grid_points_asymptotic_patch_to_mpt(impt)
  end do
!
  call iteration_helmholtz_solver_binary_test_mpt(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!
  call printout_physq_BBH_mpt(1)
  do impt = 1, nmpt
    call copy_from_mpatch_all_BBH_CF(impt)
    call copy_metric_and_matter_BHNS_test_from_mpt(impt)
    if (outdata_type.eq.'3D') call IO_output_poisson_test_3D_mpt(impt)
  end do
!
!  call printout_NS_shape_seq(iseq)
!!  call printout_debug
!  call printout_NS_shape
!
END PROGRAM Main_helmholtz_binary_test_3mpt
