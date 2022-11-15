!______________________________________________
include '../Include_file/include_modulefiles_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_modulefiles_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_subroutines_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_subroutines_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_irrot_BNS_CF_3mpt
!
  use grid_parameter_binary_excision
  use grid_points_asymptotic_patch
  use weight_midpoint_binary_excision
  use phys_constant
  use def_metric, only : psi
  use def_matter, only: vep, vepxf, vepyf, vepzf
  use def_metric_on_SFC_CF
  use def_metric_cartesian, only : psica
  use interface_modules_cartesian
  use interface_calc_gradvep
  use grid_points_binary_excision
  use trigonometry_grav_phi
  use def_binary_parameter, only : dis
  use def_matter_parameter, only : radi
  use def_matter_parameter_mpt
  implicit none
  integer :: impt, impt_ex
  character(30) :: char1
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
  character(len=2) :: cocp
!
! -- Read parameters
  call allocate_grid_parameter_mpt
  call allocate_grid_parameter_binary_excision_mpt
  call allocate_def_matter_parameter_mpt
  do impt = 1, nmpt
    call read_parameter_mpt(impt)
    call read_surf_parameter_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call read_parameter_binary_excision_mpt(impt)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
    call peos_initialize_mpt(impt)
    call copy_def_peos_parameter_to_mpt(impt)
  end do
!
! -- Allocate arrays
  call set_allocate_size_mpt
!
  call allocate_grid_points_binary_excision
  call allocate_grid_points_asymptotic_patch
  call allocate_metric_CF
  call allocate_metric_on_SFC_CF
  call allocate_matter
  call allocate_matter_velocity
  call allocate_grid_points_binary_excision_mpt
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_weight_midpoint_binary_excision
  call allocate_metric_CF_mpt
  call allocate_matter_mpt
  call allocate_matter_velocity_mpt
  call allocate_trig_grav_mphi
  call coordinate_patch_cartesian
  call allocate_metric_and_matter_cartesian
  call allocate_mpatch_all_BNS_CF
  call allocate_coordinate_patch_kit_grav_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_peos_parameter_from_mpt(impt)
    call coordinate_patch_kit_grav_grid_noGreen_mpt(3)  ! 3:r_surf is used
    call calc_parameter_binary_excision
    call calc_grid_points_binary_excision
    call calc_weight_midpoint_binary_excision
    call copy_to_mpatch_interpolation_utility(impt)
    call IO_input_converged_solution_3D_CF_irrot_NS_mpt(impt)
    call copy_def_metric_to_mpt(impt)
    call copy_def_matter_to_mpt(impt)
    call copy_grid_parameter_to_mpt(impt)
    call copy_def_peos_parameter_to_mpt(impt)
    call copy_def_binary_parameter_to_mpt(impt)
    call copy_def_matter_parameter_to_mpt(impt)
  end do
  if (nmpt.gt.2) then
    call copy_from_mpatch_interpolation_utility(nmpt)
    do impt = 1, nmpt
! -- coordinates of asymptotic patch in central patches
      call calc_grid_points_asymptotic_patch(impt,nmpt)
      call copy_grid_points_asymptotic_patch_to_mpt(impt)
    end do
  end if
!
! copy radi (i=5) of COCP1 to ARCP
  def_matter_param_real_(5,nmpt) = def_matter_param_real_(5,1)
!
  write(6,*) ""
  write(6,*) "================================================================================================"
  do impt = 1, nmpt
    call copy_from_mpatch_interpolation_utility(impt)
    call copy_def_metric_from_mpt(impt)
    call copy_def_matter_from_mpt(impt)
    if(impt<=2)  then
      cocp='ns'
    else
      cocp='bh'
    end if

    write(6,*) "Patch:", impt, ", writing xyz files..."
    call IO_output_plot_xyz_irrot_BNS_CF_mpt(impt)
    if (impt==1 .or. impt==2)  then
      call interpo_gr2fl_metric_CF
      call calc_gradvep(vep,vepxf,vepyf,vepzf)
      call reset_fluid_gradvep

!@      call interpolation_cartesian_irrot_BNS_mpt(impt)
!@      write(6,*) "Patch:", impt, ", writing contour files..."
!@      call IO_output_cartesian_contour_irrot_BNS_mpt(impt)
!@      write(6,*) "Patch:", impt, ", writing surface files..."
!@      call IO_output_surface_BNS_mpt(impt)
!@      write(6,*) "Patch:", impt, ", writing BNS shape files..."
!@      call printout_NS_shape_contour_mpt(impt)
!    
!cal      write(6,*) '-----------------------------------------------',dis
!cal      call calc_vector_x_matter(2)
!cal      call calc_vector_phi_matter(2)
!cal      call calc_vector_x_grav(2)
!cal      call calc_vector_phi_grav(2)
!cal      call excurve_CF(cocp)
!cal      call excurve_CF_gridpoint
!cal      call rotation_law
!cal      call calc_4velocity_ut_irrot_v1
!
!cal      call calc_rest_mass_peos_irrot      ! integral over star
!cal      call calc_mass_peos_irrot           ! integral over star only for Komar mass
!cal      call calc_proper_mass_peos_irrot    ! integral over star
!cal      call calc_ang_mom_peos_irrot        ! integral over star
!      call calc_ToverW_peos
!cal      call calc_radius_CF_rsurf
!      call calc_redblue_shift
!      call calc_quad_pole_peos
!cal      call calc_physq_center_peos_grid
!cal      call calc_physq_cgs_peos
!cal      call calc_enthalpy_xyzaxis
!cal      call calc_qua_loc_spin_grav
!cal      call calc_circ_line_peos_irrot(impt)
!cal      call calc_circ_surf_peos_irrot(impt)
!cal      call calc_soundspeed_peos
!cal      call IO_printout_grid_data_mpt(impt)
!cal      call printout_NS_shape_mpt(impt)
!cal      call copy_def_quantities_BNS_to_mpt(impt)
    end if
    if (impt==nmpt)  then
!cal      call calc_vector_x_grav(1)
!cal      call calc_vector_phi_grav(1)
!cal      call excurve_CF(cocp)
!cal      call excurve_CF_gridpoint_bhex   !in source_angmom_asympto
!
!cal      write(6,*) '-----------------------------------------------',dis
!cal      call calc_mass_asympto('ns')
!cal      call calc_angmom_asympto('ns')
!cal      call calc_admmom_asympto('ns')
!cal      call copy_def_quantities_BNS_to_mpt(impt)
    end if
  end do
!
!cal  char1 = 'bnsphyseq_plot.dat'
!cal  call printout_physq_BNS_all_mpt(char1)
!
END PROGRAM interpolation_contour_potential_irrot_BNS_CF_3mpt
