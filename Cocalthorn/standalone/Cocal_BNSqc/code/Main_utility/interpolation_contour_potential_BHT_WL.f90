!
include '../Include_file/include_modulefiles_BHT_WL_peos_plot.f90'
include '../Include_file/include_modulefiles_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_BHT_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_subroutines_BHT_WL_peos_plot.f90'
include '../Include_file/include_subroutines_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_BHT_WL
!
  use grid_parameter
  implicit none
!
  call read_parameter
  call read_bht_parameter
  call calc_bht_excision_radius

  call coordinate_patch_kit_grav_grid_noGreen(5)   ! 5 for BHT
!  call IO_printout_grid_data_noex
  call read_parameter_drot
  call allocate_metric_3p1_WL_CTT
  call allocate_BHT_matter_WL

  call peos_initialize

  call coordinate_patch_cartesian
  call allocate_metric_and_matter_WL_cartesian
!
  write(6,*) "Reading 3D initial data..."
  call IO_input_initial_3D_BHT
  call IO_input_initial_3D_WL_BHT
!  call IO_input_initial_Kij_3D     needs work
!
  write(6,*) "Interpolating on the cartesian grid..."
  call interpolation_cartesian_BHT_WL
  call IO_output_cartesian_contour_CF
  call IO_output_cartesian_contour_WL
  call IO_output_cartesian_contour_BHT
!
END PROGRAM interpolation_contour_potential_BHT_WL
