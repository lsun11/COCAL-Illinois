!______________________________________________
include '../Include_file/include_modulefiles_RNS_CF_peos_plot.f90'
include '../Include_file/include_modulefiles_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_RNS_CF_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_subroutines_RNS_CF_peos_plot.f90'
include '../Include_file/include_subroutines_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_main
!
  implicit none
!
  call coordinate_patch_kit_grav_noGreen
  call coordinate_patch_cartesian
  call allocate_metric_and_matter
  call allocate_metric_and_matter_cartesian
  call allocate_matter_velocity
  call IO_input_initial_3D
!
  call interpolation_cartesian
  call IO_output_cartesian_contour
  call IO_output_surface
  call IO_output_plot_xyz
  call printout_NS_shape_seq(1)
!
END PROGRAM interpolation_main
