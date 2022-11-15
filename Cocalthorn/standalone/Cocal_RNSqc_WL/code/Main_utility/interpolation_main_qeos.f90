!______________________________________________
include '../Include_file/include_modulefiles_RQS_CF_qeos_plot.f90'
include '../Include_file/include_modulefiles_analysis_RQS_CF_qeos_plot.f90'
include '../Include_file/include_interface_modulefiles_RQS_CF_qeos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_RQS_CF_qeos_plot.f90'
include '../Include_file/include_subroutines_RQS_CF_qeos_plot.f90'
include '../Include_file/include_subroutines_analysis_RQS_CF_qeos_plot.f90'
include '../Include_file/include_QEOS_modulefile.f90'
include '../Include_file/include_QEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_main_qeos
! 
  use def_matter_parameter, only : rhos_qs
  use def_qeos_parameter 
  implicit none
!
  call qeos_initialize
  rhos_qs = rhosurf_gcm1
  call coordinate_patch_kit_grav_noGreen
  call coordinate_patch_cartesian
  call allocate_metric_and_matter_qeos
  call allocate_metric_and_matter_cartesian_qeos
  call allocate_matter_velocity
  call IO_input_initial_3D_qeos
!
  call interpolation_cartesian_qeos
  call IO_output_cartesian_contour_qeos
  call IO_output_surface
  call IO_output_plot_xyz_qeos
  call printout_NS_shape_seq(1)
!
END PROGRAM interpolation_main_qeos
