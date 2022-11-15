!______________________________________________
!include '../Include_file/include_modulefiles_MRNS.f90'
!include '../Include_file/include_interface_modulefiles_MRNS.f90'
!include '../Include_file/include_PEOS_modulefile.f90'
!include '../Include_file/include_PEOS_subroutines.f90'
!include '../Include_file/include_subroutines_MRNS.f90'
!include '../Include_file/include_functions.f90'
!
include '../Include_file/include_modulefiles_plot_MRNS.f90'
include '../Include_file/include_modulefiles_analysis_plot_MRNS.f90'
include '../Include_file/include_interface_modulefiles_plot_MRNS.f90'
include '../Include_file/include_interface_modulefiles_analysis_plot_MRNS.f90'
include '../Include_file/include_subroutines_plot_MRNS.f90'
include '../Include_file/include_subroutines_analysis_plot_MRNS.f90'
include '../Include_file/include_functions.f90'
!
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM interpolation_contour_potential_MRNS
!
  use weight_midpoint_grav
  implicit none
!
  call coordinate_patch_kit_grav_noGreen
  call coordinate_patch_cartesian
  call allocate_metric_and_matter_WL_MHD
  call allocate_metric_and_matter_WL_cartesian
  call allocate_matter_velocity
  call allocate_emfield
  call allocate_emfield_cartesian
  call allocate_weight_midpoint_grav
  call weight_calc_midpoint_grav
!
  call IO_input_initial_3D
  call IO_input_initial_3D_WL
  call IO_input_initial_3D_4ve
  call IO_input_initial_3D_EMF
!
  call interpolation_cartesian_MRNS
  call IO_output_cartesian_contour_CF
  call IO_output_cartesian_contour_WL
  call IO_output_cartesian_contour_FLU
  call IO_output_cartesian_contour_EMF
  call IO_output_surface
  call printout_NS_shape_seq(1)
!
END PROGRAM interpolation_contour_potential_MRNS
