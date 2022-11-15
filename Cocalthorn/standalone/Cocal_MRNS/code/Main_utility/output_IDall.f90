!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Include_file/include_interface_utilities.f90'
include '../Include_file/include_subroutines_utilities.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!
include '../Module/def_matter_velocity.f90'
!
include '../Subroutine/allocate_metric_CF.f90'
include '../Subroutine/allocate_matter_emdrsrho.f90'
include '../Subroutine/allocate_matter_3velocity.f90'
include '../Subroutine/allocate_matter_4velocity.f90'
include '../Subroutine/calc_3velocity_corot.f90'
include '../Subroutine/calc_4velocity_corot.f90'
include '../Subroutine/calc_matter_rhof.f90'
!
include '../Subroutine/interpolation_fl2gr_IDall.f90'
include '../Subroutine/IO_input_converged_solution_3D.f90'
include '../Subroutine/IO_input_converged_solution_fluid_3D.f90'
include '../Subroutine/IO_output_solution_IDall_3D.f90'
!
include '../Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM output_ID_all
!
  use def_vector_phi
  implicit none
!
  call coordinate_patch_kit_grav_noGreen
  call peos_initialize
  call allocate_metric_CF
  call allocate_matter_emdrsrho
  call allocate_matter_3velocity
  call allocate_matter_4velocity
  call allocate_vector_phi
  call IO_input_converged_solution_3D
  call IO_input_converged_solution_fluid_3D
!
  call calc_vector_phi_matter(1)
  call calc_3velocity_corot
  call calc_4velocity_corot
  call calc_matter_rhof
  call interpolation_fl2gr_IDall
!
  call IO_output_solution_IDall_3D
!
END PROGRAM output_ID_all
