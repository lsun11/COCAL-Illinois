!______________________________________________
include '../Include_file/include_modulefiles_peos.f90'
include '../Include_file/include_interface_utilities.f90'
include '../Include_file/include_subroutines_utilities.f90'
!
include '../Module/def_matter_velocity.f90'
!
include '../Subroutine/allocate_matter_emdrsrhof.f90'
include '../Subroutine/allocate_matter_3velocity.f90'
include '../Subroutine/allocate_matter_4velocity.f90'
include '../Subroutine/calc_3velocity_corot.f90'
include '../Subroutine/calc_4velocity_corot.f90'
include '../Subroutine/calc_matter_rhof.f90'
!
include '../Subroutine/IO_input_converged_solution_fluid_3D.f90'
include '../Subroutine/IO_output_solution_fluid_3D.f90'
!
include '../Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90'
!______________________________________________
!
!              Interpolation Program
!______________________________________________
PROGRAM output_fluid_velocity
!
  use def_vector_phi
  implicit none
!
  call coordinate_patch_kit_grav_noGreen
  call peos_initialize
  call allocate_matter_emdrsrhof
  call allocate_matter_3velocity
  call allocate_matter_4velocity
  call allocate_vector_phi
  call IO_input_converged_solution_fluid_3D
!
  call calc_vector_phi_matter(1)
  call calc_3velocity_corot
  call calc_4velocity_corot
  call calc_matter_rhof
!
  call IO_output_solution_fluid_3D
!
END PROGRAM output_fluid_velocity
