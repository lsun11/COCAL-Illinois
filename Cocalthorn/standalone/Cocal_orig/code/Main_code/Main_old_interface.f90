!______________________________________________
include '../Include_file/include_modulefiles.f90'
include '../Include_file/include_subroutines.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM main
!
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max, num_sol_seq,  &
  &                          sw_mass_iter, sw_art_deform
  use def_matter, only : emd
  use def_matter_parameter, only : emdc
  implicit none
  integer :: flag_restmass, count_adj
  integer :: iseq, iter_count, total_iteration
!
  call coordinate_patch_kit_grav
  call allocate_metric_and_matter
  if (indata_type.eq.'1D') call IO_input_initial_1D
  if (indata_type.eq.'3D') call IO_input_initial_3D
  if (sw_mass_iter.eq.'y') emdc = emd(0,0,0)
!
  if (sw_art_deform.eq.'y') call artificial_deformation
!
  do iseq = 1, num_sol_seq
    total_iteration = 0
    flag_restmass   = 0 ! 0, 1 = iterate for a constant rest mass
!                       !    2 = rest mass converged
!                       !  999 = no rest mass iteration and converged
    count_adj       = 0 ! count number of adjustment of the rest mass
    do                  ! -- iteration for a constant rest mass
      call iteration(iter_count)
      total_iteration = total_iteration + iter_count
      write(6,*)'-- Total # of iteration --', total_iteration
      call calc_physical_quantities
      call printout_physq_console
      if (total_iteration.ge.iter_max) exit
      if (sw_mass_iter.ne.'y') exit
      call adjust_rest_mass(flag_restmass,count_adj)
      if (flag_restmass.eq.2) exit
    end do
    write(6,*)'== Solution sequence # == ', iseq
    if (total_iteration.ge.iter_max) then
      write(6,*)' ** Solution did not converge **'
    else 
      if (sw_mass_iter.ne.'y') flag_restmass = 999
    end if
    call calc_physical_quantities
    call printout_physq(iseq,flag_restmass)
    call printout_physq_plot(iseq,flag_restmass)
!testtest    if (sw_mass_iter.ne.'y') exit
    call printout_NS_shape_seq(iseq)
    call next_solution
  end do
  if (outdata_type.eq.'1D') call IO_output_solution_1D
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!!  call printout_debug
  call printout_NS_shape
!
END PROGRAM main
