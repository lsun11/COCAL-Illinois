!______________________________________________
!!include '../Include_file/include_modulefiles_peos.f90'
!!include '../Include_file/include_interface_modulefiles_peos.f90'
!!include '../Include_file/include_subroutines_peos.f90'
include '../Include_file/include_modulefiles_RNS_CF_peos.f90'
include '../Include_file/include_interface_modulefiles_RNS_CF_peos.f90'
include '../Include_file/include_subroutines_RNS_CF_peos.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_RNS_CF_peos
!
  use grid_parameter, only : indata_type, outdata_type, & 
  &                          iter_max, num_sol_seq, nrf_deform, &
  &                          sw_mass_iter, sw_art_deform, ntfxy, &
  &                          chrot, chgra, chope
  use def_matter, only : emd
  use def_matter_parameter, only : emdc
  use def_quantities, only : dhdr_x, dhdr_y, dhdr_z
  implicit none
  integer :: flag_restmass, count_adj
  integer :: iseq, iter_count, total_iteration, nmx
  character(40) :: char1, char2, char3, char4, char5
  character(LEN=300) :: cmd
!
  call coordinate_patch_kit_grav
  call read_parameter_drot
  call allocate_metric_and_matter
  call allocate_metric_on_SFC_CF
  call allocate_excurve_on_SFC_CF
  if (indata_type.eq.'1D') call IO_input_initial_1D
  if (indata_type.eq.'3D') call IO_input_initial_3D
  if (sw_mass_iter.eq.'y') then
    call search_emdmax_xaxis_grid(nmx)
    emdc = emd(nmx,ntfxy,0)
  end if
!
  if (sw_art_deform.eq.'y') call artificial_deformation
!
  call peos_initialize
!

  if (chrot=='p' .and. chgra=='e' .and. chope=='r') then
    iseq=1;   flag_restmass=999
    write(6,*)  "************************ ADDING A PERTURBATION ********************************"
    call input_perturbation_peos
    call iteration_peos(iter_count)
    call calc_physical_quantities_peos_v1
    call calc_quad_pole_peos
    write(6,*) "Basic quantities after perturbation is applied:"
    call printout_physq_console

    call printout_physq_peos_v1(iseq,flag_restmass)
    call printout_quad_pole(iseq)
    call printout_physq_plot(iseq,flag_restmass)
  else
    do iseq = 1, num_sol_seq
      total_iteration = 0
      flag_restmass   = 0 ! 0, 1 = iterate for a constant rest mass
!                         !    2 = rest mass converged
!                         !  999 = no rest mass iteration and converged
      count_adj       = 0 ! count number of adjustment of the rest mass
      do                  ! -- iteration for a constant rest mass
        call iteration_peos(iter_count)
        total_iteration = total_iteration + iter_count
        write(6,*)'-- Total # of iteration --', total_iteration
        call calc_physical_quantities_peos_v1
        call printout_physq_console
        if (total_iteration.ge.iter_max) exit
        if (sw_mass_iter.ne.'y') exit
        write(6,*) "Adjust mass", count_adj
        call adjust_rest_mass(flag_restmass,count_adj)
        if (flag_restmass.eq.2) exit
      end do
      write(6,*)'== Solution sequence # == ', iseq
      if (total_iteration.ge.iter_max) then
        write(6,*)' ** Solution did not converge **'
      else 
        if (sw_mass_iter.ne.'y') flag_restmass = 999
      end if
      call calc_physical_quantities_peos_v1
      call calc_quad_pole_peos
!
!      if ((dhdr_x/dhdr_z) < 0.0)  exit
!
!      call printout_physq_peos(iseq,flag_restmass)
      call printout_physq_peos_v1(iseq,flag_restmass)
      call printout_quad_pole(iseq)
      call printout_physq_plot(iseq,flag_restmass)

!     Print intermediate solution
!      if ( iseq<num_sol_seq .and. iseq>9 .and. (mod(iseq,5)==0) ) then
      if ( iseq>0 ) then
!      if ( iseq==24 ) then
        call open_directory(iseq,0)
        call IO_output_solution_3D
!        call system("rm -f rnsgrids_3D.las")
        write(char1, '(i5)') iseq;         char2 = adjustl(char1)
        write(char3, '(i5)') nrf_deform;   char4 = adjustl(char3)
        cmd = "echo 'iseq, nrf_deform : " // trim(char2) // "     " // trim(char4) // " ' > deformation.txt"
        call system(cmd)
        char5 = "iseq"//trim(char2)//trim(".dat")
        call printout_axis_CF(char5)
        call chdir('../')
      end if

      call printout_NS_shape_seq(iseq)
      call next_solution
    end do
  end if
  if (outdata_type.eq.'1D') call IO_output_solution_1D
  if (outdata_type.eq.'3D') call IO_output_solution_3D
!!  call printout_debug
  call printout_NS_shape
!
END PROGRAM Main_RNS_CF_peos
