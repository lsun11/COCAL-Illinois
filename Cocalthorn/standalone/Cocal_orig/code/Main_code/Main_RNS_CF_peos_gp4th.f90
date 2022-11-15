!______________________________________________
!!include '../Include_file/include_modulefiles_peos.f90'
!!include '../Include_file/include_interface_modulefiles_peos.f90'
!!include '../Include_file/include_subroutines_peos.f90'
include '../Include_file/include_modulefiles_RNS_CF_peos.f90'
include '../Include_file/include_interface_modulefiles_RNS_CF_peos.f90'
include '../Include_file/include_subroutines_RNS_CF_peos_gp4th.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_RNS_CF_peos_gp4th
!
  use grid_parameter 
  use def_matter, only : emd
  use def_matter_parameter, only : emdc
  use def_quantities, only : dhdr_x, dhdr_y, dhdr_z
  use def_peos_parameter, only : emdini_gcm1
  use interface_violation_midpoint_MoC_CF_peos
  use interface_violation_gridpoint_MoC_CF_peos
  use interface_violation_midpoint_HaC_CF_peos
  use interface_violation_gridpoint_HaC_CF_peos
  use interface_IO_output_3D_general
  use interface_IO_output_2D_general
  use interface_IO_output_1D_general
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
  integer :: flag_restmass, count_adj
  integer :: iseq, iter_count, total_iteration
  real(long), pointer :: pot(:,:,:), HaC_vio(:,:,:), MoC_vio(:,:,:,:)
  character(30) :: char1, char2, char3, char4, char5
!
  call coordinate_patch_kit_grav
  call read_parameter_drot
  call allocate_metric_and_matter
  call allocate_metric_on_SFC_CF
  if (indata_type.eq.'1D') call IO_input_initial_1D
  if (indata_type.eq.'3D') call IO_input_initial_3D
  if (sw_mass_iter.eq.'y') emdc = emd(0,0,0)
!
  if (sw_art_deform.eq.'y') call artificial_deformation
!
  call peos_initialize

!emdc = emdini_gcm1
!
  do iseq = 1, num_sol_seq
    total_iteration = 0
    flag_restmass   = 0 ! 0, 1 = iterate for a constant rest mass
!                       !    2 = rest mass converged
!                       !  999 = no rest mass iteration and converged
    count_adj       = 0 ! count number of adjustment of the rest mass
    do                  ! -- iteration for a constant rest mass
      call iteration_peos_gp4th(iter_count)
      total_iteration = total_iteration + iter_count
      write(6,*)'-- Total # of iteration --', total_iteration
      call calc_physical_quantities_peos
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
    call calc_physical_quantities_peos
    call calc_quad_pole_peos(iseq)
!
!    if ((dhdr_x/dhdr_z) < 0.0)  exit
!
    call printout_physq_peos(iseq,flag_restmass)
    call printout_physq_peos_v1(iseq,flag_restmass)
    call printout_quad_pole(iseq)
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
  call alloc_array4d(MoC_vio,0,nrg,0,ntg,0,npg,1,3)
  call alloc_array3d(HaC_vio,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call excurve_CF('ns')             !   3rd order midpoint from ir0=-2,...
  call excurve_CF_gridpoint         !   4th order from ir0=-2,...

  MoC_vio(0:nrg,0:ntg,0:npg,1:3) = 0.0d0
  call violation_gridpoint_MoC_CF_peos(MoC_vio)
  pot = 0.0d0
  pot(0:nrg,0:ntg,0:npg) = MoC_vio(0:nrg,0:ntg,0:npg,2)
  char1 = 'MoC_by_3D.txt'
  call IO_output_3D_general(char1,'g','g',pot)
  char1 = 'MoC_by_xy.txt'
  call IO_output_2D_general(char1,'g','g',pot,'xy')
  char1 = 'MoC_by_phi000.txt'
  call IO_output_1D_general(char1,'g','g',pot,-1,ntg/2,0)
  char1 = 'MoC_by_phi180.txt'
  call IO_output_1D_general(char1,'g','g',pot,-1,ntg/2,npg/2)

  HaC_vio = 0.0d0
  call violation_gridpoint_HaC_CF_peos(HaC_vio)
  char1 = 'HaC_3D.txt'
  call IO_output_3D_general(char1,'g','g',HaC_vio)
  char1 = 'HaC_xy.txt'
  call IO_output_2D_general(char1,'g','g',HaC_vio,'xy')
  char1 = 'HaC_phi000.txt'
  call IO_output_1D_general(char1,'g','g',HaC_vio,-1,ntg/2,0)
  char1 = 'HaC_phi180.txt'
  call IO_output_1D_general(char1,'g','g',HaC_vio,-1,ntg/2,npg/2)

  HaC_vio = 0.0d0
  call violation_midpoint_HaC_CF_peos(HaC_vio)
  char1 = 'HaC_3D_mid.txt'
  call IO_output_3D_general(char1,'g','m',HaC_vio)
  char1 = 'HaC_xy_mid.txt'
  call IO_output_2D_general(char1,'g','m',HaC_vio,'xy')
  char1 = 'HaC_phi000_mid.txt'
  call IO_output_1D_general(char1,'g','m',HaC_vio,-1,ntg/2,1)
  char1 = 'HaC_phi180_mid.txt'
  call IO_output_1D_general(char1,'g','m',HaC_vio,-1,ntg/2,npg/2)

  deallocate(MoC_vio)
  deallocate(HaC_vio)
  deallocate(pot)
!
END PROGRAM Main_RNS_CF_peos_gp4th
