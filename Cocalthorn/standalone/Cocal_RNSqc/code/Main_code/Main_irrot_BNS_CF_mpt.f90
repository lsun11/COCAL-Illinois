!______________________________________________
include '../Include_file/include_modulefiles_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_interface_modulefiles_BNS_CF_3mpt.f90'
include '../Include_file/include_subroutines_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM Main_irrot_BNS_CF_mpt
!
  use phys_constant, only : long, nmpt
  use grid_parameter
  use def_matter
  use def_matter_parameter, only : emdc, radi, ome
  use def_matter_parameter_mpt
  use def_quantities_mpt
  use def_peos_parameter, only : rhoini_cgs, emdini_gcm1
!  use def_bh_parameter, only : bh_soltype
  use grid_parameter_binary_excision
  use grid_points_binary_excision
  use grid_points_asymptotic_patch
  use grid_points_binary_in_asympto
  use weight_midpoint_binary_excision
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_dd
!  use radial_green_fn_grav_bhex_di
  use radial_green_fn_grav_bhex_nd
  use radial_green_fn_grav_bhex_dh
  use radial_green_fn_grav_bhex_nh
  use radial_green_fn_grav_bhex_sd
  use interface_copy_to_hgfn_and_gfnsf
  use interface_calc_gradvep
  use interface_calc_gradvep_from_corot_id
  implicit none

!------------------------------------------------------------------------------------------------------------!
!                                                                                                            !
! numseq=1 for a single solution with given TOV initial data or given 3D initial data                        !
! numseq>1 produces a sequence of solutions with increasing emdc using routine next_solution_BNS_mpt.f90     !
! In this case the iteration starts from TOV initial data (or 3D) computes the first irrotating solution     !
! and then increases emdc to move to the next irrotating solution.                                           !
! When computing parallel irrotating solutions starting from different TOV use numseq=1 and the script       !
! debug_BNS_CF_peos_seq_from_TOV.sh                                                                          !
!                                                                                                            !
! sw_mass_iter = 'n' implements the above scenarios                                                          !
! sw_mass_iter = 'y' finds a solution for a given rest mass.                                                 !
!                                                                                                            !
! print_sol_seq=0 for NOT printing solution data files                                                       !
!                                                                                                            !
!------------------------------------------------------------------------------------------------------------!

  integer :: impt, itg, print_sol_seq=1, numseq=1
  integer :: iseq, iter_count, total_iteration
  character(30) :: char1, char2, char3, char4, char5
  character(100) :: dircommand
  real(long) :: iter_eps
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
    if (impt==1 .or. impt==2) then
      call peos_initialize_mpt(impt)
      call copy_def_peos_parameter_to_mpt(impt)
    end if
    call copy_grid_parameter_to_mpt(impt)
  end do

!stop

!
! -- Allocate arrays
  call set_allocate_size_mpt
!
  call allocate_coordinate_patch_kit_grav_mpt
  call allocate_grid_points_binary_excision
  call allocate_grid_points_asymptotic_patch
  call allocate_grid_points_binary_in_asympto
  call allocate_weight_midpoint_binary_excision
  call allocate_hgfn_bhex
  call allocate_hgfn_bhex_nb
  call allocate_hgfn_bhex_dd
  call allocate_hgfn_bhex_nd
!  call allocate_hgfn_bhex_di
  call allocate_hgfn_bhex_sd
  call allocate_BNS_CF
  call allocate_metric_on_SFC_CF
  call allocate_BNS_CF_irrot
!
  call allocate_mpatch_all_BNS_CF
  call allocate_grid_points_asymptotic_patch_mpt
  call allocate_grid_points_binary_in_asympto_mpt
  call allocate_BNS_CF_mpt
  call allocate_BNS_CF_irrot_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_peos_parameter_from_mpt(impt)
    write(6,*) "*****************************************************************************"
    write(6,*) "* Change integer input (1,2,3,4) ALSO in following subroutines:             *"
    write(6,*) "* code/Main_utility/interpolation_contour_potential_irrot_BNS_CF_3mpt.f90   *"
    write(6,*) "* spherical_data/initial_isotropic_coord/code/Import_isotr_peos_surf.f90    *"
    write(6,*) "*****************************************************************************"
    call coordinate_patch_kit_grav_grid_mpt(3)       ! 3:r_surf is used,  4:constant dr until r=7

    call calc_parameter_binary_excision
    call IO_printout_grid_data_mpt(impt)
    if (impt.ne.nmpt) then
      call calc_grid_points_binary_excision
      call calc_grid_points_binary_in_asympto(impt,nmpt)
      call copy_grid_points_binary_in_asympto_to_mpt(impt)
    end if
    call calc_weight_midpoint_binary_excision
!      call calc_weight_midpoint_binary_excision_hybrid
    if (impt==1 .or. impt==2) then
      call calc_vector_x_grav(2)
      call calc_vector_phi_grav(2)
    else
      call calc_vector_x_grav(1)
      call calc_vector_phi_grav(1)
    end if
!    call calc_vector_x_matter(2)
!    call calc_vector_phi_matter(2)
    call copy_to_mpatch_all_BNS_CF(impt)
    call copy_def_metric_and_matter_to_mpt(impt)
    call copy_def_matter_irrot_to_mpt(impt)
  end do
  call copy_from_mpatch_all_BNS_CF(nmpt)
! -- coordinates of asymptotic patch in central patches
  do impt = 1, 2
    call calc_grid_points_asymptotic_patch(impt,nmpt)
    call copy_grid_points_asymptotic_patch_to_mpt(impt)
  end do
!
! Initialization of variables
  do impt = 1, nmpt
    call copy_from_mpatch_all_BNS_CF(impt)
    call copy_def_metric_and_matter_from_mpt(impt)
    call copy_def_matter_irrot_from_mpt(impt)
    if (indata_type.eq.'3D') then
      if(NS_shape.eq.'IR')  then
        write(6,*) "########################### Reading 3D Irrotational Initial Data ####################################"
        call IO_input_initial_3D_CF_irrot_NS_mpt(impt)
!       if (impt==1 .or. impt==2)     emdc = emdini_gcm1    ! from peos_parameter_mpt file
        if (impt==1 .or. impt==2)    then
          emdc = emd(0,0,0)
          call calc_gradvep(vep,vepxf,vepyf,vepzf)
!          ======= check initial xaxis
!          write(char1, '(i5)') 0
!          char2 = adjustl(char1)
!          write(char4, '(i5)') impt
!          char5 = adjustl(char4)
!          char3 = 'iteration' // trim(char2) // '_mpt' // trim(char5) // '.txt'
!          call printout_axis_irrot_BNS_mpt(impt,char3)
        else
          emdc = 0.0d0
        end if 
      else if(NS_shape.eq.'CO')  then
        write(6,*) "########################### Reading 3D Corotating Initial Data ####################################"
        call IO_input_initial_3D_CF_NS_mpt(impt)
        if (impt==1 .or. impt==2)   then
          emdc = emd(0,0,0)
          call calc_gradvep_from_corot_id(vep,vepxf,vepyf,vepzf)
        else
          emdc = 0.0d0
        end if
      else
        write(6,*) "**************** indata_type=3D and NS_shape is neither IR nor CO....exiting"
        stop
      end if
    else
      if (impt==1 .or. impt==2)  then
        call  IO_input_initial_1D_CF_NS_mpt(impt) 
        emdc = emd(0,0,0)
        call initial_velocity_potential_NS_mpt(impt)
      else
        emdc = 0.0d0
      end if
    end if
    write (6,'(a8,i2,a16,1p,e20.12)')  "--Patch=",impt,"           emdc=",emdc
    if (impt==1 .or. impt==2)  then
!      call initial_velocity_potential_NS_mpt(impt)
!      call initial_vrot_NS_mpt(impt)
    end if
    call copy_def_metric_and_matter_to_mpt(impt)
    call copy_def_matter_irrot_to_mpt(impt)
    call copy_to_mpatch_all_BNS_CF(impt)
  end do  ! =========>>> end initialization of all patches
!

  call copy_grid_parameter_from_mpt(1)  
  if ( sw_mass_iter=='y' ) then
    call calc_irrot_BNS_CF_mpt(total_iteration)
    call calc_physical_quantities_irrot_BNS_CF_mpt
    char3 = 'main_bnsphys_all_mpt.txt'
    call printout_physq_BNS_all_mpt(char3)
    call write_last_physq_BNS_mpt
    do impt = 1, nmpt
      call copy_from_mpatch_all_BNS_CF(impt)
      call copy_def_metric_and_matter_from_mpt(impt)
      call copy_def_matter_irrot_from_mpt(impt)
      if (outdata_type.eq.'3D')  call IO_output_solution_3D_CF_irrot_NS_mpt(impt,'.las')
      if (impt==1 .or. impt==2)  call printout_NS_shape_mpt(impt)
    end do
  else
    do iseq = 1, numseq   !num_sol_seq
      write(char1, '(i5)') iseq
      char2 = adjustl(char1)
      if (iseq < 10)  then
        char3 = 'iseq0' // trim(char2)
      else
        char3 = 'iseq'  // trim(char2)
      end if
      dircommand = 'mkdir ' // char3
      call system(dircommand)
!      dircommand = 'cp -f peos_parameter_mpt* ./' // char3
!      call system(dircommand)
      call chdir(char3)
!
      write(6,'(a50,i2,a32)') '============================== Solution sequence #',iseq,' ============================== '
      write(6,'(a14,1p,2e20.12)')  "emdc COCP1,2: ", def_matter_param_real_(2,1), def_matter_param_real_(2,2)
!     Change the argument iseq to a number bigger than 1 (ex 2) if you want to  compute hydro from first iteration. 
!     Usually this is the case when the initial data are 3D.
      iter_eps = 1.0d-02
      call iter_irrot_BNS_CF_mpt(iter_count, iseq, iter_eps)
      call calc_physical_quantities_irrot_BNS_CF_mpt
      char3 = 'main_bnsphys_all_mpt.txt'
      call printout_physq_BNS_all_mpt(char3)
      call write_last_physq_BNS_mpt
      if (print_sol_seq==1) then
        do impt = 1, nmpt
          call copy_from_mpatch_all_BNS_CF(impt)
          call copy_def_metric_and_matter_from_mpt(impt)
          call copy_def_matter_irrot_from_mpt(impt)
          if (outdata_type.eq.'3D') call IO_output_solution_3D_CF_irrot_NS_mpt(impt,'.las')
          if (impt==1 .or. impt==2)  call printout_NS_shape_mpt(impt)
        end do
      end if
      call chdir('../')   ! now we are in work_area_BNS
      if ( iseq < numseq )   call next_solution_BNS_mpt(iseq)
    end do
  end if 
!
END PROGRAM Main_irrot_BNS_CF_mpt
