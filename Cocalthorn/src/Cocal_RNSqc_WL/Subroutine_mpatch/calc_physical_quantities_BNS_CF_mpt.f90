subroutine calc_physical_quantities_BNS_CF_mpt
  use phys_constant, only :  long, nmpt
  use def_binary_parameter, only : dis
  use def_matter_parameter, only : radi
  use def_matter_parameter_mpt
  use def_quantities_mpt
  use def_binary_parameter_mpt
  use grid_parameter_binary_excision 
  implicit none
  integer :: impt
!
  do impt = 1, nmpt
    call copy_from_mpatch_all_BNS_CF(impt)
    call copy_def_metric_and_matter_from_mpt(impt)
!
!!!!    call excurve
    if (impt==1 .or. impt==2)  then
      call calc_vector_x_matter(2)
      call calc_vector_phi_matter(2)
      call calc_vector_x_grav(2)
      call calc_vector_phi_grav(2)
      call excurve_CF('ns')             !   3rd order midpoint from ir0=-2,...
      call excurve_CF_gridpoint         !   4th order from ir0=-2,...

      write(6,*) '----------------------------------------- dis, radi:', dis, radi
      call calc_rest_mass_peos      ! integral over star
      call calc_mass_peos           ! integral over star only for Komar mass
      call calc_proper_mass_peos    ! integral over star
      call calc_ang_mom_peos        ! integral over star
!      call calc_ToverW_peos
      call calc_radius_CF_rsurf
!      call calc_redblue_shift
!      call calc_quad_pole_peos
      call calc_physq_center_peos_grid
      call calc_physq_cgs_peos
      call calc_enthalpy_xyzaxis
      call calc_qua_loc_spin_grav
      call calc_circ_line_peos
      call calc_circ_surf_peos
      call calc_soundspeed_peos
      call IO_printout_grid_data_mpt(impt)
      call printout_NS_shape_mpt(impt)
!      call copy_def_quantities_BNS_to_mpt(impt)
    end if
    if (impt==nmpt) then
      call calc_vector_x_grav(1)
      call calc_vector_phi_grav(1)
      call excurve_CF('bh')             !   3rd order midpoint from ir0=0,...
      call excurve_CF_gridpoint_bhex    !   4th order from ir0=0,...

!     copy radi (i=5) of COCP1 to ARCP
      def_matter_param_real_(5,nmpt) = def_matter_param_real_(5,1)
      radi = def_matter_param_real_(5,nmpt)
      write(6,*) '----------------------------------------- dis, radi:', dis, radi
      call calc_mass_asympto('ns')
      call calc_angmom_asympto('ns')
      call calc_admmom_asympto('ns')
      def_matter_param_real_(5,nmpt) = 0.0d0
      radi = def_matter_param_real_(5,nmpt)
!      call copy_def_quantities_BNS_to_mpt(impt)
    end if
    call copy_to_mpatch_all_BNS_CF(impt)
    write(6,*) "CCCC resetmass, qc:", def_quantities_real_( 4,impt), def_quantities_real_(49,impt)
    write(6,*) "DDDD Madm asym, qc:", def_quantities_real_( 7,impt), def_quantities_real_(49,impt)
    write(6,*) "EEEE dis, qc      :", dis_(impt), def_quantities_real_(49,impt)
  end do
!
end subroutine calc_physical_quantities_BNS_CF_mpt
