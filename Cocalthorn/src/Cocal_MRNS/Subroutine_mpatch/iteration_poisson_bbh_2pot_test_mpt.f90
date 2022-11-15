subroutine iteration_poisson_bbh_2pot_test_mpt(iter_count)
  use phys_constant, only :  long, nnrg, nntg, nnpg, nmpt
  use grid_parameter
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use weight_midpoint_grav
  use def_metric, only : psi, alph
  use def_metric_mpt, only : psi_, alph_
  use copy_array_4dto3d_mpt
  use make_array_2d
  use make_array_3d
  use interface_poisson_solver
  use interface_update_grfield
  use interface_error_metric
  use interface_error_metric_type0
  use interface_interpo_fl2gr
  use interface_sourceterm_poisson_solver_test
  use interface_sourceterm_exsurf_eqm_binary
  use interface_sourceterm_surface_int
!  use interface_sourceterm_volume_int_bbh_2pot_test
  use interface_sourceterm_2pot_test
  use interface_poisson_solver_binary_bhex_homosol
  use interface_bh_boundary_test_mpt
  use interface_bh_boundary_psi_test_mpt
  use interface_bh_boundary_alph_test_mpt
  use interface_bh_boundary_nh_psi_test
  use interface_bh_boundary_nh_alph_test
  use interface_interpolation_fillup_binary_mpt
!  use interface_poisson_solver_binary_bhex
  use interface_bh_boundary_test
  implicit none
  real(long), pointer :: sou(:,:,:)
  real(long), pointer :: pot(:,:,:), psi_bak(:,:,:), alph_bak(:,:,:)
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long), pointer :: fnc_itp(:,:,:)
  real(long) :: error_psi =0.0d0, error_alph =0.0d0, error_tmp =0.0d0, count
  integer    :: iter_count, flag_psi = 0, flag_alph = 0, flag_tmp = 0
  integer    :: irf, itf, ipf, irg, itg, ipg, ihy, impt, impt_ex
!
  call set_allocate_size_mpt
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call alloc_array3d(psi_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(alph_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(fnc_itp,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_exsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_exsurf,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
!
  iter_count = 0
  do
    iter_count = iter_count + 1      
    count = dble(iter_count) 
!
    do impt = 1, nmpt
      call copy_grid_parameter_from_mpt(impt)
      conv_gra = dmin1(conv0_gra,conv_ini*count)
      conv_den = dmin1(conv0_den,conv_ini*count)
      call copy_grid_parameter_to_mpt(impt)
!
      call copy_from_mpatch_all_test(impt)
      call calc_vector_x_grav(0)
      call calc_vector_x_matter(0)
      call calc_vector_phi_grav(0)
      call calc_vector_phi_matter(0)
      call copy_to_mpatch_all_test(impt)
! --
! For psi
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
!
      if(impt.eq.1) impt_ex = 2
      if(impt.eq.2) impt_ex = 1
      call copy_from_mpatch_all_test(impt_ex)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt_ex)
      call sourceterm_exsurf_eqm_binary(psi,sou_exsurf,dsou_exsurf)
      call copy_from_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt)
! --
!      call reset_bh_boundary('n') ! 'd' for Dirichlet
!      call sourceterm_surface_int(psi,0,sou_bhsurf,dsou_bhsurf)
      call sourceterm_surface_int(psi,nrg,sou_outsurf,dsou_outsurf)
!
!-- For BH boundary substitute Dirichlet data to  sou_bhsurf
!--                            Neumann   data to dsou_bhsurf
!      call bh_boundary_psi_test_mpt(impt,'n',sou_bhsurf,dsou_bhsurf)
!      call poisson_solver_binary_bhex_homosol('nd',sou, &
      call bh_boundary_psi_test_mpt(impt,'n',sou_bhsurf,dsou_bhsurf)
      call poisson_solver_binary_bhex_homosol('nd',sou, &
      &                                        sou_exsurf,dsou_exsurf, &
      &                                        sou_bhsurf,dsou_bhsurf, &
      &                                       sou_outsurf,dsou_outsurf,pot)
!
      psi_bak(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,psi,conv_gra)
!
      call copy_grid_parameter_interpo_from_mpt(impt_ex)
      call copy_grid_parameter_binary_excision_interpo_from_mpt(impt_ex)
      call copy_coordinate_grav_extended_interpo_from_mpt(impt_ex)
      call copy_array4dto3d_mpt(impt_ex,psi_,fnc_itp, &
      &                         0,nrg_itp,0,ntg_itp,0,npg_itp)
      call interpolation_fillup_binary_mpt(psi,fnc_itp)
!
      call error_metric_type0(psi,psi_bak,error_tmp,flag_tmp,'bh')
      flag_psi = max(flag_psi,flag_tmp)
      error_psi = dmax1(error_psi,error_tmp)
      call printout_error_metric(iter_count,error_psi)
      call reset_bh_boundary('n') ! Dirichlet at rgout
      call copy_to_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_to_mpt(impt)
!
! ---
! For alpha
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
!
      if(impt.eq.1) impt_ex = 2
      if(impt.eq.2) impt_ex = 1
      call copy_from_mpatch_all_test(impt_ex)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt_ex)
      call sourceterm_exsurf_eqm_binary(alph,sou_exsurf,dsou_exsurf)
      call copy_from_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt)
!
!-- For BH boundary substitute Dirichlet data to  sou_bhsurf
!--                            Neumann   data to dsou_bhsurf
!      call reset_bh_boundary('n') ! 'd' for Dirichlet
!      call sourceterm_surface_int(alph,0,sou_bhsurf,dsou_bhsurf)
      call sourceterm_surface_int(alph,nrg,sou_outsurf,dsou_outsurf)
!      call bh_boundary_alph_test_mpt(impt,'n',sou_bhsurf,dsou_bhsurf)
      call bh_boundary_alph_test_mpt(impt,'n',sou_bhsurf,dsou_bhsurf)
!      call sourceterm_volume_int_bbh_2pot_test(sou)
      call sourceterm_2pot_test(sou)
!      call poisson_solver_binary_bhex_homosol('nd',sou, &
      call poisson_solver_binary_bhex_homosol('nd',sou, &
      &                                       sou_exsurf,dsou_exsurf, &
      &                                       sou_bhsurf,dsou_bhsurf, & 
      &                                       sou_outsurf,dsou_outsurf,pot)
      alph_bak(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,alph,conv_gra)
!
      call copy_grid_parameter_interpo_from_mpt(impt_ex)
      call copy_grid_parameter_binary_excision_interpo_from_mpt(impt_ex)
      call copy_coordinate_grav_extended_interpo_from_mpt(impt_ex)
      call copy_array4dto3d_mpt(impt_ex,alph_,fnc_itp, &
      &                         0,nrg_itp,0,ntg_itp,0,npg_itp)
      call interpolation_fillup_binary_mpt(alph,fnc_itp)
!
      call error_metric_type0(alph,alph_bak,error_tmp,flag_tmp,'bh')
      flag_alph = max(flag_alph,flag_tmp)
      error_alph = dmax1(error_alph,error_tmp)
      call printout_error_metric(iter_count,error_alph)
      call reset_bh_boundary('n') ! Dirichlet at rgout
      call copy_to_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_to_mpt(impt)
    end do
!
    if ((flag_psi==0).and.(flag_alph==0)) exit
    if (iter_count >= iter_max) exit
    if (iter_count >= 150 .and. error_psi > 1.5d0) exit
    if (iter_count >= 150 .and. error_alph > 1.5d0) exit
    flag_psi = 0
    flag_alph = 0
    error_psi  = 0.0d0
    error_alph = 0.0d0
  end do
!
  deallocate(sou)
  deallocate(pot)
  deallocate(psi_bak)
  deallocate(alph_bak)
  deallocate(fnc_itp)
  deallocate(sou_exsurf)
  deallocate(dsou_exsurf)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine iteration_poisson_bbh_2pot_test_mpt
