subroutine iteration_poisson_BHNS_test_mpt(iter_count)
  use phys_constant, only :  long, nnrg, nntg, nnpg, nmpt
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use weight_midpoint_grav
  use def_metric, only : psi
  use def_matter, only : emd, emdg
  use make_array_2d
  use make_array_3d
  use radial_green_fn_grav
  use radial_green_fn_grav_bhex_di
  use radial_green_fn_grav_bhex_nb
  use radial_green_fn_grav_bhex_nh
  use radial_green_fn_grav_bhex_dh
  use interface_poisson_solver
  use interface_update_grfield
  use interface_error_metric
  use interface_interpo_fl2gr
  use interface_sourceterm_poisson_solver_test
  use interface_sourceterm_exsurf_eqm_binary
  use interface_copy_to_hgfn_and_gfnsf
  use interface_error_metric_type0
  use interface_poisson_solver_binary
  use interface_sourceterm_surface_int
  use interface_poisson_solver_binary_bhex_homosol
!  use interface_poisson_solver_binary_bhex
  use interface_bh_boundary_BHNS_test_mpt
!  
  implicit none
  real(long), pointer :: sou(:,:,:), pot(:,:,:), psi_bak(:,:,:)
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long) :: error_psi, count
  integer    :: iter_count, flag = 0
  integer    :: irf, itf, ipf, irg, itg, ipg, ihy, impt, impt_ex
!
  call set_allocate_size_mpt
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call alloc_array3d(psi_bak,0,nrg,0,ntg,0,npg)
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
      if(impt.eq.1) impt_ex = 2
      if(impt.eq.2) impt_ex = 1
      call copy_from_mpatch_all_test(impt_ex)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt_ex)
      call sourceterm_exsurf_eqm_binary(psi,sou_exsurf,dsou_exsurf)
      call copy_from_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt)
!
      if(impt.eq.1) then
        sou(0:nrg,0:ntg,0:npg) = 0.0d0
        call reset_bh_boundary('d') ! 'd' for Dirichlet
        call sourceterm_surface_int(psi,0,sou_bhsurf,dsou_bhsurf)
        call sourceterm_surface_int(psi,nrg,sou_outsurf,dsou_outsurf)
        call bh_boundary_BHNS_test_mpt('dh',sou_exsurf,dsou_exsurf)
        call poisson_solver_binary_bhex_homosol('dh',sou, &
        &                                     sou_exsurf,dsou_exsurf, &
        &                                     sou_bhsurf,dsou_bhsurf, & 
        &                                    sou_outsurf,dsou_outsurf,pot)
! --
      end if
!
      if(impt.eq.2) then
        emdg = 0.0d0
        call interpo_fl2gr(emd, emdg)
        call sourceterm_poisson_solver_test(sou)
        call calc_hgfn_bhex_nb
        call copy_to_hgfn_and_gfnsf(hgfn_nb,gfnsf_nb)
        call poisson_solver_binary(sou,sou_exsurf,dsou_exsurf,pot)
      end if
!
      psi_bak(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,psi,conv_gra)
      call error_metric_type0(psi,psi_bak,error_psi,flag,'bh')
      call printout_error_metric(iter_count,error_psi)
      if(impt.eq.1) call reset_bh_boundary('d') ! 'd' for Dirichlet
      call copy_to_mpatch_all_test(impt)
      call copy_metric_and_matter_BHNS_test_to_mpt(impt)
    end do
!
    if (flag == 0) exit
    if (iter_count >= iter_max) exit
    if (iter_count >= 10 .and. error_psi > 1.5d0) exit
    flag = 0
  end do
!
  deallocate(sou)
  deallocate(pot)
  deallocate(psi_bak)
  deallocate(sou_exsurf)
  deallocate(dsou_exsurf)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine iteration_poisson_BHNS_test_mpt
