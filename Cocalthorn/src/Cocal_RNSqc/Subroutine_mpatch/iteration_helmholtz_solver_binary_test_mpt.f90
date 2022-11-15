subroutine iteration_helmholtz_solver_binary_test_mpt(iter_count)
  use phys_constant, only :  long, nnrg, nntg, nnpg, nmpt
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use weight_midpoint_grav
  use def_metric, only : psi
  use def_matter, only : emd, emdg
  use def_binary_parameter, only : mass_ratio
  use make_array_2d
  use make_array_3d
  use interface_update_grfield
  use interface_error_metric_type0
  use interface_interpo_fl2gr
  use interface_sourceterm_helmholtz_solver_test
  use interface_sourceterm_exsurf_binary_COCP
  use interface_sourceterm_outsurf_COCP_from_ARCP
  use interface_sourceterm_insurf_ARCP_from_COCP
  use interface_helmholtz_solver_binary
  use interface_helmholtz_solver_asymptotic_patch_homosol
  use interface_helmholtz_solver_asymptotic_patch_homosol_outgoing
  use interface_poisson_solver_binary_star_homosol
  implicit none
  real(long), pointer :: sou(:,:,:), pot(:,:,:), psi_bak(:,:,:)
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
  real(long), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long) :: error_psi, error_tmp, count
  integer    :: iter_count, iter_extra = 0, istep_niq = -1, &
  &             flag = 0, flag0 = 0, flag_param = 99
  integer    :: impt, impt_bin, impt_ex
  integer    ::  irg, itg, ipg, ihy
  character(len=2) :: chgreen, chpa, chpB
  character(len=4) :: chbe
!
  call set_allocate_size_mpt
!
  call alloc_array2d(sou_insurf,0,ntg,0,npg)
  call alloc_array2d(dsou_insurf,0,ntg,0,npg)
  call alloc_array2d(sou_exsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_exsurf,0,ntg,0,npg)
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
  call alloc_array3d(psi_bak,0,nrg,0,ntg,0,npg)
! 
  iter_count = 0
  do
    iter_count = iter_count + 1      
    count = dble(iter_count) 
!
    do impt = 1, nmpt
      write(6,'(a10,i2)')" Patch # =", impt
      call copy_grid_parameter_from_mpt(impt)
      conv_gra = dmin1(conv0_gra,conv_ini*count)
      conv_den = dmin1(conv0_den,conv_ini*count)
      call copy_grid_parameter_to_mpt(impt)
!
      call copy_from_mpatch_all_BBH_CF(impt)
      call copy_metric_and_matter_BHNS_test_from_mpt(impt)
      call calc_vector_x_grav(2)
      call calc_vector_x_matter(2)
      call calc_vector_phi_grav(2)
      call calc_vector_phi_matter(2)
      call copy_to_mpatch_all_BBH_CF(impt)
! --
!
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_insurf(0:ntg,0:npg) = 0.0d0 ; dsou_insurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
!
      if (impt.ne.nmpt) then
!!        call sourceterm_helmholtz_solver_test(sou,'he')
        call sourceterm_helmholtz_solver_test(sou,'po')
        call sourceterm_exsurf_binary_COCP(impt,'psi ','ev',sou_exsurf, &
        &                                                  dsou_exsurf)
        call sourceterm_outsurf_COCP_from_ARCP(impt,'psi ','ev',sou_outsurf, &
        &                                                      dsou_outsurf)
        chgreen = 'sd'
        call poisson_solver_binary_star_homosol(chgreen,sou, &
        &                                       sou_exsurf,dsou_exsurf, &
        &                                       sou_outsurf,dsou_outsurf,pot)
!!        call helmholtz_solver_binary(psi,sou,sou_exsurf,dsou_exsurf, &
!!        &                                    sou_outsurf,dsou_outsurf,pot)
      else if (impt.eq.nmpt) then
        call sourceterm_insurf_ARCP_from_COCP(impt,'psi ','ev',sou_insurf, &
        &                                                     dsou_insurf)
!hrethadv        call helmholtz_solver_asymptotic_patch_homosol &
!hrethadv        &              (sou,sou_insurf,dsou_insurf,pot)
        call helmholtz_solver_asymptotic_patch_homosol_outgoing &
        &                       (sou,sou_insurf,dsou_insurf,pot)
      end if
! 
! --
      psi_bak(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,psi,conv_gra)
      if (impt.eq.1) call error_metric_type0(psi,psi_bak,error_psi,flag0,'ns')
      if (impt.eq.2) call error_metric_type0(psi,psi_bak,error_psi,flag0,'ns')
      if (impt.eq.3) call error_metric_type0(psi,psi_bak,error_psi,flag0,'bh')
      flag = flag + flag0
      error_tmp = dmax1(error_psi,error_tmp)
      call printout_error_metric(iter_count,error_psi)
      call copy_to_mpatch_all_BBH_CF(impt)
      call copy_metric_and_matter_BHNS_test_to_mpt(impt)
    end do
!
! -- For ciruclar solutions
    if (mass_ratio.ne.1.0d0.and.error_tmp.le.eps) then
!         error_tmp.le.1.0d-04) then
!         error_tmp.le.1.0d-04.and.iter_extra.ge.20) then
      iter_extra = 0
      call calc_physical_quantities_helm_test_mpt
      call adjust_multi_parameter_helm_test_mpt(flag_param,istep_niq)
      flag = flag + flag_param
      call update_coordinates_mpt
    end if
!
! -- convergence of iteration is checked.
!
    if (flag == 0) exit
    if (iter_count >= iter_max) exit
!    if (iter_count >= 20 .and. error_psi > 1.5d0) exit
    flag  = 0
    flag0 = 0
    flag_param = 0
    error_psi = 0.0d0
    error_tmp = 0.0d0
  end do
!
  deallocate(sou)
  deallocate(sou_exsurf)
  deallocate(dsou_exsurf)
  deallocate(sou_insurf)
  deallocate(dsou_insurf)
  deallocate(pot)
  deallocate(psi_bak)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine iteration_helmholtz_solver_binary_test_mpt
