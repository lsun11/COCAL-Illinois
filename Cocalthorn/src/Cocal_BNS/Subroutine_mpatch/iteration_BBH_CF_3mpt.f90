subroutine iteration_BBH_CF_3mpt(iter_count)
  use phys_constant, only :  long, nmpt
  use grid_parameter
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use def_metric, only : psi, alph, alps, bvxd, bvyd, bvzd
  use def_bh_parameter, only : bh_bctype, bh_soltype, ome_bh
  use def_binary_parameter, only : dis
  use copy_array_4dto3d_mpt
  use make_array_2d
  use make_array_3d
  use make_array_4d
  use interface_compute_alps2alph
  use interface_sourceterm_HaC_CF
  use interface_sourceterm_trG_CF
  use interface_sourceterm_MoC_CF_with_divshift
  use interface_sourceterm_surface_int
  use interface_sourceterm_exsurf_binary_COCP
  use interface_sourceterm_outsurf_COCP_from_ARCP
  use interface_sourceterm_insurf_ARCP_from_COCP
  use interface_bh_boundary_CF
  use interface_outer_boundary_d_BBH_CF
  use interface_outer_boundary_d_psi
  use interface_outer_boundary_d_alps
  use interface_outer_boundary_d_bvxd
  use interface_outer_boundary_d_bvyd
  use interface_outer_boundary_d_bvzd
  use interface_poisson_solver_binary_bhex_homosol
  use interface_poisson_solver_asymptotic_patch_homosol
  use interface_interpolation_fillup_binary_COCP
  use interface_update_grfield
  use interface_error_metric
  use interface_error_metric_type0
  use interface_error_metric_type2
  implicit none
  real(long), pointer :: sou(:,:,:), pot(:,:,:)
  real(long), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:) 
  real(long), pointer :: pot_bak(:,:,:), work(:,:,:)
  real(long), pointer :: potx_bak(:,:,:), poty_bak(:,:,:), potz_bak(:,:,:) 
  real(long), pointer :: souvec(:,:,:,:)
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long), pointer :: fnc_itp(:,:,:)

  real(long) :: count
  real(long) :: error_psi  = 0.0d0, error_alph = 0.0d0, error_tmp = 0.0d0, &
  &             error_bvxd = 0.0d0, error_bvyd = 0.0d0, error_bvzd = 0.0d0, &
  &             error_ome = 2.0d0
  real(long) :: pari, ome_bh_new
  integer    :: iter_count, iter_extra = 0, istep_niq = -1, &
  &             flag_tmp = 0, flag_param = 99, &
  &             flag_psi = 0, flag_alph = 0, &
  &            flag_bvxd = 0, flag_bvyd = 0, flag_bvzd = 0
  integer    :: irg, itg, ipg, ii
  integer    :: impt, impt_ex, impt_bin
  character(len=2) :: chgreen, chpa, chpB
  character(len=4) :: chbe
!
  call set_allocate_size_mpt
  call alloc_array4d(souvec,0,nrg,0,ntg,0,npg,1,3)
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potx,0,nrg,0,ntg,0,npg)
  call alloc_array3d(poty,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potz,0,nrg,0,ntg,0,npg)
  call alloc_array3d(pot_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potx_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(poty_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(potz_bak,0,nrg,0,ntg,0,npg)
  call alloc_array3d(work,0,nrg,0,ntg,0,npg)
  call alloc_array3d(fnc_itp,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_exsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_exsurf,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_insurf,0,ntg,0,npg)
  call alloc_array2d(dsou_insurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
!
  iter_count = 0
!
  do
    iter_count = iter_count + 1
    iter_extra = iter_extra + 1
    count = dble(iter_count)
!
    do impt = 1, nmpt
      write(6,'(a10,i2)')" Patch # =", impt
      call copy_grid_parameter_from_mpt(impt)
      conv_gra = dmin1(conv0_gra,conv_ini*count)
      conv_den = dmin1(conv0_den,conv_ini*count)
      call copy_grid_parameter_to_mpt(impt)
! ---
      call copy_from_mpatch_all_BBH_CF(impt)
      call copy_def_metric_from_mpt(impt)
      call calc_vector_x_grav(2)
      call calc_vector_phi_grav(2)
      call calc_vector_bh(2)
      call excurve
      call excurve_CF_gridpoint_bhex
      if (impt.eq.nmpt) call calc_mass_asympto('bh')
      call copy_to_mpatch_all_BBH_CF(impt)
!
! --- psi
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
      call sourceterm_HaC_CF(sou)
!
      if (impt.ne.nmpt) then
        call sourceterm_exsurf_binary_COCP(impt,'psi ','ev',sou_exsurf, &
        &                                                  dsou_exsurf)
        call sourceterm_outsurf_COCP_from_ARCP(impt,'psi ','ev',sou_outsurf, &
        &                                                      dsou_outsurf)
        call bh_boundary_CF(sou_bhsurf,dsou_bhsurf,'psi ')
        if (bh_bctype.eq.'TU') chgreen = 'dd'
        if (bh_bctype.eq.'AH'.or.bh_bctype.eq.'PG') chgreen = 'nd'
        call poisson_solver_binary_bhex_homosol(chgreen,sou, &
        &                                       sou_exsurf,dsou_exsurf, &
        &                                       sou_bhsurf,dsou_bhsurf, &
        &                                       sou_outsurf,dsou_outsurf,pot)
      else if (impt.eq.nmpt) then
        call sourceterm_insurf_ARCP_from_COCP(impt,'psi ','ev',sou_insurf, &
        &                                                     dsou_insurf)
!        call outer_boundary_d_psi(sou_outsurf)
!!sou = 0.0d0
!!dsou_outsurf = 0.0d0
        call outer_boundary_d_BBH_CF(sou_outsurf,'psi ')
        call poisson_solver_asymptotic_patch_homosol('dd',sou, &
        &                                      sou_insurf,dsou_insurf, &
        &                                     sou_outsurf,dsou_outsurf,pot)
      end if
!
      pot_bak(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,psi,conv_gra)
      call reset_outer_boundary_BBH_CF('psi ')
!
      if (impt.ne.nmpt) call interpolation_fillup_binary_COCP &
      &                                (impt,'psi ','ev',psi)
      call error_metric_type0(psi,pot_bak,error_tmp,flag_tmp,'bh')
       flag_psi =  max ( flag_psi, flag_tmp)
      error_psi = dmax1(error_psi,error_tmp)
      call printout_error_metric(iter_count,error_psi)
!      if (impt.eq.nmpt) call reset_bh_boundary_AH !mpt version not available
      call copy_to_mpatch_all_BBH_CF(impt)
      call copy_def_metric_to_mpt(impt)
!
! --- alps
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
      alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)
      call sourceterm_trG_CF(sou)
!
      if (impt.ne.nmpt) then
        call sourceterm_exsurf_binary_COCP(impt,'alps','ev',sou_exsurf, &
        &                                                  dsou_exsurf)
        call sourceterm_outsurf_COCP_from_ARCP(impt,'alps','ev',sou_outsurf, &
        &                                                      dsou_outsurf)
        call bh_boundary_CF(sou_bhsurf,dsou_bhsurf,'alps')
        call poisson_solver_binary_bhex_homosol('dd',sou, &
        &                                       sou_exsurf,dsou_exsurf, &
        &                                       sou_bhsurf,dsou_bhsurf, &
        &                                       sou_outsurf,dsou_outsurf,pot)
      else if (impt.eq.nmpt) then
        call sourceterm_insurf_ARCP_from_COCP(impt,'alps','ev',sou_insurf, &
        &                                                     dsou_insurf)
!        call outer_boundary_d_alps(sou_outsurf)
        call outer_boundary_d_BBH_CF(sou_outsurf,'alps')
        call poisson_solver_asymptotic_patch_homosol('dd',sou, &
        &                                      sou_insurf,dsou_insurf, &
        &                                     sou_outsurf,dsou_outsurf,pot)
      end if
!
      call compute_alps2alph(pot,psi)
      pot_bak(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,alph,conv_gra)
      call reset_outer_boundary_BBH_CF('alph')
      if (impt.ne.nmpt) call interpolation_fillup_binary_COCP &
      &                                (impt,'alph','ev',alph)
      call error_metric_type0(alph,pot_bak,error_tmp,flag_tmp,'bh')
       flag_alph =  max ( flag_alph, flag_tmp)
      error_alph = dmax1(error_alph,error_tmp)
      call printout_error_metric(iter_count,error_alph)
!      if (impt.eq.nmpt) call reset_bh_boundary_AH !mpt version not available
      alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)
      call copy_to_mpatch_all_BBH_CF(impt)
      call copy_def_metric_to_mpt(impt)
! 
! --- shift
      souvec(0:nrg,0:ntg,0:npg,1:3)  = 0.0d0
      call sourceterm_MoC_CF_with_divshift(souvec)
!
      do ii = 1, 3
        sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
        sou_bhsurf(0:ntg,0:npg) = 0.0d0 ; dsou_bhsurf(0:ntg,0:npg) = 0.0d0
        sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
        sou(0:nrg,0:ntg,0:npg) = souvec(0:nrg,0:ntg,0:npg,ii)
        if (ii.eq.1) work(0:nrg,0:ntg,0:npg) = bvxd(0:nrg,0:ntg,0:npg)
        if (ii.eq.2) work(0:nrg,0:ntg,0:npg) = bvyd(0:nrg,0:ntg,0:npg)
        if (ii.eq.3) work(0:nrg,0:ntg,0:npg) = bvzd(0:nrg,0:ntg,0:npg)
        if (ii.eq.1) then; chbe = 'bvxd'; chpa = 'od'; end if
        if (ii.eq.2) then; chbe = 'bvyd'; chpa = 'od'; end if
        if (ii.eq.3) then; chbe = 'bvzd'; chpa = 'ev'; end if
!
        if (impt.ne.nmpt) then
          chpB = 'ev'
          if (impt.eq.2.and.ii.eq.1) chpB = 'od'
          if (impt.eq.2.and.ii.eq.2) chpB = 'od'
          call sourceterm_exsurf_binary_COCP(impt,chbe,chpa,sou_exsurf, &
          &                                                dsou_exsurf)
          call sourceterm_outsurf_COCP_from_ARCP(impt,chbe,chpB,sou_outsurf,&
          &                                                    dsou_outsurf)
          call bh_boundary_CF(sou_bhsurf,dsou_bhsurf,chbe)
          call poisson_solver_binary_bhex_homosol('dd',sou, &
          &                                     sou_exsurf,dsou_exsurf, &
          &                                     sou_bhsurf,dsou_bhsurf, &
          &                                    sou_outsurf,dsou_outsurf,pot)
        else if (impt.eq.nmpt) then
          call sourceterm_insurf_ARCP_from_COCP(impt,chbe,chpa,sou_insurf, &
          &                                                   dsou_insurf)
!          if (ii.eq.1) call outer_boundary_d_bvxd(sou_outsurf)
!          if (ii.eq.2) call outer_boundary_d_bvyd(sou_outsurf)
!          if (ii.eq.3) call outer_boundary_d_bvzd(sou_outsurf)
          if (ii.eq.1) call outer_boundary_d_BBH_CF(sou_outsurf,'bvxd')
          if (ii.eq.2) call outer_boundary_d_BBH_CF(sou_outsurf,'bvyd')
          if (ii.eq.3) call outer_boundary_d_BBH_CF(sou_outsurf,'bvzd')
!!!          call poisson_solver_asymptotic_patch_homosol('di',sou, &
!!!          call poisson_solver_asymptotic_patch_homosol('nd',sou, &
          call poisson_solver_asymptotic_patch_homosol('dd',sou, &
          &                                      sou_insurf,dsou_insurf, &
          &                                     sou_outsurf,dsou_outsurf,pot)
        end if
        if (ii.eq.1) potx(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
        if (ii.eq.2) poty(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
        if (ii.eq.3) potz(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
      end do
! 
      potx_bak(0:nrg,0:ntg,0:npg) = bvxd(0:nrg,0:ntg,0:npg)
      poty_bak(0:nrg,0:ntg,0:npg) = bvyd(0:nrg,0:ntg,0:npg)
      potz_bak(0:nrg,0:ntg,0:npg) = bvzd(0:nrg,0:ntg,0:npg)
      call update_grfield(potx,bvxd,conv_gra)
      call update_grfield(poty,bvyd,conv_gra)
      call update_grfield(potz,bvzd,conv_gra)
!
      if (impt.ne.nmpt) then 
        call interpolation_fillup_binary_COCP(impt,'bvxd','od',bvxd)
        call interpolation_fillup_binary_COCP(impt,'bvyd','od',bvyd)
        call interpolation_fillup_binary_COCP(impt,'bvzd','ev',bvzd)
      end if
!      if (bh_bctype.eq.'TU') call reset_metric_CF
!      call reset_metric_CF
      call error_metric_type2(bvxd,potx_bak,error_tmp,flag_tmp,'bh')
       flag_bvxd =  max ( flag_bvxd, flag_tmp)
      error_bvxd = dmax1(error_bvxd,error_tmp)
      call error_metric_type2(bvyd,poty_bak,error_tmp,flag_tmp,'bh')
       flag_bvyd =  max ( flag_bvyd, flag_tmp)
      error_bvyd = dmax1(error_bvyd,error_tmp)
      call error_metric_type2(bvzd,potz_bak,error_tmp,flag_tmp,'bh')
       flag_bvzd =  max ( flag_bvzd, flag_tmp)
      error_bvzd = dmax1(error_bvzd,error_tmp)
      call printout_error_metric(iter_count,error_bvxd)
      call printout_error_metric(iter_count,error_bvyd)
      call printout_error_metric(iter_count,error_bvzd)
!      if (impt.eq.nmpt) call reset_bh_boundary_AH !mpt version not available
      call copy_to_mpatch_all_BBH_CF(impt)
      call copy_def_metric_to_mpt(impt)
!
! ---
    end do
!
! -- For ciruclar solutions
    error_tmp = dmax1(error_psi,error_alph,error_bvxd, &
    &                           error_bvyd,error_bvzd)
    if ((bh_soltype.eq.'CI'.or.bh_soltype.eq.'SQ').and. &
         error_tmp.le.eps) then 
!         error_tmp.le.1.0d-04) then 
!         error_tmp.le.1.0d-04.and.iter_extra.ge.20) then 
      iter_extra = 0
      call calc_physical_quantities_BBH_CF_mpt
      call adjust_multi_parameter_ome_cm_ratio_mpt(flag_param,istep_niq)
      call update_coordinates_mpt
    end if
!
! -- convergence of iteration is checked.
!
    if ((flag_param==0).and.(flag_psi==0).and.(flag_alph==0).and. &
    &   (flag_bvxd==0).and.(flag_bvyd==0).and.(flag_bvzd==0)) exit
    if (iter_count >= iter_max) exit
!
! --
!
    flag_psi = 0
    flag_alph = 0
    flag_bvxd = 0
    flag_bvyd = 0
    flag_bvzd = 0
    flag_param= 0
    error_psi  = 0.0d0
    error_alph = 0.0d0
    error_bvxd = 0.0d0
    error_bvyd = 0.0d0
    error_bvzd = 0.0d0
  end do
!
  deallocate(souvec)
  deallocate(sou)
  deallocate(pot)
  deallocate(potx)
  deallocate(poty)
  deallocate(potz)
  deallocate(pot_bak)
  deallocate(potx_bak)
  deallocate(poty_bak)
  deallocate(potz_bak)
  deallocate(work)
  deallocate(fnc_itp)
  deallocate(sou_exsurf)
  deallocate(dsou_exsurf)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_insurf)
  deallocate(dsou_insurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine iteration_BBH_CF_3mpt
