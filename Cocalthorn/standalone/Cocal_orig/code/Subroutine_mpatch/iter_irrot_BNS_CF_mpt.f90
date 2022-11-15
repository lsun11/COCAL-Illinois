subroutine iter_irrot_BNS_CF_mpt(iter_count, iseq, iter_eps)
  use phys_constant, only :  long, nmpt
  use grid_parameter
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp
  use def_metric, only : psi, alph, alps, bvxd, bvyd, bvzd
  use coordinate_grav_r
  use grid_points_binary_excision, only : irg_exin, irg_exout
!
  use def_matter
  use def_velocity_potential
  use def_matter_parameter
  use interface_interpo_fl2gr
  use interface_sourceterm_HaC_CF_peos_irrot
  use interface_sourceterm_HaC_CF
  use interface_sourceterm_trG_CF_peos_irrot
  use interface_sourceterm_trG_CF_helical
  use interface_sourceterm_MoC_CF_with_divshift_peos_irrot
  use interface_sourceterm_MoC_CF_with_divshift_r3rd
!  use interface_violation_midpoint_MoC_CF_peos_irrot
!  use interface_violation_gridpoint_MoC_CF_peos_irrot
!  use interface_violation_midpoint_HaC_CF_peos_irrot
!  use interface_violation_gridpoint_HaC_CF_peos_irrot
  use interface_hydrostatic_eq_CF_peos_irrot
  use interface_calc_gradvep
  use interface_hydro_irbns_vep_CF_peos
  use interface_calc_surface
  use interface_update_matter
  use interface_update_parameter_BNS_irrot
  use interface_update_surface
  use interface_error_matter
  use interface_error_matter_type2
  use interface_error_velocity_potential
!  use interface_IO_output_3D_general
!  use interface_IO_output_1D_general
!
  use def_binary_parameter, only : dis
  use copy_array_4dto3d_mpt
  use make_array_2d
  use make_array_3d
  use make_array_4d
  use interface_compute_alps2alph
  use interface_sourceterm_surface_int
  use interface_sourceterm_exsurf_binary_COCP
  use interface_sourceterm_outsurf_COCP_from_ARCP
  use interface_sourceterm_insurf_ARCP_from_COCP
  use interface_outer_boundary_d_BBH_CF
  use interface_outer_boundary_d_psi
  use interface_outer_boundary_d_alps
  use interface_outer_boundary_d_bvxd
  use interface_outer_boundary_d_bvyd
  use interface_outer_boundary_d_bvzd
!  use interface_poisson_solver_binary_bhex_homosol
!  use interface_poisson_solver_binary_star_homosol_plmex
  use interface_poisson_solver_binary_star_homosol
  use interface_poisson_solver_asymptotic_patch_homosol
  use interface_interpolation_fillup_binary_COCP
  use interface_update_grfield
  use interface_error_metric
  use interface_error_metric_type0
  use interface_error_metric_type2
  use interface_error_metric_type2_mpt
  use interface_IO_output_1D_general
  implicit none
  real(long), pointer :: sou(:,:,:), pot(:,:,:), HaC_vio(:,:,:)
  real(long), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:) 
  real(long), pointer :: pot_bak(:,:,:), work(:,:,:)
  real(long), pointer :: potx_bak(:,:,:), poty_bak(:,:,:), potz_bak(:,:,:) 
  real(long), pointer :: souvec(:,:,:,:), MoC_vio(:,:,:,:)
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
!  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long), pointer :: fnc_itp(:,:,:)
  real(long), pointer :: potf(:,:,:), emd_bak(:,:,:), vep_bak(:,:,:)
  real(long), pointer :: potrs(:,:)
  real(long), pointer :: potut(:,:,:), potux(:,:,:), potuy(:,:,:), potuz(:,:,:)

  real(long) :: count, iter_eps
  real(long) :: error_psi  = 0.0d0, error_alph = 0.0d0, error_tmp  = 0.0d0, &
  &             error_bvxd = 0.0d0, error_bvyd = 0.0d0, error_bvzd = 0.0d0, &
  &             error_ome  = 2.0d0, error_emd  = 0.0d0, error_metr = 0.0d0, &
  &             error_vep  = 0.0d0, error_emd_last, error_psi_last, error_alph_last, &
  &             error_vep_last
  real(long) :: pari, ome_bh_new, work_shift, error_metric_mpt_mid(10,10)
  integer    :: iter_count, iter_extra = 0, istep_niq = -1,     &
  &             flag_tmp  = 0, flag_param = 99, flag_emd  = 0,  &
  &             flag_psi  = 0, flag_alph  = 0,  flag_vep  = 0,  &
  &             flag_bvxd = 0, flag_bvyd  = 0,  flag_bvzd = 0
  integer    :: irg, itg, ipg, ii, iseq, jj, sumjj, je
  integer    :: impt, impt_ex, impt_bin
  character(len=2) :: chgreen, chpa, chpB, cocp
  character(len=4) :: chbe
  integer    :: flag = 0, hydro_iter = 4, ishift, shift_iter=1, &
  &             flag_mid,   flag_mpt_mid(10)
  integer    :: irf, itf, ipf, ihy, npg_l, npg_r, icount
  character(30) :: char1, char2, char3, char4, char5
  real(long) :: start, finish
  character(len=1) :: char_1
!
  call set_allocate_size_mpt
  call alloc_array4d(souvec ,0,nrg,0,ntg,0,npg,1,3)
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
  call alloc_array2d(sou_insurf,0,ntg,0,npg)
  call alloc_array2d(dsou_insurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
  call alloc_array3d( potf,0,nrf,0,ntf,0,npf)
  call alloc_array3d(potut,0,nrf,0,ntf,0,npf)
  call alloc_array3d(potux,0,nrf,0,ntf,0,npf)
  call alloc_array3d(potuy,0,nrf,0,ntf,0,npf)
  call alloc_array3d(potuz,0,nrf,0,ntf,0,npf)
  call alloc_array3d(emd_bak,0,nrf,0,ntf,0,npf)
  call alloc_array3d(vep_bak,0,nrf,0,ntf,0,npf)
  call alloc_array2d(potrs,0,ntf,0,npf)
!
  iter_count = 0
  error_emd_last =1.0d0
  error_psi_last =1.0d0
  error_alph_last=1.0d0
  flag_mid=0
!
  do
    iter_count = iter_count + 1
!    iter_extra = iter_extra + 1
    count = dble(iter_count)
    flag_mpt_mid(1:10) =0    !  one flag for every patch
    write(6,'(a10,i4,a108)') "Iteration:",iter_count, &
    &  " ---psi------------alph------------bvxd------------bvyd------------bvzd------------emd------------vep-------"
!
    do impt = 1, nmpt
      write(6,'(a8,i2,a1)', ADVANCE = "NO")  " Patch #", impt, ""
      write(6,*)  " "
      call copy_grid_parameter_from_mpt(impt)
      conv_gra = dmin1(conv0_gra,conv_ini*count)
      conv_den = dmin1(conv0_den,conv_ini*count)
      conv_vep = dmin1(conv0_vep,conv_ini*count)
      call copy_grid_parameter_to_mpt(impt)
! ---
      call copy_from_mpatch_all_BNS_CF(impt)
      call copy_def_metric_and_matter_from_mpt(impt)
      call copy_def_matter_irrot_from_mpt(impt)
! ---
      if (iter_count==1)   write(6,*)  "BBBBBB emdc=", emdc

!      write(6,'(a20,i5,i5,1p,3e20.12)') "================", iter_count, impt, ome, ber, radi

!!!!      call excurve
!      emdg = 0.0d0
      if (impt==1 .or. impt==2) then
        cocp = 'ns'
        call interpo_fl2gr(emd, emdg)
        call interpo_fl2gr(vepxf, vepxg)
        call interpo_fl2gr(vepyf, vepyg)
        call interpo_fl2gr(vepzf, vepzg)
!       ********************************   CHECK   ***************************************************
        if (mod(iter_count,5)==0.and.impt==10) then
          write(char1, '(i5)') iter_count
          char2 = adjustl(char1)
          write(char4, '(i5)') impt
          char5 = adjustl(char4)

          char3 = 'vep_iteration' // trim(char2) // '_mpt' // trim(char5) // '.txt'
          open(12,file=char3,status='unknown')
          itg=12;  ipg=6
          write(12,'(a1,a2,4a16,a30)') '#', 'ir', 'rg*rs', 'vepxf', 'vepyf', 'vepzf', 'plot using every :::0::0'
          do irg=0,nrf
            write(12,'(i3,1p,4e16.6)')  irg, rg(irg)*rs(itg,ipg), vepxf(irg,itg,ipg), vepyf(irg,itg,ipg), vepzf(irg,itg,ipg)
          end do
          write(12,*) '#------------------------------------------------------------------------------------------------------'
          write(12,*) ""
          write(12,'(a1,a2,4a16,a30)') '#', 'ir', 'rg', 'vepxg', 'vepyg', 'vepzg',  'plot using every :::1::1'
          do irg=0,nrg
            write(12,'(i3,1p,4e16.6)')  irg, rg(irg), vepxg(irg,itg,ipg), vepyg(irg,itg,ipg), vepzg(irg,itg,ipg)
          end do
          write(12,*) ""
          close(12)
        endif
!       *********************************************************************************************

        call calc_vector_x_matter(2)
        call calc_vector_phi_matter(2)
        call excurve_CF(cocp)             !   3rd order midpoint from ir0=-2,...
        call excurve_CF_gridpoint         !   4th order from ir0=-2,...
        call copy_def_metric_and_matter_to_mpt(impt)
        call copy_def_matter_irrot_to_mpt(impt)
        if (iter_count==1 .and. impt==1) then
          open(2,file='iter_data.txt',status='unknown')
        end if
        write(2,'(i5,i5,1p,3e20.12)') iter_count, impt, ome, ber, radi 
      else
        cocp = 'bh'
        call excurve_CF(cocp)             !   3rd order midpoint from ir0=0,...
        call excurve_CF_gridpoint_bhex    !   4th order from ir0=0,...
        call copy_def_metric_and_matter_to_mpt(impt)
        call copy_def_matter_irrot_to_mpt(impt)
        write(2,'(i5,i5,1p,3e20.12)') iter_count, impt, ome, ber, radi 
        write(2,*)  '-----------------------------------------------------------------------'
      end if
!
      if (impt.eq.nmpt) call calc_mass_asympto('ns')
      call copy_to_mpatch_all_BNS_CF(impt)
!
! --- psi
      if(error_psi_last.le.1.0d-13) then
        write(6,'(a16)', ADVANCE = "NO")    ' ## SKIP psi  ##'
        write(6,*) " "
      else
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
      if (impt==1 .or. impt==2) then
        call sourceterm_HaC_CF_peos_irrot(sou)
!        char1 = 'psi_sou.txt'
!        call IO_output_1D_general(char1, 'g', 'm', sou, -1, ntg/2, 1)
      else
        call sourceterm_HaC_CF(sou)
      end if

      if (impt==1 .or. impt==2) then
        call sourceterm_exsurf_binary_COCP(impt,'psi ', 'ev', sou_exsurf, dsou_exsurf)
        call sourceterm_outsurf_COCP_from_ARCP(impt,'psi ', 'ev', sou_outsurf, dsou_outsurf)

        call cpu_time(start)
        call poisson_solver_binary_star_homosol('sd', sou, sou_exsurf, dsou_exsurf, &
        &                                                  sou_outsurf,dsou_outsurf,pot)

!        call poisson_solver_binary_star_homosol_plmex(impt,'sd', sou, sou_exsurf, dsou_exsurf, &
!        &                                                             sou_outsurf,dsou_outsurf,pot)
        call cpu_time(finish)
        print '("Time = ",f6.3," seconds.")',finish-start

      else if (impt.eq.nmpt) then
        call sourceterm_insurf_ARCP_from_COCP(impt,'psi ', 'ev', sou_insurf, dsou_insurf)
        call outer_boundary_d_psi(sou_outsurf)
!        call outer_boundary_d_BBH_CF(sou_outsurf,'psi ')
        call poisson_solver_asymptotic_patch_homosol('dd',sou, sou_insurf, dsou_insurf, &
        &                                                      sou_outsurf,dsou_outsurf,pot)
      end if

      pot_bak(0:nrg,0:ntg,0:npg) = psi(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,psi,conv_gra)

      if ((impt==1 .or. impt==2) .and. (iter_count>5 .or. iseq>1 .or. indata_type.eq.'3D')) &
        &          call update_parameter_BNS_irrot(conv_gra)

      if (impt.eq.nmpt)       psi(nrg,0:ntg,0:npg) = 1.0d0 
!      call reset_outer_boundary_BBH_CF('psi ')
      if (impt.ne.nmpt) call interpolation_fillup_binary_COCP(impt,'psi ','ev',psi)
!!!      call error_metric_type0(psi,pot_bak,error_tmp,flag_tmp,cocp)
      call error_metric_type2_mpt(psi,pot_bak,error_tmp,flag_tmp,cocp,impt)
      flag_psi  =  max ( flag_psi, flag_tmp)
      error_psi = dmax1(error_psi,error_tmp)
!!!      write(6,'(1p,e14.6,a2)', ADVANCE = "NO")  error_psi, "  "
      error_psi_last = error_psi
!      call printout_error_metric_v1('psi ', iter_count,error_psi)
      call copy_to_mpatch_all_BNS_CF(impt)
      call copy_def_metric_and_matter_to_mpt(impt)
      call copy_def_matter_irrot_to_mpt(impt)
      end if
!
! --- alps
      if(error_alph_last.le.1.0d-13) then
        write(6,'(a16)', ADVANCE = "NO")    ' ## SKIP alph ##'
        write(6,*) " "
      else
      sou(0:nrg,0:ntg,0:npg)  = 0.0d0
      sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
      sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
      alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)
      if (impt==1 .or. impt==2) then
        call sourceterm_trG_CF_peos_irrot(sou)
!        char1 = 'alps_sou.txt'
!        call IO_output_1D_general(char1, 'g', 'm', sou, -1, ntg/2, 1)
      else
        call sourceterm_trG_CF_helical(sou)
      end if
      if (impt==1 .or. impt==2) then
        call sourceterm_exsurf_binary_COCP(impt,'alps','ev',sou_exsurf, dsou_exsurf)
        call sourceterm_outsurf_COCP_from_ARCP(impt,'alps','ev',sou_outsurf, dsou_outsurf)

        call poisson_solver_binary_star_homosol('sd',sou, sou_exsurf, dsou_exsurf, &
        &                                                 sou_outsurf,dsou_outsurf,pot)

!        call poisson_solver_binary_star_homosol_plmex(impt,'sd',sou, sou_exsurf, dsou_exsurf, &
!        &                                                            sou_outsurf,dsou_outsurf,pot)
      else if (impt.eq.nmpt) then
        call sourceterm_insurf_ARCP_from_COCP(impt,'alps','ev',sou_insurf, dsou_insurf)
        call outer_boundary_d_alps(sou_outsurf)
!        call outer_boundary_d_BBH_CF(sou_outsurf,'alps')
        call poisson_solver_asymptotic_patch_homosol('dd',sou, sou_insurf,dsou_insurf, & 
        &                                                      sou_outsurf,dsou_outsurf,pot)
      end if

      call compute_alps2alph(pot,psi)
      pot_bak(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)
      call update_grfield(pot,alph,conv_gra)

      if ((impt==1 .or. impt==2) .and. (iter_count>5 .or. iseq>1 .or. indata_type.eq.'3D')) &
        &          call update_parameter_BNS_irrot(conv_gra)

      if (impt.eq.nmpt)       alph(nrg,0:ntg,0:npg) = 1.0d0
!      call reset_outer_boundary_BBH_CF('alph')
      if (impt.ne.nmpt) call interpolation_fillup_binary_COCP(impt,'alph','ev',alph)
!!!      call error_metric_type0(alph,pot_bak,error_tmp,flag_tmp,cocp)
      call error_metric_type2_mpt(alph,pot_bak,error_tmp,flag_tmp,cocp,impt)
      flag_alph  =  max ( flag_alph, flag_tmp)
      error_alph = dmax1(error_alph,error_tmp)
!!!      write(6,'(1p,e14.6,a2)', ADVANCE = "NO")  error_alph, "  "
      error_alph_last = error_alph
!      call printout_error_metric_v1('alph', iter_count,error_alph)
      alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)
      call copy_to_mpatch_all_BNS_CF(impt)
      call copy_def_metric_and_matter_to_mpt(impt)
      call copy_def_matter_irrot_to_mpt(impt)
      end if
! 
! --- shift
      do ishift=1,shift_iter
      souvec(0:nrg,0:ntg,0:npg,1:3)  = 0.0d0
      if (impt==1 .or. impt==2) then
        call sourceterm_MoC_CF_with_divshift_peos_irrot(souvec)
      else
        call sourceterm_MoC_CF_with_divshift_r3rd(souvec)
      end if

      do ii = 1, 3
        sou_exsurf(0:ntg,0:npg) = 0.0d0 ; dsou_exsurf(0:ntg,0:npg) = 0.0d0
        sou_outsurf(0:ntg,0:npg)= 0.0d0 ; dsou_outsurf(0:ntg,0:npg)= 0.0d0
        sou(0:nrg,0:ntg,0:npg) = souvec(0:nrg,0:ntg,0:npg,ii)

!        write(char1, '(i5)') ii
!        char2 = adjustl(char1)
!        write(char4, '(i5)') impt
!        char5 = adjustl(char4)
!        char3 = 'beta_sou' // trim(char2) // '_mpt' // trim(char5) // '.txt'
!        call IO_output_1D_general(char3,'g','m', sou, -1, ntg/2, 1)

        if (ii.eq.1) work(0:nrg,0:ntg,0:npg) = bvxd(0:nrg,0:ntg,0:npg)
        if (ii.eq.2) work(0:nrg,0:ntg,0:npg) = bvyd(0:nrg,0:ntg,0:npg)
        if (ii.eq.3) work(0:nrg,0:ntg,0:npg) = bvzd(0:nrg,0:ntg,0:npg)
        if (ii.eq.1) then; chbe = 'bvxd'; chpa = 'od'; end if
        if (ii.eq.2) then; chbe = 'bvyd'; chpa = 'od'; end if
        if (ii.eq.3) then; chbe = 'bvzd'; chpa = 'ev'; end if
        if (impt==1 .or. impt==2) then
          chpB = 'ev'
          if (impt.eq.2.and.ii.eq.1) chpB = 'od'
          if (impt.eq.2.and.ii.eq.2) chpB = 'od'
          call sourceterm_exsurf_binary_COCP(impt,chbe,chpa,sou_exsurf, dsou_exsurf)
          call sourceterm_outsurf_COCP_from_ARCP(impt,chbe,chpB,sou_outsurf, dsou_outsurf)

          call poisson_solver_binary_star_homosol('sd',sou, sou_exsurf,dsou_exsurf, &
          &                                                 sou_outsurf,dsou_outsurf,pot)

!          call poisson_solver_binary_star_homosol_plmex(impt,'sd',sou,sou_exsurf,dsou_exsurf, &
!          &                                                           sou_outsurf,dsou_outsurf,pot)
        else if (impt.eq.nmpt) then
          call sourceterm_insurf_ARCP_from_COCP(impt,chbe,chpa,sou_insurf, dsou_insurf)
          if (ii.eq.1) call outer_boundary_d_bvxd(sou_outsurf)
          if (ii.eq.2) call outer_boundary_d_bvyd(sou_outsurf)
          if (ii.eq.3) call outer_boundary_d_bvzd(sou_outsurf)
!          if (ii.eq.1) call outer_boundary_d_BBH_CF(sou_outsurf,'bvxd')
!          if (ii.eq.2) call outer_boundary_d_BBH_CF(sou_outsurf,'bvyd')
!          if (ii.eq.3) call outer_boundary_d_BBH_CF(sou_outsurf,'bvzd')
          call poisson_solver_asymptotic_patch_homosol('dd',sou, sou_insurf,dsou_insurf, &
          &                                                      sou_outsurf,dsou_outsurf,pot)
        end if
        if (ii.eq.1) potx(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
        if (ii.eq.2) poty(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
        if (ii.eq.3) potz(0:nrg,0:ntg,0:npg) = pot(0:nrg,0:ntg,0:npg)
      end do

      potx_bak(0:nrg,0:ntg,0:npg) = bvxd(0:nrg,0:ntg,0:npg)
      poty_bak(0:nrg,0:ntg,0:npg) = bvyd(0:nrg,0:ntg,0:npg)
      potz_bak(0:nrg,0:ntg,0:npg) = bvzd(0:nrg,0:ntg,0:npg)
      call update_grfield(potx,bvxd,conv_gra)
      call update_grfield(poty,bvyd,conv_gra)
      call update_grfield(potz,bvzd,conv_gra)

      if ((impt==1 .or. impt==2) .and. (iter_count>5 .or. iseq>1 .or. indata_type.eq.'3D')) &
        &          call update_parameter_BNS_irrot(conv_gra)

      if (impt.ne.nmpt) then 
        call interpolation_fillup_binary_COCP(impt,'bvxd','od',bvxd)
        call interpolation_fillup_binary_COCP(impt,'bvyd','od',bvyd)
        call interpolation_fillup_binary_COCP(impt,'bvzd','ev',bvzd)
      end if
!      call reset_metric_CF
      call error_metric_type2_mpt(bvxd,potx_bak,error_tmp,flag_tmp,cocp,impt)
      flag_bvxd  =  max ( flag_bvxd, flag_tmp)
      error_bvxd = dmax1(error_bvxd,error_tmp)
      call error_metric_type2_mpt(bvyd,poty_bak,error_tmp,flag_tmp,cocp,impt)
      flag_bvyd  =  max ( flag_bvyd, flag_tmp)
      error_bvyd = dmax1(error_bvyd,error_tmp)
      call error_metric_type2_mpt(bvzd,potz_bak,error_tmp,flag_tmp,cocp,impt)
      flag_bvzd  =  max ( flag_bvzd, flag_tmp)
      error_bvzd = dmax1(error_bvzd,error_tmp)
!!!      if(ishift==shift_iter)  write(6,'(1p,e14.6,a2)', ADVANCE = "NO")  error_bvxd, "  "
!!!      if(ishift==shift_iter)  write(6,'(1p,e14.6,a2)', ADVANCE = "NO")  error_bvyd, "  "
!!!      if(ishift==shift_iter)  write(6,'(1p,e14.6,a2)', advance = "no")  error_bvzd, "  "
      end do   ! ==========>> ishift 

!      call printout_error_metric_v1('bvxd', iter_count,error_bvxd)
!      call printout_error_metric_v1('bvyd', iter_count,error_bvyd)
!      call printout_error_metric_v1('bvzd', iter_count,error_bvzd)
      call copy_to_mpatch_all_BNS_CF(impt)
      call copy_def_metric_and_matter_to_mpt(impt)
      call copy_def_matter_irrot_to_mpt(impt)
!
      error_metr = dmax1(error_psi,error_alph,error_bvxd,error_bvyd,error_bvzd)

      if(iter_count<5 .and. iseq<=1 .and. indata_type.eq.'1D') then
        write(6,'(a14)') '#SKIP HYDRO a#' 
        error_emd_last = 1.0d0 
        cycle
      end if
! ---
! -- Hydro equations.
!      if (error_emd_last.le.eps*0.0001d0 .and. error_metr.gt.eps*2.0d0) then
      if (error_emd_last.le.eps*0.0001d0 .and. error_vep_last.le.eps*0.0001d0 ) then
        write(6,'(a14)')  '#SKIP HYDRO b#'
      else
        if (impt==1 .or. impt==2) then
          do ihy = 1, hydro_iter
            call interpo_gr2fl_metric_CF
            call hydrostatic_eq_CF_peos_irrot(potf,potut,potux,potuy,potuz)
            call calc_surface(potrs, potf)
            emd_bak = emd 
            call update_matter( potf,emd,conv_den)
!            call update_matter(potut,utf,conv_den)
!            call update_matter(potux,vepxf,conv_den)
!            call update_matter(potuy,vepyf,conv_den)
!            call update_matter(potuz,vepzf,conv_den)
            call update_surface(potrs,rs,conv_den)

            call reset_fluid_BNS
            call update_parameter_BNS_irrot(conv_den)
            call calc_vector_x_matter(2)
            call calc_vector_phi_matter(2)
!
            call hydro_irbns_vep_CF_peos(potf,iter_count,impt,10)
            vep_bak = vep
            call update_matter(potf,vep,conv_vep)
            call reset_fluid_vep
            call calc_gradvep(vep,vepxf,vepyf,vepzf)
            call reset_fluid_gradvep
            call update_parameter_BNS_irrot(conv_vep)

            if (ihy.eq.hydro_iter) then
!            if (ihy.eq.1) then
              call error_matter_type2(emd,emd_bak,error_tmp,flag_tmp)
              flag_emd  =  max ( flag_emd, flag_tmp)
              error_emd = dmax1(error_emd,error_tmp)
              error_emd_last = error_emd
!!!              write(6,'(1p,e14.6,a2)', ADVANCE = "NO")  error_emd, "  "
!              call printout_error_matter("emd ", iter_count,error_emd)
              call error_velocity_potential(vep,vep_bak,error_tmp,flag_tmp)
              flag_vep  =  max ( flag_vep, flag_tmp)
              error_vep = dmax1(error_vep,error_tmp)
              error_vep_last = error_vep
!!!              write(6,'(1p,e14.6,a2)')  error_vep, "  "        
            end if
          end do   ! end hydro
          call copy_to_mpatch_all_BNS_CF(impt)
          call copy_def_metric_and_matter_to_mpt(impt)
          call copy_def_matter_irrot_to_mpt(impt)
        else
          write(6,*)   ' ##   NO HYDRO ##'
        end if
      end if
!      error_emd_last = error_emd
!
      if(impt>5) then
        if (mod(iter_count,1)==0.or.iter_count==1) then
          write(char1, '(i5)') iter_count
          char2 = adjustl(char1)
          write(char4, '(i5)') impt
          char5 = adjustl(char4)
          char3 = 'iteration' // trim(char2) // '_mpt' // trim(char5) // '.txt'
          call printout_axis_irrot_BNS_mpt(impt,char3)
        endif
      end if  !  ====> iteration.txt

!     print once the solution when the error in shift for all patches is <0.1 
      if (flag_mid==0) then
        if(error_metr < iter_eps ) then
          flag_mpt_mid(impt) = 1
          error_metric_mpt_mid(impt,1) = error_psi
          error_metric_mpt_mid(impt,2) = error_alph
          error_metric_mpt_mid(impt,3) = error_bvxd
          error_metric_mpt_mid(impt,4) = error_bvyd
          error_metric_mpt_mid(impt,5) = error_bvzd
        end if
        sumjj=0
        do jj=1,nmpt 
          sumjj = sumjj + flag_mpt_mid(jj)
        end do
        if (sumjj==nmpt) then
          open(3,file='mid_error_data.txt',status='unknown')
          write(3,'(a11,i5)') 'Iteration =', iter_count
          write(3,'(a8,1p,5e17.6)') 'COCP #=1', (error_metric_mpt_mid(1,je), je=1,5)
          write(3,'(a8,1p,5e17.6)') 'COCP #=2', (error_metric_mpt_mid(2,je), je=1,5)
          write(3,'(a8,1p,5e17.6)') 'ARCP #  ', (error_metric_mpt_mid(3,je), je=1,5)
          close(3)
!          do jj=1,nmpt
!            call copy_from_mpatch_all_BNS_CF(jj)
!            call copy_def_metric_and_matter_from_mpt(jj)
!            call copy_def_matter_irrot_from_mpt(jj)
!            call IO_output_solution_3D_CF_irrot_NS_mpt(jj,'.mid')
!            if (jj==1 .or. jj==2)  call printout_NS_shape_mpt(jj)
!          end do
          flag_mid = 1
        end if 
      end if
    end do    !  ====================>> end all patches
!
! -- convergence of iteration is checked.
!
!    if (iter_count==6) then
!      do impt=1,nmpt
!        call copy_from_mpatch_all_BNS_CF(impt)
!        write(2,'(i5,i5,1p,3e20.12)') iter_count, impt, ome, ber, radi
!      end do
!      write(2,*) '-----------------------------------------------------------------------'
!      close(2)
!      exit
!    end if

    if (flag_mid==1) then
      do impt=1,nmpt
        call copy_from_mpatch_all_BNS_CF(impt)
!        call copy_def_metric_and_matter_from_mpt(impt)
!        call copy_def_matter_irrot_from_mpt(impt)
        write(2,'(i5,i5,1p,3e20.12)') iter_count, impt, ome, ber, radi

!        MoC_vio(0:nrg,0:ntg,0:npg,1:3) = 0.0d0
!        call violation_gridpoint_MoC_CF_peos_irrot(MoC_vio,'ns')
!        pot = 0.0d0
!        pot(0:nrg,0:ntg,0:npg) = MoC_vio(0:nrg,0:ntg,0:npg,2)   
!        char1 = 'MoC_vio_by_3D_mpt1.txt'
!        call IO_output_3D_general(char1,'g',pot)
!        char1 = 'MoC_vio_by_1D_mpt1.txt'
!        call IO_output_1D_general(char1,'g',pot,-1,ntg/2,0)

!        HaC_vio = 0.0d0
!        call violation_gridpoint_HaC_CF_peos_irrot(HaC_vio,'ns')
!        char1 = 'HaC_vio_3D_mpt1.txt'
!        call IO_output_3D_general(char1,'g',HaC_vio)
!        char1 = 'HaC_vio_1D_mpt1.txt'
!        call IO_output_1D_general(char1,'g',HaC_vio,-1,ntg/2,0)
      end do
      write(2,*) '-----------------------------------------------------------------------'
      close(2)
      exit
    end if

    if ( (flag_vep==0).and.(flag_emd==0).and.(flag_psi==0).and.(flag_alph==0).and. &
    &   (flag_bvxd==0).and.(flag_bvyd==0).and.(flag_bvzd==0))      then
      do impt=1, nmpt
        call copy_from_mpatch_all_BNS_CF(impt)
        write(2,'(i5,i5,1p,3e20.12)') iter_count, impt, ome, ber, radi
      end do
      write(2,*) '-----------------------------------------------------------------------'
      close(2)
      exit
    endif

    if (iter_count >= iter_max) then
      write(6,*) "*************** iter_count==iter_max ....exiting *******************"
      stop
!      exit
    end if
!    error_emd_last = error_emd
!
    flag_psi  = 0
    flag_alph = 0
    flag_bvxd = 0
    flag_bvyd = 0
    flag_bvzd = 0
    flag_emd  = 0
    flag_vep  = 0
    error_psi  = 0.0d0
    error_alph = 0.0d0
    error_bvxd = 0.0d0
    error_bvyd = 0.0d0
    error_bvzd = 0.0d0
    error_emd  = 0.0d0
    error_vep  = 0.0d0
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
  deallocate(sou_insurf)
  deallocate(dsou_insurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
  deallocate(potf)
  deallocate(potut)
  deallocate(potux)
  deallocate(potuy)
  deallocate(potuz)
  deallocate(emd_bak)
  deallocate(vep_bak)
  deallocate(potrs)
  
end subroutine iter_irrot_BNS_CF_mpt
