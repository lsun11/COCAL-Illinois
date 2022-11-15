!______________________________________________
include '../Include_file/include_modulefiles_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_modulefiles_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_interface_modulefiles_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_subroutines_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_subroutines_analysis_plot_BNS_CF_3mpt.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!      COCAL ID to CACTUS
!______________________________________________
PROGRAM coc2cac
!
  use grid_parameter_binary_excision
  use phys_constant
  use grid_parameter
  use interface_modules_cartesian
  use interface_calc_gradvep_export
  use trigonometry_grav_phi
  use def_binary_parameter, only : dis
  use def_matter_parameter_mpt
  use interface_IO_input_CF_grav_export
  use interface_IO_input_CF_flsp_export
  use interface_IO_input_CF_flco_export
  use interface_IO_input_CF_flir_export
  use interface_IO_input_CF_surf_export
  use interface_excurve_CF_gridpoint_export
  use interface_interpo_gr2fl_metric_CF_export
  implicit none
  integer :: impt, impt_ex, ico, irr, isp
  real(8) :: confpow, psifcacp, wxspfca, wyspfca, wzspfca
  character(30) :: char1
  character*400 :: dir_path
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
  real(8) :: rr3, rrcm, xc,yc,zc, xc_p1, yc_p1, zc_p1, xc_p2, yc_p2, zc_p2, dis_cm
  real(8) :: xc_p3, yc_p3, zc_p3
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc
  real(8) :: emdca, vepca, psica, alphca, bvxdca, bvydca, bvzdca, psi4ca, psif4ca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: vepxfca, vepyfca, vepzfca, vxu, vyu, vzu, lam_p1, lam_p2
  real(8) :: bxcor, bycor, bzcor, bvxdfca, bvydfca, bvzdfca, psifca, alphfca
  real(8) :: gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz
  real(8) :: axx, axy, axz, ayy, ayz, azz
  real(8) :: ome_p1, ber_p1, radi_p1, r_surf_p1
  real(8) :: ome_p2, ber_p2, radi_p2, r_surf_p2
!
  real(8), pointer :: emd_p1(:,:,:), vep_p1(:,:,:), rs_p1(:,:)
  real(8), pointer :: vepxf_p1(:,:,:), vepyf_p1(:,:,:), vepzf_p1(:,:,:)
  real(8), pointer :: wxspf_p1(:,:,:), wyspf_p1(:,:,:), wzspf_p1(:,:,:)
  real(8), pointer :: psif_p1(:,:,:), alphf_p1(:,:,:), bvxdf_p1(:,:,:), bvydf_p1(:,:,:), bvzdf_p1(:,:,:)
  real(8), pointer :: psi_p1(:,:,:), alph_p1(:,:,:), bvxd_p1(:,:,:), bvyd_p1(:,:,:), bvzd_p1(:,:,:)
  real(8), pointer :: axx_p1(:,:,:), axy_p1(:,:,:) , axz_p1(:,:,:) , ayy_p1(:,:,:) , ayz_p1(:,:,:), azz_p1(:,:,:)
!
  real(8), pointer :: emd_p2(:,:,:), vep_p2(:,:,:), rs_p2(:,:)
  real(8), pointer :: vepxf_p2(:,:,:), vepyf_p2(:,:,:), vepzf_p2(:,:,:)
  real(8), pointer :: wxspf_p2(:,:,:), wyspf_p2(:,:,:), wzspf_p2(:,:,:)
  real(8), pointer :: psif_p2(:,:,:), alphf_p2(:,:,:), bvxdf_p2(:,:,:), bvydf_p2(:,:,:), bvzdf_p2(:,:,:)
  real(8), pointer :: psi_p2(:,:,:), alph_p2(:,:,:), bvxd_p2(:,:,:), bvyd_p2(:,:,:), bvzd_p2(:,:,:)
  real(8), pointer :: axx_p2(:,:,:), axy_p2(:,:,:) , axz_p2(:,:,:) , ayy_p2(:,:,:) , ayz_p2(:,:,:), azz_p2(:,:,:)
!
  real(8), pointer :: psi_p3(:,:,:), alph_p3(:,:,:), bvxd_p3(:,:,:), bvyd_p3(:,:,:), bvzd_p3(:,:,:)
  real(8), pointer :: axx_p3(:,:,:), axy_p3(:,:,:) , axz_p3(:,:,:) , ayy_p3(:,:,:) , ayz_p3(:,:,:), azz_p3(:,:,:)
!
  integer :: kij_parity
  confpow = -6.0d0 
!
  gxx=0.0d0; gxy=0.0d0; gxz=0.0d0; gyy=0.0d0; gyz=0.0d0; gzz=0.0d0
  kxx=0.0d0; kxy=0.0d0; kxz=0.0d0; kyy=0.0d0; kyz=0.0d0; kzz=0.0d0
  axx=0.0d0; axy=0.0d0; axz=0.0d0; ayy=0.0d0; ayz=0.0d0; azz=0.0d0

  !TODO remove this
  !dir_path="/home/astro/mundim/tmp/ET_2014_05_wheeler/Cactus/repos/Cocal/standalone/Cocal/ID_BNS"
  !dir_path="../../standalone/Cocal/ID_BNS"
  dir_path='.'

! -- Read parameters
  call allocate_grid_parameter_mpt
  call allocate_grid_parameter_binary_excision_mpt
  call allocate_def_matter_parameter_mpt
  do impt = 1, nmpt
    call read_parameter_mpt_cactus(impt,dir_path)
    if (impt .le. 2) call read_surf_parameter_mpt_cactus(impt,dir_path)
    call copy_grid_parameter_to_mpt(impt)
    call read_parameter_binary_excision_mpt_cactus(impt,dir_path)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
    if (impt .le. 2) call peos_initialize_mpt_cactus(impt,dir_path)
    call copy_def_peos_parameter_to_mpt(impt)
  end do
!
  call set_allocate_size_mpt
  call allocate_trig_grav_mphi
  call allocate_trigonometry_grav_phi_mpt
!
  do impt = 1, nmpt
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_def_peos_parameter_from_mpt(impt)
    call coordinate_patch_kit_grav_grid_coc2cac_mpt(3)  ! 3:r_surf is used
    call calc_parameter_binary_excision
    call copy_grid_parameter_to_mpt(impt)
    call copy_grid_parameter_binary_excision_to_mpt(impt)
    call copy_coordinate_grav_extended_to_mpt(impt)
    call copy_coordinate_grav_phi_to_mpt(impt)
    call copy_coordinate_grav_r_to_mpt(impt)
    call copy_coordinate_grav_theta_to_mpt(impt)
    call copy_def_binary_parameter_to_mpt(impt)
    call copy_trigonometry_grav_phi_to_mpt(impt)
    call copy_trigonometry_grav_theta_to_mpt(impt)
  end do

  do impt = 1,nmpt
!=>    call copy_from_mpatch_interpolation_utility(impt)
    call copy_grid_parameter_from_mpt(impt)
    call copy_grid_parameter_binary_excision_from_mpt(impt)
    call copy_coordinate_grav_extended_from_mpt(impt)
    call copy_coordinate_grav_phi_from_mpt(impt)
    call copy_coordinate_grav_r_from_mpt(impt)
    call copy_coordinate_grav_theta_from_mpt(impt)
    call copy_def_binary_parameter_from_mpt(impt)
    call copy_trigonometry_grav_phi_from_mpt(impt)
    call copy_trigonometry_grav_theta_from_mpt(impt)

!    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
    if (impt==1) then
      if(NS_shape.eq.'SP')  then
        ico=0; irr=0; isp=1  
        write(6,*) "************ Spinning configuration ************"
      end if
      if(NS_shape.eq.'IR')  then
        ico=0; irr=1; isp=0
        write(6,*) "************ Irrotational configuration ************"
      end if
      if(NS_shape.eq.'CO')  then
        ico=1; irr=0; isp=0
        write(6,*) "************ Corotating configuration ************"
      end if

      rr3 = 0.7d0*(rgout - rgmid)
      dis_cm = dis
      r_surf_p1 = r_surf
      write(6,'(a13,i1,a8)') "Reading COCP-", impt, " data..."

      allocate (  emd_p1(0:nrf,0:ntf,0:npf))
      allocate (  vep_p1(0:nrf,0:ntf,0:npf))
      allocate (wxspf_p1(0:nrf,0:ntf,0:npf))
      allocate (wyspf_p1(0:nrf,0:ntf,0:npf))
      allocate (wzspf_p1(0:nrf,0:ntf,0:npf))
      allocate (vepxf_p1(0:nrf,0:ntf,0:npf))
      allocate (vepyf_p1(0:nrf,0:ntf,0:npf))
      allocate (vepzf_p1(0:nrf,0:ntf,0:npf))
      allocate ( psif_p1(0:nrf,0:ntf,0:npf))
      allocate (alphf_p1(0:nrf,0:ntf,0:npf))
      allocate (bvxdf_p1(0:nrf,0:ntf,0:npf))
      allocate (bvydf_p1(0:nrf,0:ntf,0:npf))
      allocate (bvzdf_p1(0:nrf,0:ntf,0:npf))
      allocate (   rs_p1(0:ntf,0:npf))
      allocate (  psi_p1(0:nrg,0:ntg,0:npg))
      allocate ( alph_p1(0:nrg,0:ntg,0:npg))
      allocate ( bvxd_p1(0:nrg,0:ntg,0:npg))
      allocate ( bvyd_p1(0:nrg,0:ntg,0:npg))
      allocate ( bvzd_p1(0:nrg,0:ntg,0:npg))
      allocate (  axx_p1(0:nrg,0:ntg,0:npg))
      allocate (  axy_p1(0:nrg,0:ntg,0:npg))
      allocate (  axz_p1(0:nrg,0:ntg,0:npg))
      allocate (  ayy_p1(0:nrg,0:ntg,0:npg))
      allocate (  ayz_p1(0:nrg,0:ntg,0:npg))
      allocate (  azz_p1(0:nrg,0:ntg,0:npg))
      emd_p1=0.0d0;  vep_p1 =0.0d0;  rs_p1  =0.0d0; wxspf_p1=0.0d0; wyspf_p1=0.0d0; wzspf_p1=0.0d0
      psi_p1=0.0d0;  alph_p1=0.0d0;  bvxd_p1=0.0d0;  bvyd_p1=0.0d0;  bvzd_p1=0.0d0
      axx_p1=0.0d0;  axy_p1 =0.0d0;  axz_p1 =0.0d0;   ayy_p1=0.0d0;   ayz_p1=0.0d0;   azz_p1=0.0d0

      !TODO: Change the following:
!      call IO_input_CF_irrot_BNS_cactus(impt,emd_p1,vep_p1,ome_p1,ber_p1,radi_p1,rs_p1, &
!         &    psi_p1,alph_p1,bvxd_p1,bvyd_p1,bvzd_p1,dir_path)

      call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt1.las",psi_p1,alph_p1,bvxd_p1,bvyd_p1,bvzd_p1)

      if(NS_shape.eq.'CO')   call IO_input_CF_flco_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",emd_p1,ome_p1,ber_p1,radi_p1)
      if(NS_shape.eq.'IR')   call IO_input_CF_flir_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",emd_p1,vep_p1,ome_p1,ber_p1,radi_p1)
      if(NS_shape.eq.'SP')   call IO_input_CF_flsp_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",emd_p1,vep_p1,wxspf_p1,wyspf_p1, &
                                  &   wzspf_p1,ome_p1,ber_p1,radi_p1)

      call IO_input_CF_surf_export(trim(dir_path)//"/bnssur_3D_mpt1.las",rs_p1)

      call excurve_CF_gridpoint_export(alph_p1,bvxd_p1,bvyd_p1,bvzd_p1, & 
         &    axx_p1, axy_p1, axz_p1, ayy_p1, ayz_p1, azz_p1)

      call interpo_gr2fl_metric_CF_export(alph_p1, psi_p1, bvxd_p1, bvyd_p1, bvzd_p1, &
         &    alphf_p1, psif_p1, bvxdf_p1, bvydf_p1, bvzdf_p1, rs_p1)

      if((NS_shape.eq.'IR').or.(NS_shape.eq.'SP'))  call calc_gradvep_export(vep_p1,vepxf_p1,vepyf_p1,vepzf_p1,rs_p1)
    end if
    if (impt==2) then
      if(NS_shape.eq.'SP')  then
        write(6,*) "************ Spinning configuration ************"
        if (isp.ne.1)  stop
      end if
      if(NS_shape.eq.'IR')  then
        write(6,*) "************ Irrotational configuration ************"
        if(irr.ne.1)  stop
      endif
      if(NS_shape.eq.'CO')  then
        write(6,*) "************ Corotating configuration ************"
        if(ico.ne.1)  stop
      end if
      r_surf_p2 = r_surf
      write(6,'(a13,i1,a8)') "Reading COCP-", impt, " data..."
      allocate (  emd_p2(0:nrf,0:ntf,0:npf))
      allocate (  vep_p2(0:nrf,0:ntf,0:npf))
      allocate (wxspf_p2(0:nrf,0:ntf,0:npf))
      allocate (wyspf_p2(0:nrf,0:ntf,0:npf))
      allocate (wzspf_p2(0:nrf,0:ntf,0:npf))
      allocate (vepxf_p2(0:nrf,0:ntf,0:npf))
      allocate (vepyf_p2(0:nrf,0:ntf,0:npf))
      allocate (vepzf_p2(0:nrf,0:ntf,0:npf))
      allocate ( psif_p2(0:nrf,0:ntf,0:npf))
      allocate (alphf_p2(0:nrf,0:ntf,0:npf))
      allocate (bvxdf_p2(0:nrf,0:ntf,0:npf))
      allocate (bvydf_p2(0:nrf,0:ntf,0:npf))
      allocate (bvzdf_p2(0:nrf,0:ntf,0:npf))
      allocate (   rs_p2(0:ntf,0:npf))
      allocate (  psi_p2(0:nrg,0:ntg,0:npg))
      allocate ( alph_p2(0:nrg,0:ntg,0:npg))
      allocate ( bvxd_p2(0:nrg,0:ntg,0:npg))
      allocate ( bvyd_p2(0:nrg,0:ntg,0:npg))
      allocate ( bvzd_p2(0:nrg,0:ntg,0:npg))
      allocate (  axx_p2(0:nrg,0:ntg,0:npg))
      allocate (  axy_p2(0:nrg,0:ntg,0:npg))
      allocate (  axz_p2(0:nrg,0:ntg,0:npg))
      allocate (  ayy_p2(0:nrg,0:ntg,0:npg))
      allocate (  ayz_p2(0:nrg,0:ntg,0:npg))
      allocate (  azz_p2(0:nrg,0:ntg,0:npg))
      emd_p2=0.0d0;   vep_p2=0.0d0;    rs_p2=0.0d0;  wxspf_p2=0.0d0; wyspf_p2=0.0d0; wzspf_p2=0.0d0
      psi_p2=0.0d0;  alph_p2=0.0d0;  bvxd_p2=0.0d0;   bvyd_p2=0.0d0;  bvzd_p2=0.0d0
      axx_p2=0.0d0;   axy_p2=0.0d0;   axz_p2=0.0d0;    ayy_p2=0.0d0;   ayz_p2=0.0d0;   azz_p2=0.0d0

      !TODO: Change the following:
!      call IO_input_CF_irrot_BNS_cactus(impt, emd_p2,vep_p2,ome_p2,ber_p2,radi_p2,rs_p2, &
!         &    psi_p2,alph_p2,bvxd_p2,bvyd_p2,bvzd_p2,dir_path)

      call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt2.las",psi_p2,alph_p2,bvxd_p2,bvyd_p2,bvzd_p2)

      if(NS_shape.eq.'CO')  call IO_input_CF_flco_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",emd_p2,ome_p2,ber_p2,radi_p2)
      if(NS_shape.eq.'IR')  call IO_input_CF_flir_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",emd_p2,vep_p2,ome_p2,ber_p2,radi_p2)
      if(NS_shape.eq.'SP')  call IO_input_CF_flsp_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",emd_p2,vep_p2,wxspf_p2,wyspf_p2, &
                                                      &   wzspf_p2,ome_p2,ber_p2,radi_p2)

      call IO_input_CF_surf_export(trim(dir_path)//"/bnssur_3D_mpt2.las",rs_p2)

      call excurve_CF_gridpoint_export(alph_p2,bvxd_p2,bvyd_p2,bvzd_p2,  &
         &    axx_p2, axy_p2, axz_p2, ayy_p2, ayz_p2, azz_p2)

      call interpo_gr2fl_metric_CF_export(alph_p2, psi_p2, bvxd_p2, bvyd_p2, bvzd_p2, &
         &    alphf_p2, psif_p2, bvxdf_p2, bvydf_p2, bvzdf_p2, rs_p2)

      if((NS_shape.eq.'IR').or.(NS_shape.eq.'SP'))  call calc_gradvep_export(vep_p2,vepxf_p2,vepyf_p2,vepzf_p2,rs_p2)

      vepxf_p2 = - vepxf_p2
      vepyf_p2 = - vepyf_p2
      bvxdf_p2 = - bvxdf_p2
      bvydf_p2 = - bvydf_p2
      bvxd_p2  = - bvxd_p2
      bvyd_p2  = - bvyd_p2
    end if
    if (impt==3) then
      write(6,'(a13,i1,a8)') "Reading ARCP-", impt, " data..."
      allocate ( psi_p3(0:nrg,0:ntg,0:npg))
      allocate (alph_p3(0:nrg,0:ntg,0:npg))
      allocate (bvxd_p3(0:nrg,0:ntg,0:npg))
      allocate (bvyd_p3(0:nrg,0:ntg,0:npg))
      allocate (bvzd_p3(0:nrg,0:ntg,0:npg))
      allocate ( axx_p3(0:nrg,0:ntg,0:npg))
      allocate ( axy_p3(0:nrg,0:ntg,0:npg))
      allocate ( axz_p3(0:nrg,0:ntg,0:npg))
      allocate ( ayy_p3(0:nrg,0:ntg,0:npg))
      allocate ( ayz_p3(0:nrg,0:ntg,0:npg))
      allocate ( azz_p3(0:nrg,0:ntg,0:npg))
      psi_p3=0.0d0;  alph_p3=0.0d0;  bvxd_p3=0.0d0;  bvyd_p3=0.0d0; bvzd_p3=0.0d0
      axx_p3=0.0d0;   axy_p3=0.0d0;   axz_p3=0.0d0;   ayy_p3=0.0d0;  ayz_p3=0.0d0; azz_p3=0.0d0

      call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt3.las",psi_p3,alph_p3,bvxd_p3,bvyd_p3,bvzd_p3)
      
      call excurve_CF_gridpoint_export(alph_p3,bvxd_p3,bvyd_p3,bvzd_p3, &
         &    axx_p3, axy_p3, axz_p3, ayy_p3, ayz_p3, azz_p3)
    end if
  end do
  write(6,'(2e20.12)') emd_p1(0,0,0), emd_p1(58,0,0)
  write(6,'(3e20.12)') ome_p1, ber_p1, radi_p1
  write(6,'(e20.12)') dis_cm
!
  write(6,'(a56)', ADVANCE = "NO") "Give cartesian coordinates (x,y,z) separated by a space:"
  read(5,*) xcac,ycac,zcac
  write(6,'(a23,3e20.12)') "Point given wrt CACTUS:", xcac,ycac,zcac
  write(6,'(a38,2e20.12)') "Cocal radius scale in COCP-1, COCP-2 :", radi_p1, radi_p2
  write(6,'(a38,2e20.12)') "Cocal surface scale in COCP-1, COCP-2:", r_surf_p1, r_surf_p2
  xcoc = xcac/(radi_p1)
  ycoc = ycac/(radi_p1)
  zcoc = zcac/(radi_p1)
  write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc

  rrcm = sqrt(xcoc**2 + ycoc**2 + zcoc**2)
  write(6,*)  rrcm, rr3
  if (rrcm > rr3) then
!=>    call copy_from_mpatch_interpolation_utility(3)
    call copy_grid_parameter_from_mpt(3)
    call copy_grid_parameter_binary_excision_from_mpt(3)
    call copy_coordinate_grav_extended_from_mpt(3)
    call copy_coordinate_grav_phi_from_mpt(3)
    call copy_coordinate_grav_r_from_mpt(3)
    call copy_coordinate_grav_theta_from_mpt(3)
    call copy_def_binary_parameter_from_mpt(3)
    call copy_trigonometry_grav_phi_from_mpt(3)
    call copy_trigonometry_grav_theta_from_mpt(3)

    kij_parity = 1
    xc_p3 = xcoc
    yc_p3 = ycoc
    zc_p3 = zcoc
    write(6,'(a23,3e20.12)') "Point given wrt COCP-3:", xc_p3,yc_p3,zc_p3

    call interpo_gr2cgr_4th(psi_p3 , psica , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(alph_p3, alphca, xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(bvxd_p3, bvxdca, xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(bvyd_p3, bvydca, xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(bvzd_p3, bvzdca, xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(axx_p3 , axx   , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(axy_p3 , axy   , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(axz_p3 , axz   , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(ayy_p3 , ayy   , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(ayz_p3 , ayz   , xc_p3, yc_p3, zc_p3)
    call interpo_gr2cgr_4th(azz_p3 , azz   , xc_p3, yc_p3, zc_p3)
    psi4ca = psica**4
    vxu = 0.0d0
    vyu = 0.0d0
    vzu = 0.0d0
    emdca = 0.0d0
  else
    if (xcoc<=0.0d0) then
!=>      call copy_from_mpatch_interpolation_utility(1)
      call copy_grid_parameter_from_mpt(1)
      call copy_grid_parameter_binary_excision_from_mpt(1)
      call copy_coordinate_grav_extended_from_mpt(1)
      call copy_coordinate_grav_phi_from_mpt(1)
      call copy_coordinate_grav_r_from_mpt(1)
      call copy_coordinate_grav_theta_from_mpt(1)
      call copy_def_binary_parameter_from_mpt(1)
      call copy_trigonometry_grav_phi_from_mpt(1)
      call copy_trigonometry_grav_theta_from_mpt(1)

      kij_parity = 1
      xc_p1 = xcoc + dis_cm
      yc_p1 = ycoc
      zc_p1 = zcoc
      write(6,'(a23,3e20.12)') "Point given wrt COCP-1:", xc_p1,yc_p1,zc_p1  
      call interpo_gr2cgr_4th(psi_p1 , psica , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(alph_p1, alphca, xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(bvxd_p1, bvxdca, xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(bvyd_p1, bvydca, xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(bvzd_p1, bvzdca, xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(axx_p1 , axx   , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(axy_p1 , axy   , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(axz_p1 , axz   , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(ayy_p1 , ayy   , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(ayz_p1 , ayz   , xc_p1, yc_p1, zc_p1)
      call interpo_gr2cgr_4th(azz_p1 , azz   , xc_p1, yc_p1, zc_p1)
      call interpo_fl2cgr_4th_export(emd_p1  , emdca   , xc_p1, yc_p1, zc_p1, rs_p1)
!     call interpo_fl2cgr_4th_export(vep_p1  , vepca   , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(vepxf_p1, vepxfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(vepyf_p1, vepyfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(vepzf_p1, vepzfca , xc_p1, yc_p1, zc_p1, rs_p1)

      call interpo_fl2cgr_4th_export(wxspf_p1, wxspfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(wyspf_p1, wyspfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(wzspf_p1, wzspfca , xc_p1, yc_p1, zc_p1, rs_p1)

      call interpo_fl2cgr_4th_export(psif_p1 , psifca  , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(alphf_p1, alphfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(bvxdf_p1, bvxdfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(bvydf_p1, bvydfca , xc_p1, yc_p1, zc_p1, rs_p1)
      call interpo_fl2cgr_4th_export(bvzdf_p1, bvzdfca , xc_p1, yc_p1, zc_p1, rs_p1)
      bxcor = bvxdfca + ome_p1*(-ycoc)
      bycor = bvydfca + ome_p1*(xcoc)
      bzcor = bvzdfca
      lam_p1 = ber_p1 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca 
      psi4ca = psica**4
      psif4ca = psifca**4
      psifcacp= psifca**confpow

      if (dabs(emdca) > 1.0d-14) then
        vxu = ico*( ome_p1*(-ycoc) )                 +  &
        &     irr*( alphfca*vepxfca/psif4ca/lam_p1 ) +  &
        &     isp*( alphfca*(vepxfca/psif4ca + psifcacp*wxspfca)/lam_p1 )

        vyu = ico*( ome_p1*( xcoc) )                 +  &
        &     irr*( alphfca*vepyfca/psif4ca/lam_p1 ) +  &
        &     isp*( alphfca*(vepyfca/psif4ca + psifcacp*wyspfca)/lam_p1 ) 

        vzu = ico*( 0.0d0 )                          +  &
        &     irr*( alphfca*vepzfca/psif4ca/lam_p1 ) +  &
        &     isp*( alphfca*(vepzfca/psif4ca + psifcacp*wzspfca)/lam_p1 ) 
      else
        vxu=0.0d0;  vyu=0.0d0;  vzu=0.0d0
      end if
    else
!=>      call copy_from_mpatch_interpolation_utility(2)
      call copy_grid_parameter_from_mpt(2)
      call copy_grid_parameter_binary_excision_from_mpt(2)
      call copy_coordinate_grav_extended_from_mpt(2)
      call copy_coordinate_grav_phi_from_mpt(2)
      call copy_coordinate_grav_r_from_mpt(2)
      call copy_coordinate_grav_theta_from_mpt(2)
      call copy_def_binary_parameter_from_mpt(2)
      call copy_trigonometry_grav_phi_from_mpt(2)
      call copy_trigonometry_grav_theta_from_mpt(2)

      kij_parity = -1
      xc_p2 = -(xcoc - dis_cm)
      yc_p2 = -ycoc
      zc_p2 =  zcoc
      write(6,'(a23,3e20.12)') "Point given wrt COCP-2:", xc_p2,yc_p2,zc_p2
      call interpo_gr2cgr_4th(psi_p2 , psica , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(alph_p2, alphca, xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(bvxd_p2, bvxdca, xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(bvyd_p2, bvydca, xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(bvzd_p2, bvzdca, xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(axx_p2 , axx   , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(axy_p2 , axy   , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(axz_p2 , axz   , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(ayy_p2 , ayy   , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(ayz_p2 , ayz   , xc_p2, yc_p2, zc_p2)
      call interpo_gr2cgr_4th(azz_p2 , azz   , xc_p2, yc_p2, zc_p2)
      call interpo_fl2cgr_4th_export(emd_p2  , emdca   , xc_p2, yc_p2, zc_p2, rs_p2)
!     call interpo_fl2cgr_4th_export(vep_p2  , vepca   , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(vepxf_p2, vepxfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(vepyf_p2, vepyfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(vepzf_p2, vepzfca , xc_p2, yc_p2, zc_p2, rs_p2)

      call interpo_fl2cgr_4th_export(wxspf_p2, wxspfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(wyspf_p2, wyspfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(wzspf_p2, wzspfca , xc_p2, yc_p2, zc_p2, rs_p2)

      call interpo_fl2cgr_4th_export(psif_p2 , psifca  , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(alphf_p2, alphfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(bvxdf_p2, bvxdfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(bvydf_p2, bvydfca , xc_p2, yc_p2, zc_p2, rs_p2)
      call interpo_fl2cgr_4th_export(bvzdf_p2, bvzdfca , xc_p2, yc_p2, zc_p2, rs_p2)
      bxcor = bvxdfca + ome_p2*(-ycoc)
      bycor = bvydfca + ome_p2*(xcoc)
      bzcor = bvzdfca
      lam_p2 = ber_p2 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca
      psi4ca = psica**4 
      psif4ca = psifca**4
      psifcacp= psifca**confpow 

      if (dabs(emdca) > 1.0d-14) then
        vxu = ico*( ome_p2*(-ycoc) )                 +  &
        &     irr*( alphfca*vepxfca/psif4ca/lam_p2 ) +  &
        &     isp*( alphfca*(vepxfca/psif4ca + psifcacp*wxspfca)/lam_p2 ) 

        vyu = ico*( ome_p2*( xcoc) )                 +  &
        &     irr*( alphfca*vepyfca/psif4ca/lam_p2 ) +  &
        &     isp*( alphfca*(vepyfca/psif4ca + psifcacp*wyspfca)/lam_p2 ) 

        vzu = ico*( 0.0d0 )                          +  &
        &     irr*( alphfca*vepzfca/psif4ca/lam_p2 ) +  &
        &     isp*( alphfca*(vepzfca/psif4ca + psifcacp*wzspfca)/lam_p2 )  
      else
        vxu=0.0d0;  vyu=0.0d0;  vzu=0.0d0
      end if
    end if
  end if

  gxx = psi4ca
  gxy = 0.0d0
  gxz = 0.0d0
  gyy = psi4ca
  gyz = 0.0d0
  gzz = psi4ca

  kxx = psi4ca*axx/(radi_p1)
  kxy = psi4ca*axy/(radi_p1)
  kxz = kij_parity*psi4ca*axz/(radi_p1)
  kyy = psi4ca*ayy/(radi_p1)
  kyz = kij_parity*psi4ca*ayz/(radi_p1)
  kzz = psi4ca*azz/(radi_p1)

  call peos_q2hprho(emdca, hca, preca, rhoca, eneca)

  epsca = eneca/rhoca - 1.0d0

  write(6,'(a6,e20.12)') "psi  =", psica
  write(6,'(a6,e20.12)') "alph =", alphca
  write(6,'(a6,e20.12)') "bvxd =", bvxdca
  write(6,'(a6,e20.12)') "bvyd =", bvydca
  write(6,'(a6,e20.12)') "bvzd =", bvzdca
  write(6,'(a6,e20.12)') "Radi =", r_surf_p1*radi_p1
  write(6,'(a6,e20.12)') "Omeg =", ome_p1/radi_p1
  write(6,'(a6,e20.12)') "emd  =", emdca
  write(6,'(a6,e20.12)') "h    =", hca
  write(6,'(a6,e20.12)') "pre  =", preca
  write(6,'(a6,e20.12)') "rho  =", rhoca
  write(6,'(a6,e20.12)') "ene  =", eneca
  write(6,'(a6,e20.12)') "eps  =", epsca
  write(6,'(a6,e20.12)') "vepx =", vepxfca
  write(6,'(a6,e20.12)') "vepy =", vepyfca
  write(6,'(a6,e20.12)') "vepz =", vepzfca
!
  write(6,'(a18)') "Kij at gridpoints:"
  write(6,'(3e20.12)') kxx, kxy, kxz
  write(6,'(3e20.12)') kxy, kyy, kyz
  write(6,'(3e20.12)') kxz, kyz, kzz

  write(6,'(a13)') "v^i eulerian:"
  write(6,'(a6,e20.12)') "vxu  =", vxu
  write(6,'(a6,e20.12)') "vyu  =", vyu
  write(6,'(a6,e20.12)') "vzu  =", vzu

  write(6,'(a16)') "Deallocating...."
  deallocate(  emd_p1);  deallocate(  emd_p2); 
  deallocate(  vep_p1);  deallocate(  vep_p2); 
  deallocate(wxspf_p1);  deallocate(wxspf_p2);  
  deallocate(wyspf_p1);  deallocate(wyspf_p2);  
  deallocate(wzspf_p1);  deallocate(wzspf_p2);  
  deallocate(vepxf_p1);  deallocate(vepxf_p2);  
  deallocate(vepyf_p1);  deallocate(vepyf_p2);  
  deallocate(vepzf_p1);  deallocate(vepzf_p2);  
  deallocate( psif_p1);  deallocate( psif_p2);  
  deallocate(alphf_p1);  deallocate(alphf_p2);  
  deallocate(bvxdf_p1);  deallocate(bvxdf_p2);  
  deallocate(bvydf_p1);  deallocate(bvydf_p2);  
  deallocate(bvzdf_p1);  deallocate(bvzdf_p2);  
  deallocate(   rs_p1);  deallocate(   rs_p2);  
  deallocate(  psi_p1);  deallocate(  psi_p2);  deallocate(psi_p3)
  deallocate( alph_p1);  deallocate( alph_p2);  deallocate(alph_p3)
  deallocate( bvxd_p1);  deallocate( bvxd_p2);  deallocate(bvxd_p3)
  deallocate( bvyd_p1);  deallocate( bvyd_p2);  deallocate(bvyd_p3)
  deallocate( bvzd_p1);  deallocate( bvzd_p2);  deallocate(bvzd_p3)
  deallocate(  axx_p1);  deallocate(  axx_p2);  deallocate(axx_p3)
  deallocate(  axy_p1);  deallocate(  axy_p2);  deallocate(axy_p3)
  deallocate(  axz_p1);  deallocate(  axz_p2);  deallocate(axz_p3)
  deallocate(  ayy_p1);  deallocate(  ayy_p2);  deallocate(ayy_p3)
  deallocate(  ayz_p1);  deallocate(  ayz_p2);  deallocate(ayz_p3)
  deallocate(  azz_p1);  deallocate(  azz_p2);  deallocate(azz_p3)
!

END PROGRAM coc2cac
