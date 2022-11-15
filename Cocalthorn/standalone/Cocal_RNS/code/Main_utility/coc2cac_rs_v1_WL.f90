!______________________________________________
include '../Include_file/include_modulefiles_RNS_WL_peos_plot.f90'
include '../Include_file/include_modulefiles_analysis_RNS_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_RNS_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_RNS_WL_peos_plot.f90'
include '../Include_file/include_subroutines_RNS_WL_peos_plot.f90'
include '../Include_file/include_subroutines_analysis_RNS_WL_peos_plot.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!    ROTATING STAR COCAL Waveless ID to CACTUS
!______________________________________________
PROGRAM coc2cac
!
  use phys_constant
!  use def_matter_parameter
  use grid_parameter
  use interface_modules_cartesian
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use interface_IO_input_CF_grav_export
  use interface_IO_input_WL_grav_export_hij
  use interface_IO_input_grav_export_Aij
  use interface_IO_input_CF_star_export
  use interface_invhij_WL_export
  use interface_index_vec_down2up_export
  use interface_interpo_gr2fl_metric_CF_export
  implicit none
  character(30) :: char1
  character*400 :: dir_path
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc
  real(8) :: emdca, omefca, psica, alphca, bvxdca, bvydca, bvzdca, psi4ca, psif4ca
  real(8) :: hxxdca, hxydca,hxzdca,hyydca,hyzdca,hzzdca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: axxca, axyca, axzca, ayyca, ayzca, azzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor, bvxufca, bvyufca, bvzufca, psifca, alphfca
  real(8) :: gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz
  real(8) :: ome, ber, radi
!
  real(8), pointer :: emd(:,:,:),   omef(:,:,:),  rs(:,:)
  real(8), pointer :: psif(:,:,:),  alphf(:,:,:), bvxuf(:,:,:), bvyuf(:,:,:), bvzuf(:,:,:)
  real(8), pointer :: psi(:,:,:) ,  alph(:,:,:) , bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:)
  real(8), pointer :: bvxu(:,:,:) , bvyu(:,:,:) , bvzu(:,:,:)
  real(8), pointer :: hxxd(:,:,:) , hxyd(:,:,:) , hxzd(:,:,:) , hyyd(:,:,:) , hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:) , hxyu(:,:,:) , hxzu(:,:,:) , hyyu(:,:,:) , hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: axx(:,:,:),   axy(:,:,:)  , axz(:,:,:)  , ayy(:,:,:)  , ayz(:,:,:) , azz(:,:,:)
!
  gxx=0.0d0; gxy=0.0d0; gxz=0.0d0; gyy=0.0d0; gyz=0.0d0; gzz=0.0d0
  kxx=0.0d0; kxy=0.0d0; kxz=0.0d0; kyy=0.0d0; kyz=0.0d0; kzz=0.0d0
  axxca=0.0d0; axyca=0.0d0; axzca=0.0d0; ayyca=0.0d0; ayzca=0.0d0; azzca=0.0d0

  !TODO remove this
  !dir_path="/home/astro/mundim/tmp/ET_2014_05_wheeler/Cactus/repos/Cocal/standalone/Cocal/ID_BNS"
  !dir_path="../../standalone/Cocal/ID_BNS"
  dir_path='.'

! -- Read parameters
  call read_parameter_cactus(dir_path)
  call peos_initialize_cactus(dir_path)
  call grid_r
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call grid_extended
!
!    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
!  rr3 = 0.7d0*(rgout - rgmid)
!  dis_cm = dis

  allocate (  emd(0:nrf,0:ntf,0:npf))
  allocate ( omef(0:nrf,0:ntf,0:npf))
  allocate ( psif(0:nrf,0:ntf,0:npf))
  allocate (alphf(0:nrf,0:ntf,0:npf))
  allocate (bvxuf(0:nrf,0:ntf,0:npf))
  allocate (bvyuf(0:nrf,0:ntf,0:npf))
  allocate (bvzuf(0:nrf,0:ntf,0:npf))
  allocate (   rs(0:ntf,0:npf))
  allocate (  psi(0:nrg,0:ntg,0:npg))
  allocate ( alph(0:nrg,0:ntg,0:npg))
  allocate ( bvxd(0:nrg,0:ntg,0:npg))
  allocate ( bvyd(0:nrg,0:ntg,0:npg))
  allocate ( bvzd(0:nrg,0:ntg,0:npg))
  allocate ( bvxu(0:nrg,0:ntg,0:npg))
  allocate ( bvyu(0:nrg,0:ntg,0:npg))
  allocate ( bvzu(0:nrg,0:ntg,0:npg))
  allocate ( hxxd(0:nrg,0:ntg,0:npg))
  allocate ( hxyd(0:nrg,0:ntg,0:npg))
  allocate ( hxzd(0:nrg,0:ntg,0:npg))
  allocate ( hyyd(0:nrg,0:ntg,0:npg))
  allocate ( hyzd(0:nrg,0:ntg,0:npg))
  allocate ( hzzd(0:nrg,0:ntg,0:npg))
  allocate ( hxxu(0:nrg,0:ntg,0:npg))
  allocate ( hxyu(0:nrg,0:ntg,0:npg))
  allocate ( hxzu(0:nrg,0:ntg,0:npg))
  allocate ( hyyu(0:nrg,0:ntg,0:npg))
  allocate ( hyzu(0:nrg,0:ntg,0:npg))
  allocate ( hzzu(0:nrg,0:ntg,0:npg))
  allocate (  axx(0:nrg,0:ntg,0:npg))
  allocate (  axy(0:nrg,0:ntg,0:npg))
  allocate (  axz(0:nrg,0:ntg,0:npg))
  allocate (  ayy(0:nrg,0:ntg,0:npg))
  allocate (  ayz(0:nrg,0:ntg,0:npg))
  allocate (  azz(0:nrg,0:ntg,0:npg))
  emd=0.0d0;  rs  =0.0d0;  omef=0.0d0
  psi=0.0d0;  alph=0.0d0;  bvxd=0.0d0;  bvyd=0.0d0;  bvzd=0.0d0
  bvxu=0.0d0; bvyu=0.0d0;  bvzu=0.0d0
  axx=0.0d0;  axy =0.0d0;  axz =0.0d0;   ayy=0.0d0;   ayz=0.0d0;   azz=0.0d0
  hxxd=0.0d0; hxyd=0.0d0;  hxzd=0.0d0;  hyyd=0.0d0;  hyzd=0.0d0;  hzzd=0.0d0;
  hxxu=0.0d0; hxyu=0.0d0;  hxzu=0.0d0;  hyyu=0.0d0;  hyzu=0.0d0;  hzzu=0.0d0;

  call IO_input_CF_grav_export(trim(dir_path)//"/rnsgra_3D.las",psi,alph,bvxd,bvyd,bvzd)

  call IO_input_WL_grav_export_hij(trim(dir_path)//"/rnsgra_hij_3D.las",hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)

  call IO_input_grav_export_Aij(trim(dir_path)//"/rnsgra_Aij_3D.las",axx,axy,axz,ayy,ayz,azz)

  call IO_input_CF_star_export(trim(dir_path)//"/rnsflu_3D.las",emd,rs,omef,ome,ber,radi)

  call invhij_WL_export(hxxd,hxyd,hxzd,hyyd,hyzd,hzzd,hxxu,hxyu,hxzu,hyyu,hyzu,hzzu)

  call index_vec_down2up_export(hxxu,hxyu,hxzu,hyyu,hyzu,hzzu,bvxu,bvyu,bvzu,bvxd,bvyd,bvzd)

!  call excurve_CF_gridpoint_export(alph,bvxd,bvyd,bvzd, & 
!     &    axx, axy, axz, ayy, ayz, azz)

  call interpo_gr2fl_metric_CF_export(alph, psi, bvxu, bvyu, bvzu, &
        &    alphf, psif, bvxuf, bvyuf, bvzuf, rs)


  write(6,'(2e20.12)') emd(0,0,0), omef(0,0,0)
  write(6,'(3e20.12)') ome, ber, radi
!
  write(6,'(a56)', ADVANCE = "NO") "Give cartesian coordinates (x,y,z) separated by a space:"
  read(5,*) xcac,ycac,zcac
  write(6,'(a23,3e20.12)') "Point given wrt CACTUS:", xcac,ycac,zcac
  write(6,'(a20,1e20.12)') "Cocal radius scale :", radi
  xcoc = xcac/(radi)
  ycoc = ycac/(radi)
  zcoc = zcac/(radi)
  write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc


  call interpo_gr2cgr_4th(psi , psica , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(alph, alphca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hxxd, hxxdca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hxyd, hxydca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hxzd, hxzdca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hyyd, hyydca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hyzd, hyzdca, xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(hzzd, hzzdca, xcoc, ycoc, zcoc)
  
  call interpo_gr2cgr_4th(axx , axxca , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(axy , axyca , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(axz , axzca , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(ayy , ayyca , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(ayz , ayzca , xcoc, ycoc, zcoc)
  call interpo_gr2cgr_4th(azz , azzca , xcoc, ycoc, zcoc)

  call interpo_fl2cgr_4th_export(emd  , emdca   , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(omef , omefca  , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(psif , psifca  , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(alphf, alphfca , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(bvxuf, bvxufca , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(bvyuf, bvyufca , xcoc, ycoc, zcoc, rs)
  call interpo_fl2cgr_4th_export(bvzuf, bvzufca , xcoc, ycoc, zcoc, rs)

  bxcor = bvxufca + omefca*(-ycoc)
  bycor = bvyufca + omefca*(xcoc)
  bzcor = bvzufca
  psi4ca = psica**4
  psif4ca = psifca**4

  if (dabs(emdca) > 1.0d-14) then
    vxu = bxcor/alphfca 
    vyu = bycor/alphfca
    vzu = bzcor/alphfca
  else
    emdca=0.0d0
    vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
  end if

  gxx = psi4ca*(1.0d0+hxxdca)
  gxy = psi4ca*(0.0d0+hxydca)
  gxz = psi4ca*(0.0d0+hxzdca)
  gyy = psi4ca*(1.0d0+hyydca)
  gyz = psi4ca*(0.0d0+hyzdca)
  gzz = psi4ca*(1.0d0+hzzdca)

  kxx = psi4ca*axxca/(radi)
  kxy = psi4ca*axyca/(radi)
  kxz = psi4ca*axzca/(radi)
  kyy = psi4ca*ayyca/(radi)
  kyz = psi4ca*ayzca/(radi)
  kzz = psi4ca*azzca/(radi)

  call peos_q2hprho(emdca, hca, preca, rhoca, eneca)

  epsca = eneca/rhoca - 1.0d0

  write(6,'(a6,e20.12)') "psi  =", psica
  write(6,'(a6,e20.12)') "alph =", alphca
  write(6,'(a6,e20.12)') "Radi =", radi
  write(6,'(a6,e20.12)') "Omeg =", ome/radi
  write(6,'(a6,e20.12)') "emd  =", emdca
  write(6,'(a6,e20.12)') "h    =", hca
  write(6,'(a6,e20.12)') "pre  =", preca
  write(6,'(a6,e20.12)') "rho  =", rhoca
  write(6,'(a6,e20.12)') "ene  =", eneca
  write(6,'(a6,e20.12)') "eps  =", epsca
!
  write(6,'(a18)') "gij at gridpoints:"
  write(6,'(3e20.12)') gxx, gxy, gxz
  write(6,'(3e20.12)') gxy, gyy, gyz
  write(6,'(3e20.12)') gxz, gyz, gzz

  write(6,'(a18)') "Kij at gridpoints:"
  write(6,'(3e20.12)') kxx, kxy, kxz
  write(6,'(3e20.12)') kxy, kyy, kyz
  write(6,'(3e20.12)') kxz, kyz, kzz

  write(6,'(a13)') "v^i Eulerian:"
  write(6,'(a6,e20.12)') "vxu  =", vxu
  write(6,'(a6,e20.12)') "vyu  =", vyu
  write(6,'(a6,e20.12)') "vzu  =", vzu

  write(6,'(a16)') "Deallocating...."
  deallocate(  emd);  deallocate( omef);  deallocate( psif);  deallocate(alphf);    
  deallocate(bvxuf);  deallocate(bvyuf);  deallocate(bvzuf);  deallocate(   rs);    
  deallocate(  psi);  deallocate( alph);  deallocate( bvxd);  deallocate( bvyd);  
  deallocate( bvzd);  deallocate( bvxu);  deallocate( bvyu);  deallocate( bvzu);
  deallocate( hxxd);  deallocate( hxyd);  deallocate( hxzd);  deallocate( hyyd);  
  deallocate( hyzd);  deallocate( hzzd);  deallocate( hxxu);  deallocate( hxyu);  
  deallocate( hxzu);  deallocate( hyyu);  deallocate( hyzu);  deallocate( hzzu);  
  deallocate(  axx);  deallocate(  axy);  deallocate(  axz);  deallocate(  ayy);  
  deallocate(  ayz);  deallocate(  azz);  
!
END PROGRAM coc2cac
