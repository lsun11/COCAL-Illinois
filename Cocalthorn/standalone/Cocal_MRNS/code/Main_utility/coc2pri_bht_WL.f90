!______________________________________________
include '../Include_file/include_modulefiles_BHT_WL_peos_plot.f90'
include '../Include_file/include_modulefiles_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_BHT_WL_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_subroutines_BHT_WL_peos_plot.f90'
include '../Include_file/include_subroutines_analysis_BHT_WL_peos_plot.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!
!___________________________________________________
!
!    Black hole-torus COCAL Waveless ID to CACTUS 
!___________________________________________________

SUBROUTINE coc2pri_bht_wl(dir_path,Xb,Yb,Zb,ntot,&
     lapseb,Psib,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
     gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
     Pb,rhob,vbx,vby,vbz,path_len)


  use phys_constant
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
  use interface_IO_input_grav_export_Kij
  use interface_IO_input_matter_BHT_export
  use interface_invhij_WL_export
  use interface_index_vec_down2up_export
  implicit none
  character(30) :: char1
  character*400 :: dir_path
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc, rcoc
  real(8) :: emdgca, omegca
  real(8) :: psica,  alphca, psi4ca
  real(8) :: bvxdca, bvydca, bvzdca, bvxuca, bvyuca, bvzuca
  real(8) :: hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca
  real(8) :: hxxuca, hxyuca, hxzuca, hyyuca, hyzuca, hzzuca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: kxxca, kxyca, kxzca, kyyca, kyzca, kzzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor
  real(8) :: gxx1, gxy1, gxz1, gyy1, gyz1, gzz1, kxx1, kxy1, kxz1, kyy1, kyz1, kzz1
  real(8) :: ome, ber, radi, rexc
!
  real(8), pointer :: emdg(:,:,:) , omeg(:,:,:)
  real(8), pointer :: psi(:,:,:)  , alph(:,:,:)
  real(8), pointer :: bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:) , bvxu(:,:,:) , bvyu(:,:,:), bvzu(:,:,:)
  real(8), pointer :: hxxd(:,:,:) , hxyd(:,:,:) , hxzd(:,:,:) , hyyd(:,:,:) , hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:) , hxyu(:,:,:) , hxzu(:,:,:) , hyyu(:,:,:) , hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: kxx(:,:,:),   kxy(:,:,:)  , kxz(:,:,:)  , kyy(:,:,:)  , kyz(:,:,:) , kzz(:,:,:)
!
  integer :: ntot,inx,iny,inz,path_len
  real(8) :: lapseb(ntot),Psib(ntot)
  real(8) :: gbtt(ntot),gbtx(ntot),gbty(ntot),gbtz(ntot)
  real(8) :: gbxx(ntot),gbxy(ntot),gbxz(ntot),gbyy(ntot),gbyz(ntot),gbzz(ntot)
  real(8) :: gbxx_t(ntot),gbxy_t(ntot),gbxz_t(ntot),gbyy_t(ntot),gbyz_t(ntot),gbzz_t(ntot)
  real(8) :: Pb(ntot),rhob(ntot),vbx(ntot),vby(ntot),vbz(ntot),Xb(ntot),Yb(ntot),Zb(ntot)
  real(8) :: beta_x,beta_y,beta_z, detg

  gxx1=0.0d0; gxy1=0.0d0; gxz1=0.0d0; gyy1=0.0d0; gyz1=0.0d0; gzz1=0.0d0
  kxx1=0.0d0; kxy1=0.0d0; kxz1=0.0d0; kyy1=0.0d0; kyz1=0.0d0; kzz1=0.0d0
  kxxca=0.0d0; kxyca=0.0d0; kxzca=0.0d0; kyyca=0.0d0; kyzca=0.0d0; kzzca=0.0d0

  dir_path=dir_path(1:path_len)
  write(*,*) "In coc2pri_bht_WL: dir_path=",dir_path

! -- Read parameters
  call read_parameter_cactus(dir_path)
  call read_bht_parameter_cactus(dir_path)
  call calc_bht_excision_radius
  call peos_initialize_cactus(dir_path)
  call grid_r_bht('eBH')
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call grid_extended
!
  rexc = rg(0)
  write(6,'(a14,1p,1e23.15)') "Excision at r=", rexc
!    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
!  rr3 = 0.7d0*(rgout - rgmid)
!  dis_cm = dis

  allocate ( emdg(0:nrg,0:ntg,0:npg))
  allocate ( omeg(0:nrg,0:ntg,0:npg))
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
  allocate (  kxx(0:nrg,0:ntg,0:npg))
  allocate (  kxy(0:nrg,0:ntg,0:npg))
  allocate (  kxz(0:nrg,0:ntg,0:npg))
  allocate (  kyy(0:nrg,0:ntg,0:npg))
  allocate (  kyz(0:nrg,0:ntg,0:npg))
  allocate (  kzz(0:nrg,0:ntg,0:npg))
  emdg=0.0d0; omeg=0.0d0
  psi=0.0d0;  alph=0.0d0
  bvxd=0.0d0; bvyd=0.0d0;  bvzd=0.0d0;  bvxu=0.0d0;  bvyu=0.0d0;  bvzu=0.0d0
  kxx=0.0d0;  kxy =0.0d0;  kxz =0.0d0;   kyy=0.0d0;   kyz=0.0d0;   kzz=0.0d0
  hxxd=0.0d0; hxyd=0.0d0;  hxzd=0.0d0;  hyyd=0.0d0;  hyzd=0.0d0;  hzzd=0.0d0;
  hxxu=0.0d0; hxyu=0.0d0;  hxzu=0.0d0;  hyyu=0.0d0;  hyzu=0.0d0;  hzzu=0.0d0;


  call IO_input_CF_grav_export(trim(dir_path)//'/'//'rnsgra_3D.las',psi,alph,bvxd,bvyd,bvzd)

  call IO_input_WL_grav_export_hij(trim(dir_path)//'/'//'rnsgra_hij_3D.las',hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)

  call IO_input_grav_export_Kij(trim(dir_path)//'/'//'rnsgra_Kij_3D.las',kxx,kxy,kxz,kyy,kyz,kzz)

  call IO_input_matter_BHT_export(trim(dir_path)//"/rnsflu_3D.las",emdg,omeg,ome,ber,radi)

  call invhij_WL_export(hxxd,hxyd,hxzd,hyyd,hyzd,hzzd,hxxu,hxyu,hxzu,hyyu,hyzu,hzzu)

  call index_vec_down2up_export(hxxu,hxyu,hxzu,hyyu,hyzu,hzzu,bvxu,bvyu,bvzu,bvxd,bvyd,bvzd)


!  write(6,'(2e20.12)') emd(0,0,0), omef(0,0,0)
  write(6,'(3e20.12,1I4.2)') ome, ber, radi, ntot
!  
  
        do inx=1,ntot

           xcac=Xb(inx)
           ycac=Yb(inx)
           zcac=Zb(inx)

           xcoc = xcac/(radi)
           ycoc = ycac/(radi)
           zcoc = zcac/(radi)
           rcoc = dsqrt(xcoc*xcoc + ycoc*ycoc + zcoc*zcoc)

!           if (rcoc < 1.1d0*rexc)  cycle
           if (rcoc < rexc)  cycle
!           write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc
!           write(6,'(a23,3e20.12)') "Point given:", xcac,ycac,zcac

           call interpo_gr2cgr_4th(psi , psica , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(alph, alphca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvxd, bvxdca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvyd, bvydca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvzd, bvzdca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvxu, bvxuca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvyu, bvyuca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(bvzu, bvzuca, xcoc, ycoc, zcoc)

           call interpo_gr2cgr_4th(hxxd, hxxdca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(hxyd, hxydca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(hxzd, hxzdca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(hyyd, hyydca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(hyzd, hyzdca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(hzzd, hzzdca, xcoc, ycoc, zcoc)

           call interpo_gr2cgr_4th(kxx , kxxca , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(kxy , kxyca , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(kxz , kxzca , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(kyy , kyyca , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(kyz , kyzca , xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(kzz , kzzca , xcoc, ycoc, zcoc)

           call interpo_gr2cgr_4th(emdg, emdgca, xcoc, ycoc, zcoc)
           call interpo_gr2cgr_4th(omeg, omegca, xcoc, ycoc, zcoc)

           bxcor = bvxuca + omegca*(-ycoc)
           bycor = bvyuca + omegca*(xcoc)
           bzcor = bvzuca
           psi4ca = psica**4

           if (dabs(emdgca) > 1.0d-14) then
              vxu = bxcor/alphca 
              vyu = bycor/alphca
              vzu = bzcor/alphca
           else
              emdgca=0.0d0
              vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
           end if

           ! Setting 3-metric and assigning Kij to gbij_t
           gbxx(inx) = psi4ca*(1.0d0+hxxdca)
           gbxy(inx) = psi4ca*(      hxydca)
           gbxz(inx) = psi4ca*(      hxzdca)
           gbyy(inx) = psi4ca*(1.0d0+hyydca)
           gbyz(inx) = psi4ca*(      hyzdca)
           gbzz(inx) = psi4ca*(1.0d0+hzzdca)

           gbxx_t(inx) = psi4ca*kxxca/(radi)
           gbxy_t(inx) = psi4ca*kxyca/(radi)
           gbxz_t(inx) = psi4ca*kxzca/(radi)
           gbyy_t(inx) = psi4ca*kyyca/(radi)
           gbyz_t(inx) = psi4ca*kyzca/(radi)
           gbzz_t(inx) = psi4ca*kzzca/(radi)

           beta_x = psi4ca*bvxdca     ! b_x lower index
           beta_y = psi4ca*bvydca
           beta_z = psi4ca*bvzdca

           
           gbtx(inx)=beta_x
           gbty(inx)=beta_y
           gbtz(inx)=beta_z

           lapseb(inx) = alphca       ! set alpha
           
           detg = -gbxz(inx)*gbxz(inx)*gbyy(inx) + 2.0*gbxy(inx)*gbxz(inx)*gbyz(inx)&
                - gbxx(inx)*gbyz(inx)*gbyz(inx) - gbxy(inx)*gbxy(inx)*gbzz(inx)&
                + gbxx(inx)*gbyy(inx)*gbzz(inx);
           Psib(inx)   = detg**(1.d0/12.d0)        ! set BSSN conformal factor
           

           
           call peos_q2hprho(emdgca, hca, preca, rhoca, eneca)

           epsca = eneca/rhoca - 1.0d0

           ! Setting pressure
           Pb(inx)=preca
           ! Setting rest-mass density
           rhob(inx)=rhoca
           ! Setting Eulerian velocity
           vbx(inx)=vxu
           vby(inx)=vyu
           vbz(inx)=vzu


           !   write(*,'(6e20.12)') Xcac, gbxx(inx), gbxy(inx), gbxz(inx), gbyy(inx), gbyz(inx), gbzz(inx)
           !   write(6,'(a6,e20.12)') "alph =", alphca
           !   write(6,'(a6,e20.12)') "bvxd =", bvxdca
           !   write(6,'(a6,e20.12)') "bvyd =", bvydca
           !   write(6,'(a6,e20.12)') "bvzd =", bvzdca
           !   write(6,'(a6,e20.12)') "Radi =", radi
           !   write(6,'(a6,e20.12)') "Omeg =", ome/radi
           !   write(6,'(a6,e20.12)') "emd  =", emdca
           !   write(6,'(a6,e20.12)') "h    =", hca
           !   write(6,'(a6,e20.12)') "pre  =", preca
           !   write(6,'(a6,e20.12)') "rho  =", rhoca
           !   write(6,'(a6,e20.12)') "ene  =", eneca
           !   write(6,'(a6,e20.12)') "eps  =", epsca
           ! !
           !   write(6,'(a13)') "v^i eulerian:"
           !   write(6,'(a6,e20.12)') "vxu  =", vxu
           !   write(6,'(a6,e20.12)') "vyu  =", vyu
           !   write(6,'(a6,e20.12)') "vzu  =", vzu
        end do
!     end do
!  end do
  
           
  write(6,'(a16)') "Deallocating...."
  deallocate( emdg);  deallocate( omeg);
  deallocate(  psi);  deallocate( alph);  deallocate( bvxd);  deallocate( bvyd);
  deallocate( bvzd);  deallocate( bvxu);  deallocate( bvyu);  deallocate( bvzu);
  deallocate( hxxd);  deallocate( hxyd);  deallocate( hxzd);  deallocate( hyyd);
  deallocate( hyzd);  deallocate( hzzd);  deallocate( hxxu);  deallocate( hxyu);
  deallocate( hxzu);  deallocate( hyyu);  deallocate( hyzu);  deallocate( hzzu);
  deallocate(  kxx);  deallocate(  kxy);  deallocate(  kxz);  deallocate(  kyy);
  deallocate(  kyz);  deallocate(  kzz);

END SUBROUTINE coc2pri_bht_wl
