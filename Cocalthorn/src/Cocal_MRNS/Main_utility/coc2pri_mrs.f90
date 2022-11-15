!______________________________________________
include '../Include_file/include_modulefiles_MRNS.f90'
include '../Include_file/include_modulefiles_analysis_plot_MRNS.f90'
include '../Include_file/include_interface_modulefiles_plot_MRNS.f90'
include '../Include_file/include_interface_modulefiles_analysis_plot_MRNS.f90'
include '../Include_file/include_subroutines_plot_MRNS.f90'
include '../Include_file/include_subroutines_analysis_plot_MRNS.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!
!_____________________________________________________________
!
!    MAGNETIZED ROTATING STAR COCAL Waveless ID to CACTUS
!_____________________________________________________________
SUBROUTINE coc2pri_mrs(dir_path,Xb,Yb,Zb,Nx,Ny,Nz,&
     gbtt,gbtx,gbty,gbtz,gbxx,gbxy,gbxz,gbyy,gbyz,gbzz,&
     gbxx_t,gbxy_t,gbxz_t,gbyy_t,gbyz_t,gbzz_t,&
     Pb,rhob,vbx,vby,vbz,Ax,Ay,Az,iAB,path_len)
! If iAB=1 it exports magnetic potential Ai. Otherwise it exports magnetic field Bi.

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
  use interface_IO_input_CF_star_export
  use interface_invhij_WL_export
  use interface_index_vec_down2up_export
  use interface_interpo_gr2fl_metric_CF_export
  use interface_IO_input_grav_export_Ai
  use interface_IO_input_grav_export_Faraday
  use interface_IO_input_star4ve_export
  implicit none
  character(30) :: char1
  character*400 :: dir_path
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc
  real(8) :: emdca,  omefca, psica,  alphca, psi4ca, psif4ca
  real(8) :: bvxdca, bvydca, bvzdca, bvxuca, bvyuca, bvzuca
  real(8) :: hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca
  real(8) :: hxxuca, hxyuca, hxzuca, hyyuca, hyzuca, hzzuca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: kxxca, kxyca, kxzca, kyyca, kyzca, kzzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor, bvxufca, bvyufca, bvzufca, psifca, alphfca
  real(8) :: gxx1, gxy1, gxz1, gyy1, gyz1, gzz1, kxx1, kxy1, kxz1, kyy1, kyz1, kzz1
  real(8) :: ome, ber, radi
  real(8) :: va1, vaxd1, vayd1, vazd1, fxd1, fyd1, fzd1, fxyd1, fxzd1, fyzd1
  real(8) :: vaca, vaxdca, vaydca, vazdca, fxdca, fydca, fzdca, fxydca, fxzdca, fyzdca
  real(8) :: utfca, uxfca, uyfca, uzfca
!
  real(8), pointer :: emd(:,:,:), omef(:,:,:), rs(:,:)
  real(8), pointer :: utf(:,:,:)  ,  uxf(:,:,:) ,  uyf(:,:,:) ,  uzf(:,:,:)
  real(8), pointer :: psif(:,:,:) , alphf(:,:,:), bvxuf(:,:,:), bvyuf(:,:,:), bvzuf(:,:,:)
  real(8), pointer :: psi(:,:,:)  , alph(:,:,:)
  real(8), pointer :: bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:) , bvxu(:,:,:) , bvyu(:,:,:), bvzu(:,:,:)
  real(8), pointer :: hxxd(:,:,:) , hxyd(:,:,:) , hxzd(:,:,:) , hyyd(:,:,:) , hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:) , hxyu(:,:,:) , hxzu(:,:,:) , hyyu(:,:,:) , hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: kxx(:,:,:)  , kxy(:,:,:)  , kxz(:,:,:)  , kyy(:,:,:)  , kyz(:,:,:) , kzz(:,:,:)
  real(8), pointer ::  va(:,:,:)  , vaxd(:,:,:) , vayd(:,:,:) , vazd(:,:,:)
  real(8), pointer :: fxd(:,:,:)  ,  fyd(:,:,:) ,  fzd(:,:,:) , fxyd(:,:,:) , fxzd(:,:,:), fyzd(:,:,:)

!
  integer :: Nx,Ny,Nz,inx,iny,inz,path_len, iAB
  real(8) :: gbtt(Nx,Ny,Nz),gbtx(Nx,Ny,Nz),gbty(Nx,Ny,Nz),gbtz(Nx,Ny,Nz)
  real(8) :: gbxx(Nx,Ny,Nz),gbxy(Nx,Ny,Nz),gbxz(Nx,Ny,Nz),gbyy(Nx,Ny,Nz),gbyz(Nx,Ny,Nz),gbzz(Nx,Ny,Nz)
  real(8) :: gbxx_t(Nx,Ny,Nz),gbxy_t(Nx,Ny,Nz),gbxz_t(Nx,Ny,Nz),gbyy_t(Nx,Ny,Nz),gbyz_t(Nx,Ny,Nz),gbzz_t(Nx,Ny,Nz)
  real(8) :: Pb(Nx,Ny,Nz),rhob(Nx,Ny,Nz),vbx(Nx,Ny,Nz),vby(Nx,Ny,Nz),vbz(Nx,Ny,Nz),Xb(Nx),Yb(Ny),Zb(Nz)
  real(8) :: beta_x,beta_y,beta_z
  real(8) :: Ax(Nx,Ny,Nz),Ay(Nx,Ny,Nz),Az(Nx,Ny,Nz)

  gxx1=0.0d0; gxy1=0.0d0; gxz1=0.0d0; gyy1=0.0d0; gyz1=0.0d0; gzz1=0.0d0
  kxx1=0.0d0; kxy1=0.0d0; kxz1=0.0d0; kyy1=0.0d0; kyz1=0.0d0; kzz1=0.0d0
  kxxca=0.0d0; kxyca=0.0d0; kxzca=0.0d0; kyyca=0.0d0; kyzca=0.0d0; kzzca=0.0d0
  vaca=0.0d0;  vaxdca=0.0d0;  vaydca=0.0d0;  vazdca=0.0d0;
  fxdca=0.0d0; fydca=0.0d0; fzdca=0.0d0; fxydca=0.0d0; fxzdca=0.0d0; fyzdca=0.0d0;
  fxd1=0.0d0;  fyd1=0.0d0;  fzd1=0.0d0;  fxyd1=0.0d0;  fxzd1=0.0d0; fyzd1=0.0d0;
  va1=0.0d0;  vaxd1=0.0d0; vayd1=0.0d0; vazd1=0.0d0; 
  utfca=0.0d0;  uxfca=0.0d0;  uyfca=0.0d0;  uzfca=0.0d0;

  dir_path=dir_path(1:path_len)
  write(*,*) "In coc2pri_mrs: dir_path=",dir_path
  if (iAB .eq. 1) then
    write(*,*) "Exporting Ai"
  else
    write(*,*) "Exporting Bi" 
  end if

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
  allocate (  utf(0:nrf,0:ntf,0:npf))
  allocate (  uxf(0:nrf,0:ntf,0:npf))
  allocate (  uyf(0:nrf,0:ntf,0:npf))
  allocate (  uzf(0:nrf,0:ntf,0:npf))

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
  if (iAB .eq. 1) then
    allocate (   va(0:nrg,0:ntg,0:npg))
    allocate ( vaxd(0:nrg,0:ntg,0:npg))
    allocate ( vayd(0:nrg,0:ntg,0:npg))
    allocate ( vazd(0:nrg,0:ntg,0:npg))
    va=0.0d0;    vaxd=0.0d0;  vayd=0.0d0;   vazd=0.0d0;
  else
    allocate (  fxd(0:nrg,0:ntg,0:npg))
    allocate (  fyd(0:nrg,0:ntg,0:npg))
    allocate (  fzd(0:nrg,0:ntg,0:npg))
    allocate ( fxyd(0:nrg,0:ntg,0:npg))
    allocate ( fxzd(0:nrg,0:ntg,0:npg))
    allocate ( fyzd(0:nrg,0:ntg,0:npg))
    fxd=0.0d0;   fyd=0.0d0;   fzd=0.0d0;    fxyd=0.0d0;   fxzd=0.0d0;  fyzd=0.0d0;
  end if

  emd=0.0d0;   rs  =0.0d0;  omef=0.0d0
  utf=0.0d0;   uxf=0.0d0;   uyf=0.0d0;    uzf=0.0d0;
  psif=0.0d0;  alphf=0.0d0; bvxuf=0.0d0;  bvyuf=0.0d0;  bvzuf=0.0d0
  psi=0.0d0;   alph=0.0d0;  bvxu=0.0d0;   bvyu=0.0d0;  bvzu=0.0d0
  bvxd=0.0d0;  bvyd=0.0d0;  bvzd=0.0d0; 
  kxx=0.0d0;   kxy=0.0d0;   kxz=0.0d0;    kyy=0.0d0;    kyz=0.0d0;   kzz=0.0d0
  hxxd=0.0d0;  hxyd=0.0d0;  hxzd=0.0d0;   hyyd=0.0d0;   hyzd=0.0d0;  hzzd=0.0d0;
  hxxu=0.0d0;  hxyu=0.0d0;  hxzu=0.0d0;   hyyu=0.0d0;   hyzu=0.0d0;  hzzu=0.0d0;

  call IO_input_CF_grav_export(trim(dir_path)//"/rnsgra_3D.las",psi,alph,bvxd,bvyd,bvzd)

  call IO_input_WL_grav_export_hij(trim(dir_path)//"/rnsgra_hij_3D.las",hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)

  call IO_input_grav_export_Kij(trim(dir_path)//"/rnsgra_Kij_3D.las",kxx,kxy,kxz,kyy,kyz,kzz)

  if (iAB .eq. 1) then
    call IO_input_grav_export_Ai(trim(dir_path)//"/rnsEMF_3D.las",va,vaxd,vayd,vazd)
  else
    call IO_input_grav_export_Faraday(trim(dir_path)//"/rnsEMF_faraday_3D.las",fxd,fyd,fzd,fxyd,fxzd,fyzd)
  end if

  call IO_input_CF_star_export(trim(dir_path)//"/rnsflu_3D.las",emd,rs,omef,ome,ber,radi)

  call IO_input_star4ve_export(trim(dir_path)//"/rns4ve_3D.las",utf,uxf,uyf,uzf)

  call invhij_WL_export(hxxd,hxyd,hxzd,hyyd,hyzd,hzzd,hxxu,hxyu,hxzu,hyyu,hyzu,hzzu)

  call index_vec_down2up_export(hxxu,hxyu,hxzu,hyyu,hyzu,hzzu,bvxu,bvyu,bvzu,bvxd,bvyd,bvzd)

  call interpo_gr2fl_metric_CF_export(alph, psi, bvxu, bvyu, bvzu, &
        &    alphf, psif, bvxuf, bvyuf, bvzuf, rs)
!
!
  write(6,'(2e20.12)') emd(0,0,0), omef(0,0,0)
  write(6,'(3e20.12)') ome, ber, radi
!
!
!
  do inz=1,Nz
     zcac=Zb(inz)
     do iny=1,Ny
        ycac=Yb(iny)
        do inx=1,Nx
           xcac=Xb(inx)

           xcoc = xcac/(radi)
           ycoc = ycac/(radi)
           zcoc = zcac/(radi)
           !  write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc

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

           call interpo_fl2cgr_4th_export(emd  , emdca   , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(omef , omefca  , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(psif , psifca  , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(alphf, alphfca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(bvxuf, bvxufca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(bvyuf, bvyufca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(bvzuf, bvzufca , xcoc, ycoc, zcoc, rs)

           call interpo_fl2cgr_4th_export(utf, utfca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(uxf, uxfca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(uyf, uyfca , xcoc, ycoc, zcoc, rs)
           call interpo_fl2cgr_4th_export(uzf, uzfca , xcoc, ycoc, zcoc, rs)
           
           if (iAB .eq. 1) then
             call interpo_gr2cgr_4th(  va,   vaca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(vaxd, vaxdca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(vayd, vaydca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(vazd, vazdca, xcoc, ycoc, zcoc)
           else
             call interpo_gr2cgr_4th(fxd,   fxdca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(fyd,   fydca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(fzd,   fzdca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(fxyd, fxydca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(fxzd, fxzdca, xcoc, ycoc, zcoc)
             call interpo_gr2cgr_4th(fyzd, fyzdca, xcoc, ycoc, zcoc)
           end if

           bxcor = bvxufca + omefca*(-ycoc)
           bycor = bvyufca + omefca*(xcoc)
           bzcor = bvzufca
           psi4ca = psica**4
           psif4ca = psifca**4

           if (dabs(emdca) > 1.0d-14) then
!              vxu = bxcor/alphfca 
!              vyu = bycor/alphfca
!              vzu = bzcor/alphfca
             vxu = (uxfca/utfca + bvxufca)/alphfca
             vyu = (uyfca/utfca + bvyufca)/alphfca
             vzu = (uzfca/utfca + bvzufca)/alphfca
           else
              emdca=0.0d0
              vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
           end if

           ! Setting 4-metric and assigning Kij to gbij_t
           gbxx(inx,iny,inz) = psi4ca*(1.0d0+hxxdca)
           gbxy(inx,iny,inz) = psi4ca*(      hxydca)
           gbxz(inx,iny,inz) = psi4ca*(      hxzdca)
           gbyy(inx,iny,inz) = psi4ca*(1.0d0+hyydca)
           gbyz(inx,iny,inz) = psi4ca*(      hyzdca)
           gbzz(inx,iny,inz) = psi4ca*(1.0d0+hzzdca)
           
           gbxx_t(inx,iny,inz) = psi4ca*kxxca/(radi)
           gbxy_t(inx,iny,inz) = psi4ca*kxyca/(radi)
           gbxz_t(inx,iny,inz) = psi4ca*kxzca/(radi)
           gbyy_t(inx,iny,inz) = psi4ca*kyyca/(radi)
           gbyz_t(inx,iny,inz) = psi4ca*kyzca/(radi)
           gbzz_t(inx,iny,inz) = psi4ca*kzzca/(radi)

           beta_x = psi4ca*bvxdca     ! b_x lower index
           beta_y = psi4ca*bvydca
           beta_z = psi4ca*bvzdca
           
           gbtt(inx,iny,inz)=-alphca*alphca + &
                psi4ca*(bvxdca*bvxuca+bvydca*bvyuca+bvzdca*bvzuca)
           gbtx(inx,iny,inz)=beta_x
           gbty(inx,iny,inz)=beta_y
           gbtz(inx,iny,inz)=beta_z
           
           call peos_q2hprho(emdca, hca, preca, rhoca, eneca)

           epsca = eneca/rhoca - 1.0d0

           ! Setting pressure
           Pb(inx,iny,inz)=preca
           ! Setting rest-mass density
           rhob(inx,iny,inz)=rhoca
           ! Setting Eulerian velocity
           vbx(inx,iny,inz)=vxu
           vby(inx,iny,inz)=vyu
           vbz(inx,iny,inz)=vzu

           if (iAB .eq. 1) then
             Ax(inx,iny,inz) = vaxdca           ! A_i
             Ay(inx,iny,inz) = vaydca
             Az(inx,iny,inz) = vazdca
           else
             Ax(inx,iny,inz) =  fyzdca/radi
             Ay(inx,iny,inz) = -fxzdca/radi
             Az(inx,iny,inz) =  fxydca/radi      ! B_i   
           end if
           !   write(6,'(a6,e20.12)') "psi  =", psica
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
     end do
  end do
!  
!           
  write(6,'(a16)') "Deallocating...."
  deallocate(  emd);  deallocate( omef);  deallocate( psif);  deallocate(alphf);
  deallocate(  utf);  deallocate(  uxf);  deallocate(  uyf);  deallocate(  uzf);
  deallocate(bvxuf);  deallocate(bvyuf);  deallocate(bvzuf);  deallocate(   rs);
  deallocate(  psi);  deallocate( alph);  deallocate( bvxd);  deallocate( bvyd);
  deallocate( bvzd);  deallocate( bvxu);  deallocate( bvyu);  deallocate( bvzu);
  deallocate( hxxd);  deallocate( hxyd);  deallocate( hxzd);  deallocate( hyyd);
  deallocate( hyzd);  deallocate( hzzd);  deallocate( hxxu);  deallocate( hxyu);
  deallocate( hxzu);  deallocate( hyyu);  deallocate( hyzu);  deallocate( hzzu);
  deallocate(  kxx);  deallocate(  kxy);  deallocate(  kxz);  deallocate(  kyy);
  deallocate(  kyz);  deallocate(  kzz);
  if (iAB .eq. 1) then
    deallocate(   va);  deallocate( vaxd);  deallocate( vayd);  deallocate( vazd);
  else
    deallocate(  fxd);  deallocate(  fyd);  deallocate(  fzd);
    deallocate( fxyd);  deallocate( fxzd);  deallocate( fyzd);
  end if
!
END SUBROUTINE coc2pri_mrs
