#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!______________________________________________
include '../Include_file/include_modulefiles_RNS_CF_peos_plot.f90'
include '../Include_file/include_modulefiles_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_RNS_CF_peos_plot.f90'
include '../Include_file/include_interface_modulefiles_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_subroutines_RNS_CF_peos_plot.f90'
include '../Include_file/include_subroutines_analysis_RNS_CF_peos_plot.f90'
include '../Include_file/include_PEOS_modulefile.f90'
include '../Include_file/include_PEOS_subroutines.f90'
include '../Include_file/include_functions.f90'
!______________________________________________
!
!     ROTATING STAR COCAL ID to CACTUS
!______________________________________________
SUBROUTINE coc2cac_sub(CCTK_ARGUMENTS)
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
  use interface_IO_input_CF_star_export
  use interface_excurve_CF_gridpoint_export
  use interface_interpo_gr2fl_metric_CF_export

  use cctk
!  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer :: i, j, k, nx, ny, nz
  character(30) :: char1
  character*400 :: dir_path
  character*200 :: message
  integer :: dir_path_len
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc
  real(8) :: emdca, omefca, psica, alphca, bvxdca, bvydca, bvzdca, psi4ca, psif4ca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: axxca, axyca, axzca, ayyca, ayzca, azzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor, bvxdfca, bvydfca, bvzdfca, psifca, alphfca
  real(8) :: ome, ber, radi
!
  real(8), pointer :: emd(:,:,:), omef(:,:,:), rs(:,:)
  real(8), pointer :: psif(:,:,:), alphf(:,:,:), bvxdf(:,:,:), bvydf(:,:,:), bvzdf(:,:,:)
  real(8), pointer :: psi(:,:,:) , alph(:,:,:) , bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:)
  real(8), pointer :: axx(:,:,:), axy(:,:,:) , axz(:,:,:) , ayy(:,:,:) , ayz(:,:,:), azz(:,:,:)
!
  axxca=0.0d0; axyca=0.0d0; axzca=0.0d0; ayyca=0.0d0; ayzca=0.0d0; azzca=0.0d0
!
  !TODO remove this
  !dir_path='/home/astro/galeazzi/source/ET/wheeler/Cactus/repos/Cocal/standalone/Cocal/ID_BNS'
  !dir_path='.'
  call CCTK_FortranString(dir_path_len,dir_path_RNS_ID,dir_path)
!
! -- Read parameters
  call CCTK_INFO("Cocal: Reading parameters...")

! -- Read parameters
  call read_parameter_cactus(dir_path(1:dir_path_len))
  call peos_initialize_cactus(dir_path(1:dir_path_len))
  call grid_r
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  call grid_extended

! -- Allocate arrays
  call CCTK_INFO("Cocal: Allocating local arrays...")

  allocate (  emd(0:nrf,0:ntf,0:npf))
  allocate ( omef(0:nrf,0:ntf,0:npf))
  allocate ( psif(0:nrf,0:ntf,0:npf))
  allocate (alphf(0:nrf,0:ntf,0:npf))
  allocate (bvxdf(0:nrf,0:ntf,0:npf))
  allocate (bvydf(0:nrf,0:ntf,0:npf))
  allocate (bvzdf(0:nrf,0:ntf,0:npf))
  allocate (   rs(0:ntf,0:npf))
  allocate (  psi(0:nrg,0:ntg,0:npg))
  allocate ( alph(0:nrg,0:ntg,0:npg))
  allocate ( bvxd(0:nrg,0:ntg,0:npg))
  allocate ( bvyd(0:nrg,0:ntg,0:npg))
  allocate ( bvzd(0:nrg,0:ntg,0:npg))
  allocate (  axx(0:nrg,0:ntg,0:npg))
  allocate (  axy(0:nrg,0:ntg,0:npg))
  allocate (  axz(0:nrg,0:ntg,0:npg))
  allocate (  ayy(0:nrg,0:ntg,0:npg))
  allocate (  ayz(0:nrg,0:ntg,0:npg))
  allocate (  azz(0:nrg,0:ntg,0:npg))
  emd=0.0d0;  rs  =0.0d0;  omef=0.0d0
  psi=0.0d0;  alph=0.0d0;  bvxd=0.0d0;  bvyd=0.0d0;  bvzd=0.0d0
  axx=0.0d0;  axy =0.0d0;  axz =0.0d0;   ayy=0.0d0;   ayz=0.0d0;   azz=0.0d0

  call IO_input_CF_grav_export(dir_path(1:dir_path_len)//"/rnsgra_3D.las",psi,alph,bvxd,bvyd,bvzd)

  call IO_input_CF_star_export(dir_path(1:dir_path_len)//"/rnsflu_3D.las",emd,rs,omef,ome,ber,radi)

  call excurve_CF_gridpoint_export(alph,bvxd,bvyd,bvzd, &
     &    axx, axy, axz, ayy, ayz, azz)

  call interpo_gr2fl_metric_CF_export(alph, psi, bvxd, bvyd, bvzd, &
        &    alphf, psif, bvxdf, bvydf, bvzdf, rs)

  write(6,'(2e20.12)') emd(0,0,0), omef(0,0,0)
  write(6,'(3e20.12)') ome, ber, radi
!

!  call CCTK_INFO("Cocal: Looping over local cartesian grid:")
  write(6,*)"Cocal: Looping over local cartesian grid:"

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  write(6,*)"nx, ny, nz: ", nx, ny, nz
  
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        xcac = x(i,j,k)
        ycac = y(i,j,k)
        zcac = z(i,j,k)
!        write(6,*)' i, j, k, xcac, ycac, zcac', i, j, k, xcac, ycac, zcac
        xcoc = xcac/(radi)
        ycoc = ycac/(radi)
        zcoc = zcac/(radi)
!        write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc

        call interpo_gr2cgr_4th(psi , psica , xcoc, ycoc, zcoc)
        call interpo_gr2cgr_4th(alph, alphca, xcoc, ycoc, zcoc)
        call interpo_gr2cgr_4th(bvxd, bvxdca, xcoc, ycoc, zcoc)
        call interpo_gr2cgr_4th(bvyd, bvydca, xcoc, ycoc, zcoc)
        call interpo_gr2cgr_4th(bvzd, bvzdca, xcoc, ycoc, zcoc)
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
        call interpo_fl2cgr_4th_export(bvxdf, bvxdfca , xcoc, ycoc, zcoc, rs)
        call interpo_fl2cgr_4th_export(bvydf, bvydfca , xcoc, ycoc, zcoc, rs)
        call interpo_fl2cgr_4th_export(bvzdf, bvzdfca , xcoc, ycoc, zcoc, rs)

        bxcor = bvxdfca + omefca*(-ycoc)
        bycor = bvydfca + omefca*(xcoc)
        bzcor = bvzdfca
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

        gxx(i,j,k) = psi4ca
        gxy(i,j,k) = 0.0d0
        gxz(i,j,k) = 0.0d0
        gyy(i,j,k) = psi4ca
        gyz(i,j,k) = 0.0d0
        gzz(i,j,k) = psi4ca
      
        kxx(i,j,k) = psi4ca*axxca/(radi)
        kxy(i,j,k) = psi4ca*axyca/(radi)
        kxz(i,j,k) = psi4ca*axzca/(radi)
        kyy(i,j,k) = psi4ca*ayyca/(radi)
        kyz(i,j,k) = psi4ca*ayzca/(radi)
        kzz(i,j,k) = psi4ca*azzca/(radi)
      
        call peos_q2hprho(emdca, hca, preca, rhoca, eneca)
        
        alp(i,j,k) = alphca
        betax(i,j,k) = bvxdca
        betay(i,j,k) = bvydca
        betaz(i,j,k) = bvzdca

        rho(i,j,k) = rhoca
        press(i,j,k) = preca
        eps(i,j,k) = eneca/rhoca - 1.0d0
        vel(i,j,k,1) = vxu
        vel(i,j,k,2) = vyu
        vel(i,j,k,3) = vzu
        !TODO: Set w_lorentz in terms of 4-vel components

        !write(6,*) 'eps: ', eps(i,j,k), 'rho: ', rho(i,j,k), 'press: ', press(i,j,k), 'vel: ', vel(i,j,k,1), vel(i,j,k,2), vel(i,j,k,3)

!        write(6,'(a6,e20.12)') "psi  =", psica
!        write(6,'(a6,e20.12)') "alph =", alphca
!        write(6,'(a6,e20.12)') "bvxd =", bvxdca
!        write(6,'(a6,e20.12)') "bvyd =", bvydca
!        write(6,'(a6,e20.12)') "bvzd =", bvzdca
!        write(6,'(a6,e20.12)') "Radi =", r_surf*radi
!        write(6,'(a6,e20.12)') "Omeg =", ome/radi
!        write(6,'(a6,e20.12)') "emd  =", emdca
!        write(6,'(a6,e20.12)') "h    =", hca
!        write(6,'(a6,e20.12)') "pre  =", preca
!        write(6,'(a6,e20.12)') "rho  =", rhoca
!        write(6,'(a6,e20.12)') "ene  =", eneca
!      !
!        write(6,'(a18)') "kij at gridpoints:"
!        write(6,'(3e20.12)') coc_kxx, coc_kxy, coc_kxz
!        write(6,'(3e20.12)') coc_kxy, coc_kyy, coc_kyz
!        write(6,'(3e20.12)') coc_kxz, coc_kyz, coc_kzz
!      
!        write(6,'(a13)') "v^i eulerian:"
!        write(6,'(a6,e20.12)') "vxu  =", vxu
!        write(6,'(a6,e20.12)') "vyu  =", vyu
!        write(6,'(a6,e20.12)') "vzu  =", vzu

      end do
    end do
  end do

  if( CCTK_EQUALS(initial_shift,"zero") ) then
    betax = 0.0d0
    betay = 0.0d0
    betaz = 0.0d0
  end if

  write(6,*)"Cocal: Finished looping over local cartesian grid:"
  write(6,'(a16)') "Deallocating...."
  deallocate(  emd);  
  deallocate( psif);  
  deallocate(alphf);  
  deallocate(bvxdf);  
  deallocate(bvydf); 
  deallocate(bvzdf);  
  deallocate(   rs);  
  deallocate(  psi); 
  deallocate( alph);  
  deallocate( bvxd); 
  deallocate( bvyd);  
  deallocate( bvzd);  
  deallocate(  axx);
  deallocate(  axy);  
  deallocate(  axz);  
  deallocate(  ayy);  
  deallocate(  ayz); 
  deallocate(  azz); 

!TODO: Output to standout ome, ber and radi

! write (message, '("Orbital angular Velocity: ", f24.16)') ome
! call CCTK_INFO(message)
! write (message, '("Constant of Euler integral: ", f24.16)') ber
! call CCTK_INFO(message)
! write (message, '("Maximum radius of neutron star: ", f24.16)') radi
! call CCTK_INFO(message)

END SUBROUTINE coc2cac_sub

