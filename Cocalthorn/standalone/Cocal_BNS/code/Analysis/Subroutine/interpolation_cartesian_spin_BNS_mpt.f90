subroutine interpolation_cartesian_spin_BNS_mpt(impt)
  use phys_constant, only : long
  use grid_parameter, only :  nrf, ntf, npf, &
  &                          ntfeq, npfxzp, npfxzm, NS_shape
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use grid_parameter_cartesian, only : nx, ny, nz
  use grid_parameter_binary_excision, only : ex_rgmid, ex_radius 
  use coordinate_grav_xyz, only : x, y, z
  use def_binary_parameter, only : dis
  use def_metric
  use def_matter, only : emd, rs, vep, vepxf, vepyf, vepzf
  use def_velocity_rot
  use def_matter_parameter, only : confpow
  use def_metric_on_SFC_CF, only : psif
  use def_metric_cartesian
  use def_matter_cartesian
  use def_matter_velocity
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint
  use interface_modules_cartesian
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: fncca(:,:,:)
  real(long) :: par,dxvep,dyvep,dzvep,psif4,psifc
  integer :: ix, iy, iz, impt, impt_ex, ir, it, ip
!
  if(impt.eq.1) then
    impt_ex = 2
    par = 1.0d0
  end if
  if(impt.eq.2) then
    impt_ex = 1
    par = -1.0d0
  end if
!
  call interpolation_fillup_cartesian_parity_BNS_mpt(psi, psica,  impt, impt_ex, 1.0d0)
  call interpolation_fillup_cartesian_parity_BNS_mpt(alph,alphca, impt, impt_ex, 1.0d0)
  call interpolation_fillup_cartesian_parity_BNS_mpt(bvxd,bvxdca, impt, impt_ex,-1.0d0)
  call interpolation_fillup_cartesian_parity_BNS_mpt(bvyd,bvydca, impt, impt_ex,-1.0d0)
  call interpolation_fillup_cartesian_parity_BNS_mpt(bvzd,bvzdca, impt, impt_ex, 1.0d0)
!
  write(6,*) ' Spinning model.'
  write(6,*) " confpow = ", confpow
! origin
  vepxf(0,0:ntf,0:npf) = 0.0d0
  vepzf(0,0:ntf,0:npf) = 0.0d0

! x-axis
  vepxf(0:nrf,ntfeq,0)      = 0.0d0
  vepxf(0:nrf,ntfeq,npf)    = 0.0d0
  vepxf(0:nrf,ntfeq,npfxzm) = 0.0d0
  vepzf(0:nrf,ntfeq,0)      = 0.0d0
  vepzf(0:nrf,ntfeq,npf)    = 0.0d0
  vepzf(0:nrf,ntfeq,npfxzm) = 0.0d0

! xy plane
  do ir = 1, nrf
    do ip = 0, npf
      vepzf(ir,ntfeq,ip) = 0.0d0
    end do
  end do
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        psif4 = psif(ir,it,ip)**4
        psifc = psif(ir,it,ip)**confpow

        vxu(ir,it,ip) = (vepxf(ir,it,ip)/psif4 + psifc*wxspf(ir,it,ip))*par
        vyu(ir,it,ip) = (vepyf(ir,it,ip)/psif4 + psifc*wyspf(ir,it,ip))*par
        vzu(ir,it,ip) =  vepzf(ir,it,ip)/psif4 + psifc*wzspf(ir,it,ip)
      end do
    end do
  end do

  call interpolation_matter(emd,emdca)
  call interpolation_matter(vep,vepca)
  call interpolation_matter(vxu, vxca)
  call interpolation_matter(vyu, vyca)
  call interpolation_matter(vzu, vzca)
  call interpolation_matter(wxspf, wxspca)
  call interpolation_matter(wyspf, wyspca)
  call interpolation_matter(wzspf, wzspca)
!
end subroutine interpolation_cartesian_spin_BNS_mpt
