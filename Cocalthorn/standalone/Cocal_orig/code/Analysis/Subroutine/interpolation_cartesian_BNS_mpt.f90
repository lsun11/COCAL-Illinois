subroutine interpolation_cartesian_BNS_mpt(impt)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use grid_parameter_cartesian, only : nx, ny, nz
  use grid_parameter_binary_excision, only : ex_rgmid, ex_radius 
  use coordinate_grav_xyz, only : x, y, z
  use def_binary_parameter, only : dis
  use def_metric
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome
  use def_metric_cartesian
  use def_matter_cartesian
  use def_matter_velocity
  use interface_modules_cartesian
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: fncca(:,:,:)
  real(long) :: xxxx, yyyy, par
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
  write(6,*) ' Co-rotational model.'
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        xxxx = rg(ir)*rs(it,ip)*sinthg(it)*cosphig(ip)
        yyyy = rg(ir)*rs(it,ip)*sinthg(it)*sinphig(ip)
        vxu(ir,it,ip) = - ome*yyyy*par
        vyu(ir,it,ip) =   ome*(xxxx - dis )*par
        vzu(ir,it,ip) = 0.0d0
      end do
    end do
  end do
!
  call interpolation_matter(emd,emdca)
  call interpolation_matter(vxu, vxca)
  call interpolation_matter(vyu, vyca)
  call interpolation_matter(vzu, vzca)
!
end subroutine interpolation_cartesian_BNS_mpt
