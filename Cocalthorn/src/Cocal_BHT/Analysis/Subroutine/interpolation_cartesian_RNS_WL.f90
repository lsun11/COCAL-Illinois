subroutine interpolation_cartesian_RNS_WL
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_metric
  use def_metric_hij
  use def_matter, only : emd, rs, omef
  use def_matter_parameter, only : ome
  use def_metric_cartesian
  use def_metric_hij_cartesian
  use def_matter_cartesian
  use def_matter_velocity
  use interface_modules_cartesian
  implicit none
  real(long) :: xxxx, yyyy, omew
  integer :: ir, it, ip
!
  call interpolation_metric(psi,psica)
  call interpolation_metric(alph,alphca)
  call interpolation_metric(bvxd,bvxdca)
  call interpolation_metric(bvyd,bvydca)
  call interpolation_metric(bvzd,bvzdca)
  call interpolation_metric(hxxd,hxxdca)
  call interpolation_metric(hxyd,hxydca)
  call interpolation_metric(hxzd,hxzdca)
  call interpolation_metric(hyyd,hyydca)
  call interpolation_metric(hyzd,hyzdca)
  call interpolation_metric(hzzd,hzzdca)
!
call invhij
call calc_shift_down2up
  call interpolation_metric(bvxu,bvxuca)
  call interpolation_metric(bvyu,bvyuca)
  call interpolation_metric(bvzu,bvzuca)
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        xxxx = rg(ir)*rs(it,ip)*sinthg(it)*cosphig(ip)
        yyyy = rg(ir)*rs(it,ip)*sinthg(it)*sinphig(ip)
        omew = omef(ir,it,ip)
        vxu(ir,it,ip) = - omew*yyyy
        vyu(ir,it,ip) =   omew*xxxx
        vzu(ir,it,ip) = 0.0d0
      end do
    end do
  end do
!
  call interpolation_matter(emd,emdca)
  call interpolation_matter(vxu,vxca)
  call interpolation_matter(vyu,vyca)
  call interpolation_matter(vzu,vzca)
!
end subroutine interpolation_cartesian_RNS_WL
