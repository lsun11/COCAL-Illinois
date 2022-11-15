subroutine interpolation_cartesian
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_metric
  use def_matter, only : emd, rs, omef
  use def_metric_cartesian
  use def_matter_cartesian
  use def_matter_velocity
  use interface_modules_cartesian
  implicit none
  real(long) :: xxxx, yyyy, omega
  integer :: ir, it, ip
!
  call interpolation_metric(psi,psica)
  call interpolation_metric(alph,alphca)
  call interpolation_metric(bvxd,bvxdca)
  call interpolation_metric(bvyd,bvydca)
  call interpolation_metric(bvzd,bvzdca)
!
  call allocate_matter_velocity
!
  write(6,*) ' Flow field is circular.'
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        omega = omef(ir,it,ip)
        xxxx = rg(ir)*rs(it,ip)*sinthg(it)*cosphig(ip)
        yyyy = rg(ir)*rs(it,ip)*sinthg(it)*sinphig(ip)
        vxu(ir,it,ip) = - omega*yyyy
        vyu(ir,it,ip) =   omega*xxxx
        vzu(ir,it,ip) = 0.0d0
      end do
    end do
  end do
!
  call interpolation_matter(emd,emdca)
  call interpolation_matter(omef,omeca)
  call interpolation_matter(vxu,vxca)
  call interpolation_matter(vyu,vyca)
  call interpolation_matter(vzu,vzca)
!
end subroutine interpolation_cartesian
