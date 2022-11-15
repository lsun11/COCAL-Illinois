subroutine source_vep_surface_CF_peos(vpot_v,surp)
  use phys_constant, only :  long, nmpt
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC_CF
  use def_matter, only : rs, vep
  use def_velocity_potential
  use interface_interpo_linear_surface_type0
  use interface_flgrad_midpoint_surface_type0
  use interface_calc_surface_normal_midpoint
  implicit none
  real(long), pointer :: surp(:,:), vpot_v(:,:,:)
  real(long) :: dxvpot_v, dyvpot_v, dzvpot_v, rnx, rny, rnz
  integer    :: ir, it, ip
!
  ir = nrf
  do ip = 1, npf
    do it = 1, ntf
      call flgrad_midpoint_surface_type0(vpot_v,dxvpot_v,dyvpot_v,dzvpot_v,it,ip)
      call calc_surface_normal_midpoint(rs,rnx,rny,rnz,it,ip)
      surp(it,ip) = dxvpot_v*rnx + dyvpot_v*rny + dzvpot_v*rnz
    end do
  end do
!
end subroutine source_vep_surface_CF_peos
