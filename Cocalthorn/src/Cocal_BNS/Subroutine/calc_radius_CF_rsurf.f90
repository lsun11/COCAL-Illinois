subroutine calc_radius_CF_rsurf
  use phys_constant, only  : long, pi
  use grid_parameter, only : r_surf, nrf, ntf, npf, &
  &                          ntfpolp, ntfeq, ntfxy, npfxzp, npfyzp
  use def_metric_on_SFC_CF
  use def_matter, only : rs
  use def_matter_parameter, only : radi
  use def_quantities, only : proper_radius_x,proper_radius_y,proper_radius_z, &
  &                          coord_radius_x, coord_radius_y, coord_radius_z
  use make_array_1d
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_radial_int_fluid
  implicit none
  integer             :: it, ip, ir
  real(long)          :: radius
  real(long), pointer :: sou(:)
!
  call alloc_array1d(sou, 0, nrf)
!
  it = ntfxy ; ip = npfxzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_x = r_surf*radi*radius
  coord_radius_x  = r_surf*radi*rs(it,ip)
!
  it = ntfxy ; ip = npfyzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_y = r_surf*radi*radius
  coord_radius_y  = r_surf*radi*rs(it,ip)
!
  it = ntfpolp ; ip = npfxzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_z = r_surf*radi*radius
  coord_radius_z  = r_surf*radi*rs(it,ip)
!
  deallocate(sou)
end subroutine calc_radius_CF_rsurf
