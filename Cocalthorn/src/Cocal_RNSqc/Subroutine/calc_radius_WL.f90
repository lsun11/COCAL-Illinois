subroutine calc_radius_WL
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfpolp, ntfeq, ntfxy, npfxzp, npfyzp
  use def_metric, only : psi
  use def_metric_hij, only : hxxd, hyyd, hzzd
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
  real(long), pointer :: psif(:,:,:)
  real(long), pointer :: hxxdf(:,:,:), hyydf(:,:,:), hzzdf(:,:,:)
!
  call alloc_array1d(sou, 0, nrf)
  call alloc_array3d(psif, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hxxdf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hyydf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(hzzdf, 0, nrf, 0, ntf, 0, npf)
!
  call interpo_gr2fl(psi, psif)
  call interpo_gr2fl(hxxd, hxxdf)
  call interpo_gr2fl(hyyd, hyydf)
  call interpo_gr2fl(hzzd, hzzdf)
!
  it = ntfxy ; ip = npfxzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2*sqrt(1.0d0 + hxxdf(0:nrf,it,ip))
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_x = radi*radius
  coord_radius_x  = radi*rs(it,ip)
!
  it = ntfxy ; ip = npfyzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2*sqrt(1.0d0 + hyydf(0:nrf,it,ip))
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_y = radi*radius
  coord_radius_y  = radi*rs(it,ip)
!
  it = ntfpolp ; ip = npfxzp
  sou(0:nrf) = psif(0:nrf,it,ip)**2*sqrt(1.0d0 + hzzdf(0:nrf,it,ip))
  call radial_int_fluid(sou,radius,it,ip)
  proper_radius_z = radi*radius
  coord_radius_z  = radi*rs(it,ip)
!
  deallocate(sou)
  deallocate(psif)
end subroutine calc_radius_WL
