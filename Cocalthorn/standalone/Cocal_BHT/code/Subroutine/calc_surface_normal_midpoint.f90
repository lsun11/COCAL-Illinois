subroutine calc_surface_normal_midpoint(rs,rnx,rny,rnz,it,ip)
  use phys_constant, only :  long, nmpt
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use interface_flgrad_midpoint_type0_2Dsurf_sph
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: rs(:,:)
  real(long), intent(out) :: rnx, rny, rnz
  real(long) :: rnr, rnth, rnphi, dtrs, dprs, hrs, hrsinv
  real(long) :: erx, ery, erz, ethx, ethy, ethz, ephix, ephiy, ephiz
  integer    :: ir, it, ip
!
  call flgrad_midpoint_type0_2Dsurf_sph(rs,dtrs,dprs,it,ip)
  call interpo_linear_type0_2Dsurf(hrs,rs,it,ip)
  hrsinv = 1.0d0/hrs
!
  rnr = 1.0d0
  rnth = - hrsinv*dtrs
  rnphi = - hrsinv*hcosecthg(it)*dprs
  erx = hsinthg(it)*hcosphig(ip)
  ery = hsinthg(it)*hsinphig(ip)
  erz = hcosthg(it)
  ethx = hcosthg(it)*hcosphig(ip)
  ethy = hcosthg(it)*hsinphig(ip)
  ethz = - hsinthg(it)
  ephix = - hsinphig(ip)
  ephiy =   hcosphig(ip)
  ephiz = 0.0d0
  rnx = rnr*erx + rnth*ethx + rnphi*ephix
  rny = rnr*ery + rnth*ethy + rnphi*ephiy
  rnz = rnr*erz + rnth*ethz + rnphi*ephiz
!
end subroutine calc_surface_normal_midpoint
