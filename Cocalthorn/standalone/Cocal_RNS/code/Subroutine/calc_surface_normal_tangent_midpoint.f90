subroutine calc_surface_normal_tangent_midpoint(rs,nv,tv,it,ip)
  use phys_constant, only :  long, nmpt
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use interface_flgrad_midpoint_type0_2Dsurf_sph
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: rs(:,:), nv(:), tv(:)
!  real(long), intent(out) :: rnx, rny, rnz
  real(long) :: rnr, rnth, rnphi, dtrs, dprs, hrs, hrsinv, nvnorm, tvnorm
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

  nvnorm = sqrt(rnr*rnr + rnth*rnth + rnphi*rnphi)

  erx = hsinthg(it)*hcosphig(ip)
  ery = hsinthg(it)*hsinphig(ip)
  erz = hcosthg(it)
  ethx = hcosthg(it)*hcosphig(ip)
  ethy = hcosthg(it)*hsinphig(ip)
  ethz = - hsinthg(it)
  ephix = - hsinphig(ip)
  ephiy =   hcosphig(ip)
  ephiz = 0.0d0
  nv(1) = (rnr*erx + rnth*ethx + rnphi*ephix)/nvnorm
  nv(2) = (rnr*ery + rnth*ethy + rnphi*ephiy)/nvnorm
  nv(3) = (rnr*erz + rnth*ethz + rnphi*ephiz)/nvnorm
!
! As tangent we get dxdphi. Cartesian coordinates.
  rnr = dprs
  rnth = 0.0d0 
  rnphi = hrs*hsinthg(it)  

!  tvnorm = sqrt(rnr*rnr + rnth*rnth + rnphi*rnphi)
  tvnorm = 1.0d0

  tv(1) = (rnr*erx + rnth*ethx + rnphi*ephix)/tvnorm
  tv(2) = (rnr*ery + rnth*ethy + rnphi*ephiy)/tvnorm
  tv(3) = (rnr*erz + rnth*ethz + rnphi*ephiz)/tvnorm
!
end subroutine calc_surface_normal_tangent_midpoint
