subroutine flgrad_midpoint_type0_parallel(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
  use phys_constant, only : long
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer :: fnc(:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) :: gr1, gr2, gr3, hrs, hrsinv
  real(long) :: dfncdr, dfncdth, dfncdphi, drsdth, drsdphi
  integer, intent(in) :: irf, itf, ipf
!
! --- Compute the xy gradient of a function at fixed theta gridpoint.
! --- The gradient is evaluated at mid points.
!
! --- r, phi derivatives. Theta is fixed
!
  dfncdr   = 0.5d0 &
  &     *(fnc(irf,ipf  ) - fnc(irf-1,ipf  ) &
  &     + fnc(irf,ipf-1) - fnc(irf-1,ipf-1))*drginv(irf)
  dfncdphi = 0.5d0  &
  &     *(fnc(irf  ,ipf) - fnc(irf  ,ipf-1) &
  &     + fnc(irf-1,ipf) - fnc(irf-1,ipf-1))*dphiginv
  
  drsdphi = (rs(itf,ipf) - rs(itf,ipf-1))*dphiginv

  hrs = 0.5d0*(rs(itf,ipf) + rs(itf,ipf-1))

  hrsinv = 1.0d0/hrs
!
! --- To cartesian component.
!
  gr1  = dfncdr*hrsinv
  gr2  = 0.0d0
  gr3  = dfncdphi*hrginv(irf)*hrsinv*cosecthg(itf) &
  &    - drsdphi*hrsinv*cosecthg(itf)*gr1
!
  dfdx = gr1 * sinthg(itf) * hcosphig(ipf) &
  &    + gr2 * costhg(itf) * hcosphig(ipf) &
  &    - gr3 * hsinphig(ipf)
  dfdy = gr1 * sinthg(itf) * hsinphig(ipf) &
  &    + gr2 * costhg(itf) * hsinphig(ipf) &
  &    + gr3 * hcosphig(ipf)
  dfdz = gr1 * costhg(itf)  &
  &    - gr2 * sinthg(itf)
!
end subroutine flgrad_midpoint_type0_parallel
