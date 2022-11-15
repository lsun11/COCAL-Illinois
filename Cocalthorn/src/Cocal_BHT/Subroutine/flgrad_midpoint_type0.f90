subroutine flgrad_midpoint_type0(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
  use phys_constant, only : long
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) :: gr1, gr2, gr3, hrs, hrsinv
  real(long) :: dfncdr, dfncdth, dfncdphi, drsdth, drsdphi
  integer, intent(in) :: irf, itf, ipf
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives.
!
  dfncdr   = 0.25d0 &
  &     *(fnc(irf,itf  ,ipf  ) - fnc(irf-1,itf  ,ipf  ) &
  &     + fnc(irf,itf-1,ipf  ) - fnc(irf-1,itf-1,ipf  ) &
  &     + fnc(irf,itf  ,ipf-1) - fnc(irf-1,itf  ,ipf-1) &
  &     + fnc(irf,itf-1,ipf-1) - fnc(irf-1,itf-1,ipf-1))*drginv(irf)
  dfncdth  = 0.25d0 &
  &     *(fnc(irf  ,itf,ipf  ) - fnc(irf  ,itf-1,ipf  ) &
  &     + fnc(irf-1,itf,ipf  ) - fnc(irf-1,itf-1,ipf  ) &
  &     + fnc(irf  ,itf,ipf-1) - fnc(irf  ,itf-1,ipf-1) &
  &     + fnc(irf-1,itf,ipf-1) - fnc(irf-1,itf-1,ipf-1))*dthginv
  dfncdphi = 0.25d0  &
  &     *(fnc(irf  ,itf  ,ipf) - fnc(irf  ,itf  ,ipf-1) &
  &     + fnc(irf-1,itf  ,ipf) - fnc(irf-1,itf  ,ipf-1) &
  &     + fnc(irf  ,itf-1,ipf) - fnc(irf  ,itf-1,ipf-1) &
  &     + fnc(irf-1,itf-1,ipf) - fnc(irf-1,itf-1,ipf-1))*dphiginv
  drsdth  = 0.5d0 &
  &       *(rs(itf,ipf  ) - rs(itf-1,ipf  ) &
  &       + rs(itf,ipf-1) - rs(itf-1,ipf-1))*dthginv
  drsdphi = 0.5d0 &
  &       *(rs(itf  ,ipf) - rs(itf  ,ipf-1) &
  &       + rs(itf-1,ipf) - rs(itf-1,ipf-1))*dphiginv
  hrs = 0.25d0 &
  &       *(rs(itf,ipf  ) + rs(itf-1,ipf  ) &
  &       + rs(itf,ipf-1) + rs(itf-1,ipf-1))
  hrsinv = 1.0d0/hrs
!
! --- To cartesian component.
!
  gr1  = dfncdr*hrsinv
  gr2  = dfncdth*hrginv(irf)*hrsinv &
  &    - drsdth*hrsinv*gr1
  gr3  = dfncdphi*hrginv(irf)*hrsinv*hcosecthg(itf) &
  &    - drsdphi*hrsinv*hcosecthg(itf)*gr1
!
  dfdx = gr1 * hsinthg(itf) * hcosphig(ipf) &
  &    + gr2 * hcosthg(itf) * hcosphig(ipf) &
  &    - gr3 * hsinphig(ipf)
  dfdy = gr1 * hsinthg(itf) * hsinphig(ipf) &
  &    + gr2 * hcosthg(itf) * hsinphig(ipf) &
  &    + gr3 * hcosphig(ipf)
  dfdz = gr1 * hcosthg(itf)  &
  &    - gr2 * hsinthg(itf)
!
end subroutine flgrad_midpoint_type0
