subroutine flgrad_type0(fnc,fnc_rs,dfdx,dfdy,dfdz,irf,itf,ipf)
  use phys_constant, only : long
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  implicit none
  real(long), intent(in)  :: fnc(0:1,0:1,0:1), fnc_rs(0:1,0:1)
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
  &        *(fnc(1,1,1) - fnc(0,1,1) &
  &        + fnc(1,0,1) - fnc(0,0,1) &
  &        + fnc(1,1,0) - fnc(0,1,0) &
  &        + fnc(1,0,0) - fnc(0,0,0))*drginv(irf)
  dfncdth  = 0.25d0 &
  &        *(fnc(1,1,1) - fnc(1,0,1) &
  &        + fnc(0,1,1) - fnc(0,0,1) &
  &        + fnc(1,1,0) - fnc(1,0,0) &
  &        + fnc(0,1,0) - fnc(0,0,0))*dthginv
  dfncdphi = 0.25d0  &
  &        *(fnc(1,1,1) - fnc(1,1,0) &
  &        + fnc(0,1,1) - fnc(0,1,0) &
  &        + fnc(1,0,1) - fnc(1,0,0) &
  &        + fnc(0,0,1) - fnc(0,0,0))*dphiginv
  drsdth  = 0.5d0 &
  &       *(fnc_rs(1,1) - fnc_rs(0,1) &
  &       + fnc_rs(1,0) - fnc_rs(0,0))*dthginv
  drsdphi = 0.5d0  &
  &       *(fnc_rs(1,1) - fnc_rs(1,0) &
  &       + fnc_rs(0,1) - fnc_rs(0,0))*dphiginv
  hrs = 0.25d0 &
  &   *(fnc_rs(1,1) + fnc_rs(0,1) &
  &   + fnc_rs(1,0) + fnc_rs(0,0))
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
end subroutine flgrad_type0
