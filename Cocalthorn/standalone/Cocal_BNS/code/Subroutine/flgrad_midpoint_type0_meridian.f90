subroutine flgrad_midpoint_type0_meridian(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
  use phys_constant, only : long
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer :: fnc(:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) :: gr1, gr2, gr3, hrs, hrsinv
  real(long) :: dfncdr, dfncdth, dfncdphi, drsdth, drsdphi
  integer, intent(in) :: irf, itf, ipf
!
! --- Compute the gradient of a function at fixed phi gridpoint.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, derivatives. Phi is fixed
!
  dfncdr   = 0.5d0 &
  &     *(fnc(irf,itf  ) - fnc(irf-1,itf  ) &
  &     + fnc(irf,itf-1) - fnc(irf-1,itf-1))*drginv(irf)
  dfncdth  = 0.5d0 &
  &     *(fnc(irf  ,itf) - fnc(irf  ,itf-1) &
  &     + fnc(irf-1,itf) - fnc(irf-1,itf-1))*dthginv
 
  drsdth = (rs(itf,ipf) - rs(itf-1,ipf))*dthginv

  hrs    = 0.5d0*(rs(itf,ipf) + rs(itf-1,ipf))

  hrsinv = 1.0d0/hrs
!
! --- To cartesian component.
!
  gr1  = dfncdr*hrsinv
  gr2  = dfncdth*hrginv(irf)*hrsinv &
  &    - drsdth*hrsinv*gr1
  gr3  = 0.0d0
!
  dfdx = gr1 * hsinthg(itf) * cosphig(ipf) &
  &    + gr2 * hcosthg(itf) * cosphig(ipf) &
  &    - gr3 * sinphig(ipf)
  dfdy = gr1 * hsinthg(itf) * sinphig(ipf) &
  &    + gr2 * hcosthg(itf) * sinphig(ipf) &
  &    + gr3 * cosphig(ipf)
  dfdz = gr1 * hcosthg(itf)  &
  &    - gr2 * hsinthg(itf)
!
end subroutine flgrad_midpoint_type0_meridian
