subroutine flgrad_midpoint_surface_type0(fnc,dfdx,dfdy,dfdz,itf,ipf)
  use phys_constant, only : long
  use grid_parameter, only : nrf
  use coordinate_grav_r, only : rg, rginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) :: gr1, gr2, gr3, rv, hrs, hrsinv
  real(long) :: dfncdr, dfncdth, dfncdphi, drsdth, drsdphi
  real(long) :: r3(3), fr3(3), dfncdr4(0:1,0:1)
  integer, intent(in) :: itf, ipf
  integer :: irf0, itf0, ipf0, ii, it0, ip0, irf
  real(long), external :: dfdx_2nd
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives.
!
  irf = nrf
!
  do ip0 = 0, 1
    do it0 = 0, 1
      itf0 = itf-1 + it0
      ipf0 = ipf-1 + ip0
      do ii = 1, 3
        irf0  = irf-3 + ii
        r3(ii) = rg(irf0)
        fr3(ii) = fnc(irf0,itf0,ipf0)
      end do
      rv = rg(irf)
      dfncdr4(it0,ip0) = dfdx_2nd(r3,fr3,rv)
    end do
  end do
!
  dfncdr = 0.25d0*(dfncdr4(0,0) + dfncdr4(0,1) &
  &              + dfncdr4(1,0) + dfncdr4(1,1))
  dfncdth  = 0.5d0 &
  &     *(fnc(irf  ,itf,ipf  ) - fnc(irf  ,itf-1,ipf  ) &
  &     + fnc(irf  ,itf,ipf-1) - fnc(irf  ,itf-1,ipf-1))*dthginv
  dfncdphi = 0.5d0  &
  &     *(fnc(irf  ,itf  ,ipf) - fnc(irf  ,itf  ,ipf-1) &
  &     + fnc(irf  ,itf-1,ipf) - fnc(irf  ,itf-1,ipf-1))*dphiginv
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
  gr2  = dfncdth*rginv(irf)*hrsinv &
  &    - drsdth*hrsinv*gr1
  gr3  = dfncdphi*rginv(irf)*hrsinv*hcosecthg(itf) &
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
end subroutine flgrad_midpoint_surface_type0
