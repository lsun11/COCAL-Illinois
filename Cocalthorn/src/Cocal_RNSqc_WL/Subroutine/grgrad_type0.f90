subroutine grgrad_type0(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  implicit none
  real(long) :: fnc(0:1,0:1,0:1)
  real(long) :: dfdx
  real(long) :: dfdy
  real(long) :: dfdz
  real(long) ::  gr1, gr2, gr3, dfncdr, dfncdth, dfncdphi
  integer, intent(in) :: irg, itg, ipg
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives.
!
  dfncdr = 0.25d0 &
   &     *(fnc(1,1,1) - fnc(0,1,1) &
   &     + fnc(1,0,1) - fnc(0,0,1) &
   &     + fnc(1,1,0) - fnc(0,1,0) &
   &     + fnc(1,0,0) - fnc(0,0,0))*drginv(irg)
  dfncdth = 0.25d0 &
    &     *(fnc(1,1,1) - fnc(1,0,1) &
    &     + fnc(0,1,1) - fnc(0,0,1) &
    &     + fnc(1,1,0) - fnc(1,0,0) &
    &     + fnc(0,1,0) - fnc(0,0,0))*dthginv
  dfncdphi = 0.25d0  &
     &     *(fnc(1,1,1) - fnc(1,1,0) &
     &     + fnc(0,1,1) - fnc(0,1,0) &
     &     + fnc(1,0,1) - fnc(1,0,0) &
     &     + fnc(0,0,1) - fnc(0,0,0))*dphiginv
!
! --- To cartesian component.
!
  gr1  = dfncdr
  gr2  = dfncdth*hrginv(irg)
  gr3  = dfncdphi*hrginv(irg)*hcosecthg(itg)
  dfdx = gr1 * hsinthg(itg) * hcosphig(ipg) &
&      + gr2 * hcosthg(itg) * hcosphig(ipg) &
&      - gr3 * hsinphig(ipg)
  dfdy = gr1 * hsinthg(itg) * hsinphig(ipg) &
&      + gr2 * hcosthg(itg) * hsinphig(ipg) &
&      + gr3 * hcosphig(ipg)
  dfdz = gr1 * hcosthg(itg)  &
&      - gr2 * hsinthg(itg)
!
end subroutine grgrad_type0
