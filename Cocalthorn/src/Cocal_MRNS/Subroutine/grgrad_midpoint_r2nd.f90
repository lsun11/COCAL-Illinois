subroutine grgrad_midpoint_r2nd(fnc,dfdx,dfdy,dfdz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : hrginv, drginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use trigonometry_grav_theta, only : hsinthg, hcosthg, hcosecthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(long) ::  gr1, gr2, gr3, dfncdr, dfncdth, dfncdphi
  integer :: irg, itg, ipg
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives.
!
  do irg = 1, nrg
    do itg = 1, ntg
      do ipg = 1, npg
        dfncdr = 0.25d0 &
      &     *(fnc(irg,itg  ,ipg  ) - fnc(irg-1,itg  ,ipg  ) &
      &     + fnc(irg,itg-1,ipg  ) - fnc(irg-1,itg-1,ipg  ) &
      &     + fnc(irg,itg  ,ipg-1) - fnc(irg-1,itg  ,ipg-1) &
      &     + fnc(irg,itg-1,ipg-1) - fnc(irg-1,itg-1,ipg-1))*drginv(irg)
        dfncdth = 0.25d0 &
      &     *(fnc(irg  ,itg,ipg  ) - fnc(irg  ,itg-1,ipg  ) &
      &     + fnc(irg-1,itg,ipg  ) - fnc(irg-1,itg-1,ipg  ) &
      &     + fnc(irg  ,itg,ipg-1) - fnc(irg  ,itg-1,ipg-1) &
      &     + fnc(irg-1,itg,ipg-1) - fnc(irg-1,itg-1,ipg-1))*dthginv
        dfncdphi = 0.25d0  &
      &     *(fnc(irg  ,itg  ,ipg) - fnc(irg  ,itg  ,ipg-1) &
      &     + fnc(irg-1,itg  ,ipg) - fnc(irg-1,itg  ,ipg-1) &
      &     + fnc(irg  ,itg-1,ipg) - fnc(irg  ,itg-1,ipg-1) &
      &     + fnc(irg-1,itg-1,ipg) - fnc(irg-1,itg-1,ipg-1))*dphiginv
!
! --- To cartesian component.
!
        gr1  = dfncdr
        gr2  = dfncdth*hrginv(irg)
        gr3  = dfncdphi*hrginv(irg)*hcosecthg(itg)
        dfdx(irg,itg,ipg) = gr1 * hsinthg(itg) * hcosphig(ipg) &
      &                   + gr2 * hcosthg(itg) * hcosphig(ipg) &
      &                   - gr3 * hsinphig(ipg)
        dfdy(irg,itg,ipg) = gr1 * hsinthg(itg) * hsinphig(ipg) &
      &                   + gr2 * hcosthg(itg) * hsinphig(ipg) &
      &                   + gr3 * hcosphig(ipg)
        dfdz(irg,itg,ipg) = gr1 * hcosthg(itg)  &
      &                   - gr2 * hsinthg(itg)
      end do
    end do
  end do
end subroutine grgrad_midpoint_r2nd
