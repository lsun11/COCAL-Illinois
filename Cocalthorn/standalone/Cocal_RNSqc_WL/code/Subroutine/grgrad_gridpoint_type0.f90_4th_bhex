subroutine grgrad_gridpoint_type0(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg, ntgeq, npgyzp
  use coordinate_grav_r, only : rginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use coordinate_grav_extended
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv
  real(long) ::  r5(5), th5(5), phi5(5), fr5(5), ft5(5), fp5(5)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii
  real(long), external :: dfdx_4th
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(max0(0,irg-2),nrg-4)
  it0 = itg-2
  ip0 = ipg-2
!
  do ii = 1, 5
    irg0 = ir0 + ii - 1
    itg0 = it0 + ii - 1
    ipg0 = ip0 + ii - 1
!
    r5(ii) = rgex(irg0)
    th5(ii) = thgex(itg0)
    phi5(ii) = phigex(ipg0)
!
    if (itg.eq.0.or.itg.eq.ntg) then
      irgex = irg
      itgex = itgex_th(itg0)
      ipgex = ipgex_th(0,itg0)
      fr5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irg
      itgex = itgex_th(itg0)
      ipgex = ipgex_th(npgyzp,itg0)
      ft5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irg0)
      itgex = itgex_r(itg,irg0)
      ipgex = 0
      fp5(ii) = fnc(irgex,itgex,ipgex)
    else
      irgex = irgex_r(irg0)
      itgex = itgex_r(itg,irg0)
      ipgex = ipgex_r(ipg,irg0)
      fr5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irg
      itgex = itgex_th(itg0)
      ipgex = ipgex_th(ipg,itg0)
      ft5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irg
      itgex = itg
      ipgex = ipgex_phi(ipg0)
      fp5(ii) = fnc(irgex,itgex,ipgex)
    end if
  end do
!
! --- To cartesian component.
!
  rv = rg(irg)
  tv = thg(itg)
  pv = phig(ipg)
!
  if (itg.eq.0) then
    dfdx = dfdx_4th(th5,fr5,tv)*rginv(irg)
    dfdy = dfdx_4th(th5,ft5,tv)*rginv(irg)
    dfdz = dfdx_4th(r5,fp5,rv)
  else if (itg.eq.ntg) then
    dfdx = - dfdx_4th(th5,fr5,tv)*rginv(irg)
    dfdy = - dfdx_4th(th5,ft5,tv)*rginv(irg)
    dfdz = - dfdx_4th(r5,fp5,rv)
  else
    gr1  = dfdx_4th(r5,fr5,rv)
    gr2  = dfdx_4th(th5,ft5,tv)*rginv(irg)
    gr3  = dfdx_4th(phi5,fp5,pv)*rginv(irg)*cosecthg(itg)
    dfdx = gr1*sinthg(itg)*cosphig(ipg) &
  &      + gr2*costhg(itg)*cosphig(ipg) &
  &      - gr3*sinphig(ipg)
    dfdy = gr1*sinthg(itg)*sinphig(ipg) &
  &      + gr2*costhg(itg)*sinphig(ipg) &
  &      + gr3*cosphig(ipg)
    dfdz = gr1*costhg(itg)  &
  &      - gr2*sinthg(itg)
  end if
end subroutine grgrad_gridpoint_type0
