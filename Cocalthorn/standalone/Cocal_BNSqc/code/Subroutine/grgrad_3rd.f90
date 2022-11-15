subroutine grgrad_3rd(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, ntgxy, npgyzp
  use coordinate_grav_r, only : rginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use coordinate_grav_extended
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv
  real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii
  real(long), external :: dfdx_3rd
!
! --- Compute the gradient of a function in 3rd order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irg-2,nrg-3)
  it0 = itg-2
  ip0 = ipg-2
!
  do ii = 1, 4
    irg0 = ir0 + ii - 1
    itg0 = it0 + ii - 1
    ipg0 = ip0 + ii - 1
!
    r4(ii) = rgex(irg0)
    th4(ii) = thgex(itg0)
    phi4(ii) = phigex(ipg0)
!
    if (irg.eq.0) then
      irgex = irgex_r(irg0)
      itgex = ntgxy
      ipgex = ipgex_r(0,irg0)
      fr4(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irg0)
      itgex = ntgxy
      ipgex = ipgex_r(npgyzp,irg0)
      ft4(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irg0)
      itgex = itgex_r(0,irg0)
      ipgex = 0
      fp4(ii) = fnc(irgex,itgex,ipgex)
    else
      if (itg.eq.0.or.itg.eq.ntg) then
        irgex = irg
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(0,itg0)
        fr4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(npgyzp,itg0)
        ft4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irgex_r(irg0)
        itgex = itgex_r(itg,irg0)
        ipgex = 0
        fp4(ii) = fnc(irgex,itgex,ipgex)
      else
        irgex = irgex_r(irg0)
        itgex = itgex_r(itg,irg0)
        ipgex = ipgex_r(ipg,irg0)
        fr4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(ipg,itg0)
        ft4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg
        itgex = itg
        ipgex = ipgex_phi(ipg0)
        fp4(ii) = fnc(irgex,itgex,ipgex)
      end if
    end if
  end do
!
! --- To cartesian component.
!
  rv = rg(irg)
  tv = thg(itg)
  pv = phig(ipg)
!
  if (irg.eq.0) then
    dfdx = dfdx_3rd(r4,fr4,rv)
    dfdy = dfdx_3rd(r4,ft4,rv)
    dfdz = dfdx_3rd(r4,fp4,rv)
  else
    if (itg.eq.0) then
      dfdx = dfdx_3rd(th4,fr4,tv)*rginv(irg)
      dfdy = dfdx_3rd(th4,ft4,tv)*rginv(irg)
      dfdz = dfdx_3rd(r4,fp4,rv)
    else if (itg.eq.ntg) then
      dfdx = - dfdx_3rd(th4,fr4,tv)*rginv(irg)
      dfdy = - dfdx_3rd(th4,ft4,tv)*rginv(irg)
      dfdz = - dfdx_3rd(r4,fp4,rv)
    else
      gr1  = dfdx_3rd(r4,fr4,rv)
      gr2  = dfdx_3rd(th4,ft4,tv)*rginv(irg)
      gr3  = dfdx_3rd(phi4,fp4,pv)*rginv(irg)*cosecthg(itg)
      dfdx = gr1*sinthg(itg)*cosphig(ipg) &
    &      + gr2*costhg(itg)*cosphig(ipg) &
    &      - gr3*sinphig(ipg)
      dfdy = gr1*sinthg(itg)*sinphig(ipg) &
    &      + gr2*costhg(itg)*sinphig(ipg) &
    &      + gr3*cosphig(ipg)
      dfdz = gr1*costhg(itg)  &
    &      - gr2*sinthg(itg)
    end if
  end if
end subroutine grgrad_3rd
