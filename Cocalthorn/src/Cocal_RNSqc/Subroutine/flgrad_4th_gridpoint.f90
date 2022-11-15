subroutine flgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, ntfeq, npfyzp
  use coordinate_grav_r, only : rginv
  use coordinate_grav_extended
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv, rsinv, dfdrR
  real(long) ::  r5(5), th5(5), phi5(5), fr5(5), ft5(5), fp5(5)
  real(long) ::  rs_t5(5), rs_p5(5)
  integer :: irf, itf, ipf, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irf0 , itf0 , ipf0, ii
  real(long), external :: dfdx_4th
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irf-2,nrf-4)
  it0 = itf-2
  ip0 = ipf-2
  rsinv = 1.0d0/rs(itf,ipf)
!
  do ii = 1, 5
    irf0 = ir0 + ii - 1
    itf0 = it0 + ii - 1
    ipf0 = ip0 + ii - 1
!
    th5(ii) = thgex(itf0)
    phi5(ii) = phigex(ipf0)
!
    if (irf.eq.0) then
      irgex = irgex_r(irf0)
      itgex = ntfeq
      ipgex = ipgex_r(0,irf0)
      r5(ii)  = rs(itgex,ipgex)*rgex(irf0)
      fr5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irf0)
      itgex = ntfeq
      ipgex = ipgex_r(npfyzp,irf0)
      th5(ii) = rs(itgex,ipgex)*rgex(irf0)
      ft5(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irf0)
      itgex = itgex_r(0,irf0)
      ipgex = 0
      phi5(ii)= rs(itgex,ipgex)*rgex(irf0)
      fp5(ii) = fnc(irgex,itgex,ipgex)
    else
      if (itf.eq.0.or.itf.eq.ntf) then
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(0,itf0)
        fr5(ii) = fnc(irgex,itgex,ipgex)
        rs_t5(ii) = rs(itgex,ipgex)
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(npfyzp,itf0)
        ft5(ii) = fnc(irgex,itgex,ipgex)
        rs_p5(ii) = rs(itgex,ipgex)
        irgex = irgex_r(irf0)
        itgex = itgex_r(itf,irf0)
        ipgex = 0
        r5(ii)  = rs(itgex,ipgex)*rgex(irf0)
        fp5(ii) = fnc(irgex,itgex,ipgex)
      else
        irgex = irgex_r(irf0)
        itgex = itgex_r(itf,irf0)
        ipgex = ipgex_r(ipf,irf0)
        r5(ii)  = rs(itgex,ipgex)*rgex(irf0)
        fr5(ii) = fnc(irgex,itgex,ipgex)
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(ipf,itf0)
        ft5(ii) = fnc(irgex,itgex,ipgex)
        rs_t5(ii) = rs(itgex,ipgex)
        irgex = irf
        itgex = itf
        ipgex = ipgex_phi(ipf0)
        fp5(ii) = fnc(irgex,itgex,ipgex)
        rs_p5(ii) = rs(itgex,ipgex)
      end if
    end if
  end do
!
! --- To cartesian component.
!
  rv = rg(irf)*rs(itf,ipf)
  tv = thg(itf)
  pv = phig(ipf)
!
  if (irf.eq.0) then
    dfdx = dfdx_4th(r5,fr5,rv)
    dfdy = dfdx_4th(r5,ft5,rv)
    dfdz = dfdx_4th(r5,fp5,rv)
  else
    if (itf.eq.0.or.itf.eq.ntf) then
      dfdrR= dfdx_4th(r5,fp5,rv)
      dfdx = dfdx_4th(th5,fr5,tv)*rginv(irf)*rsinv &
      &    - dfdx_4th(th5,rs_t5,tv)*rsinv*dfdrR
      dfdy = dfdx_4th(th5,ft5,tv)*rginv(irf)*rsinv &
      &    - dfdx_4th(th5,rs_p5,tv)*rsinv*dfdrR
      dfdz = dfdrR
      if (itf.eq.ntf) then
        dfdx = - dfdx
        dfdy = - dfdy
        dfdz = - dfdz
      end if
    else
      gr1  = dfdx_4th(r5,fr5,rv)
      gr2  = dfdx_4th(th5,ft5,tv)*rginv(irf)*rsinv &
      &    - dfdx_4th(th5,rs_t5,tv)*rsinv*gr1
      gr3  = dfdx_4th(phi5,fp5,pv)*rginv(irf)*rsinv*cosecthg(itf) &
      &    - dfdx_4th(phi5,rs_p5,pv)*rsinv*cosecthg(itf)*gr1
      dfdx = gr1*sinthg(itf)*cosphig(ipf) &
    &      + gr2*costhg(itf)*cosphig(ipf) &
    &      - gr3*sinphig(ipf)
      dfdy = gr1*sinthg(itf)*sinphig(ipf) &
    &      + gr2*costhg(itf)*sinphig(ipf) &
    &      + gr3*cosphig(ipf)
      dfdz = gr1*costhg(itf)  &
    &      - gr2*sinthg(itf)
    end if
  end if
end subroutine flgrad_4th_gridpoint
