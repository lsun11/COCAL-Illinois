subroutine flgrad_2nd_gridpoint_export(fnc,dfdx,dfdy,dfdz,irf,itf,ipf,rs)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, ntfeq, npfyzp
  use coordinate_grav_r, only : rginv
  use coordinate_grav_extended
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
!  use def_matter, only : rs
  implicit none
  real(long), pointer     :: fnc(:,:,:), rs(:,:)
  real(long), intent(out) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv, rsinv, dfdrR
  real(long) ::  r3(3), th3(3), phi3(3), fr3(3), ft3(3), fp3(3)
  real(long) ::  rs_t3(3), rs_p3(3)
  integer :: irf, itf, ipf, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irf0 , itf0 , ipf0, ii
  real(long), external :: dfdx_2nd
!
! --- Compute the gradient of a function in 2nd order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irf-1,nrf-2)
  it0 = itf-1
  ip0 = ipf-1
  rsinv = 1.0d0/rs(itf,ipf)
!
  do ii = 1, 3
    irf0 = ir0 + ii - 1
    itf0 = it0 + ii - 1
    ipf0 = ip0 + ii - 1
!
    th3(ii) = thgex(itf0)
    phi3(ii) = phigex(ipf0)
!
    if (irf.eq.0) then
      irgex = irgex_r(irf0)
      itgex = ntfeq
      ipgex = ipgex_r(0,irf0)
      r3(ii)  = rs(itgex,ipgex)*rgex(irf0)
      fr3(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irf0)
      itgex = ntfeq
      ipgex = ipgex_r(npfyzp,irf0)
      th3(ii) = rs(itgex,ipgex)*rgex(irf0)
      ft3(ii) = fnc(irgex,itgex,ipgex)
      irgex = irgex_r(irf0)
      itgex = itgex_r(0,irf0)
      ipgex = 0
      phi3(ii)= rs(itgex,ipgex)*rgex(irf0)
      fp3(ii) = fnc(irgex,itgex,ipgex)
    else
      if (itf.eq.0.or.itf.eq.ntf) then
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(0,itf0)
        fr3(ii) = fnc(irgex,itgex,ipgex)
        rs_t3(ii) = rs(itgex,ipgex)
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(npfyzp,itf0)
        ft3(ii) = fnc(irgex,itgex,ipgex)
        rs_p3(ii) = rs(itgex,ipgex)
        irgex = irgex_r(irf0)
        itgex = itgex_r(itf,irf0)
        ipgex = 0
        r3(ii)  = rs(itgex,ipgex)*rgex(irf0)
        fp3(ii) = fnc(irgex,itgex,ipgex)
      else
        irgex = irgex_r(irf0)
        itgex = itgex_r(itf,irf0)
        ipgex = ipgex_r(ipf,irf0)
        r3(ii)  = rs(itgex,ipgex)*rgex(irf0)
        fr3(ii) = fnc(irgex,itgex,ipgex)
        irgex = irf
        itgex = itgex_th(itf0)
        ipgex = ipgex_th(ipf,itf0)
        ft3(ii) = fnc(irgex,itgex,ipgex)
        rs_t3(ii) = rs(itgex,ipgex)
        irgex = irf
        itgex = itf
        ipgex = ipgex_phi(ipf0)
        fp3(ii) = fnc(irgex,itgex,ipgex)
        rs_p3(ii) = rs(itgex,ipgex)
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
    dfdx = dfdx_2nd(r3,fr3,rv)
    dfdy = dfdx_2nd(th3,ft3,rv)
    dfdz = dfdx_2nd(phi3,fp3,rv)
  else
    if (itf.eq.0.or.itf.eq.ntf) then
      dfdrR= dfdx_2nd(r3,fp3,rv)
      dfdx = dfdx_2nd(th3,fr3,tv)*rginv(irf)*rsinv &
      &    - dfdx_2nd(th3,rs_t3,tv)*rsinv*dfdrR  
      dfdy = dfdx_2nd(th3,ft3,tv)*rginv(irf)*rsinv &
      &    - dfdx_2nd(th3,rs_p3,tv)*rsinv*dfdrR  
      dfdz = dfdrR
      if (itf.eq.ntf) then
        dfdx = - dfdx
        dfdy = - dfdy
        dfdz = - dfdz
      end if
    else
      gr1  = dfdx_2nd(r3,fr3,rv)
      gr2  = dfdx_2nd(th3,ft3,tv)*rginv(irf)*rsinv &
      &    - dfdx_2nd(th3,rs_t3,tv)*rsinv*gr1
      gr3  = dfdx_2nd(phi3,fp3,pv)*rginv(irf)*rsinv*cosecthg(itf) &
      &    - dfdx_2nd(phi3,rs_p3,pv)*rsinv*cosecthg(itf)*gr1
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
end subroutine flgrad_2nd_gridpoint_export
