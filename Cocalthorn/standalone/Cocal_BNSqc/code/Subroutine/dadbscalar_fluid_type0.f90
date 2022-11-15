subroutine dadbscalar_type0(fnc,d2f,irf,itf,ipf)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, ntfxy, npfyzp
  use coordinate_grav_r, only : rginv
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  use coordinate_grav_extended
  use trigonometry_grav_theta, only : sinthg, costhg, cosecthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_matter, only : rs
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long) :: d2f(1:3,1:3)
  real(long) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv
  real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  real(long) ::  rs_t4(4), rs_p4(4)
  real(long) :: dfncdx(0:1,0:1,0:1), dfncdy(0:1,0:1,0:1), dfncdz(0:1,0:1,0:1)
  real(long) :: fnc_rs(0:1,0:1)
  real(long) :: d2fdxdx, d2fdxdy, d2fdxdz, &
              & d2fdydx, d2fdydy, d2fdydz, &
              & d2fdzdx, d2fdzdy, d2fdzdz
  integer :: irf, itf, ipf, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irf0 , itf0 , ipf0, ii
  integer :: ip, ir, it, irf8, itf8, ipf8
  real(long), external :: dfdx_3rd
!
! --- Compute DaDb fnc of a function in 3rd+2nd order.
! --- The DaDb fnc is evaluated on the mid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irf-2,nrf-3)
  it0 = itf-2
  ip0 = ipf-2
!
  do ip = 0, 1
    ipf8 = ipf - 1 + ip
    do it = 0, 1
      itf8 = itf - 1 + it
      do ir = 0, 1
        irf8 = irf - 1 + ir
!
        do ii = 1, 4
          irf0 = ir0 + ii - 1
          itf0 = it0 + ii - 1
          ipf0 = ip0 + ii - 1
!
          r4(ii) = rgex(irf0)
          th4(ii) = thgex(itf0)
          phi4(ii) = phigex(ipf0)
!
          if (irf8.eq.0) then      
            irgex = irgex_r(irf0)
            itgex = ntfxy
            ipgex = ipgex_r(0,irf0)
            r4(ii)  = rs(itgex,ipgex)*rgex(irf0)
            fr4(ii) = fnc(irgex,itgex,ipgex)
            irgex = irgex_r(irf0)
            itgex = ntfxy
            ipgex = ipgex_r(npfyzp,irf0)
            th4(ii) = rs(itgex,ipgex)*rgex(irf0)
            ft4(ii) = fnc(irgex,itgex,ipgex)
            irgex = irgex_r(irf0)
            itgex = itgex_r(0,irf0)
            ipgex = 0
            phi4(ii)= rs(itgex,ipgex)*rgex(irf0)
            fp4(ii) = fnc(irgex,itgex,ipgex)
          else
            if (itf8.eq.0.or.itf8.eq.ntf) then
              irgex = irf8
              itgex = itgex_th(itf0)
              ipgex = ipgex_th(0,itf0)
              fr4(ii) = fnc(irgex,itgex,ipgex)
              rs_t4(ii) = rs(itgex,ipgex)
              irgex = irf8
              itgex = itgex_th(itf0)
              ipgex = ipgex_th(npfyzp,itf0)
              ft4(ii) = fnc(irgex,itgex,ipgex)
              rs_p4(ii) = rs(itgex,ipgex)
              irgex = irgex_r(irf0)
              itgex = itgex_r(itf,irf0)
              ipgex = 0
              r4(ii)= rs(itgex,ipgex)*rgex(irf0)
              fp4(ii) = fnc(irgex,itgex,ipgex)
            else
              irgex = irgex_r(irf0)
              itgex = itgex_r(itf8,irf0)
              ipgex = ipgex_r(ipf8,irf0)
              r4(ii)= rs(itgex,ipgex)*rgex(irf0)
              fr4(ii) = fnc(irgex,itgex,ipgex)
              irgex = irf8
              itgex = itgex_th(itf0)
              ipgex = ipgex_th(ipf8,itf0)
              ft4(ii) = fnc(irgex,itgex,ipgex)
              rs_t4(ii) = rs(itgex,ipgex)
              irgex = irf8
              itgex = itf8
              ipgex = ipgex_phi(ipf0)
              fp4(ii) = fnc(irgex,itgex,ipgex)
              rs_p4(ii) = rs(itgex,ipgex)
            end if
          end if
        end do
!
! --- To cartesian component.
!
        rv = rg(irf8)*rs(itf8,ipf8)
        tv = thg(itf8)
        pv = phig(ipf8)
        rsinv = 1.0d0/rs(itf8,ipf8)
!
        if (irf8.eq.0) then
          dfdx = dfdx_3rd(r4,fr4,rv)
          dfdy = dfdx_3rd(r4,ft4,rv)
          dfdz = dfdx_3rd(r4,fp4,rv)
        else
          if (itf8.eq.0.or.itf8.eq.ntf) then
            dfdrR= dfdx_3rd(r4,fp4,rv)
            dfdx = dfdx_3rd(th4,fr4,tv)*rginv(irf8)*rsinv &
            &    - dfdx_3rd(th4,rs_t4,tv)*rsinv*dfdrR
            dfdy = dfdx_3rd(th4,ft4,tv)*rginv(irf8)*rsinv &
            &    - dfdx_3rd(th4,rs_p4,tv)*rsinv*dfdrR
            dfdz = dfdrR
            if (itf8.eq.ntf) then
              dfdx = - dfdx
              dfdy = - dfdy
              dfdz = - dfdz
            end if
          else
            gr1  = dfdx_3rd(r4,fr4,rv)
            gr2  = dfdx_3rd(th4,ft4,tv)*rginv(irf8)*rsinv &
            &    - dfdx_3rd(th4,rs_t4,tv)*rsinv*gr1
            gr3  = dfdx_3rd(phi4,fp4,pv)*rginv(irf8)*cosecthg(itf8) &
            &    - dfdx_3rd(phi4,rs_p4,pv)*rsinv*cosecthg(itf8)*gr1
            dfdx = gr1*sinthg(itf8)*cosphig(ipf8) &
            &    + gr2*costhg(itf8)*cosphig(ipf8) &
            &    - gr3*sinphig(ipf8)
            dfdy = gr1*sinthg(itf8)*sinphig(ipf8) &
            &    + gr2*costhg(itf8)*sinphig(ipf8) &
            &    + gr3*cosphig(ipf8)
            dfdz = gr1*costhg(itf8)  &
            &    - gr2*sinthg(itf8)
          end if
        end if
!
        dfncdx(ir,it,ip) = dfdx
        dfncdy(ir,it,ip) = dfdy
        dfncdz(ir,it,ip) = dfdz
      end do
      fnc_rs(it,ip) = rs(itf8,ipf8)
    end do
  end do
!
  call flgrad_type0(dfncdx,fnc_rs,d2fdxdx,d2fdxdy,d2fdxdz,irf,itf,ipf)
  call flgrad_type0(dfncdy,fnc_rs,d2fdydx,d2fdydy,d2fdydz,irf,itf,ipf)
  call flgrad_type0(dfncdz,fnc_rs,d2fdzdx,d2fdzdy,d2fdzdz,irf,itf,ipf)
  d2f(1,1) = d2fdxdx
  d2f(1,2) = d2fdxdy
  d2f(1,3) = d2fdxdz
  d2f(2,1) = d2fdydx
  d2f(2,2) = d2fdydy
  d2f(2,3) = d2fdydz
  d2f(3,1) = d2fdzdx
  d2f(3,2) = d2fdzdy
  d2f(3,3) = d2fdzdz
!
end subroutine dadbscalar_fluid_type0
