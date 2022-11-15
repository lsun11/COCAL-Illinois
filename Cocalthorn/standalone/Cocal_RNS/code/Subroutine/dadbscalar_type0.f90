subroutine dadbscalar_type0(fnc,d2f,irg,itg,ipg)
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
  real(long) :: d2f(1:3,1:3)
  real(long) :: dfdx, dfdy, dfdz
  real(long) ::  gr1, gr2, gr3, rv, tv, pv
  real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  real(long) :: dfncdx(0:1,0:1,0:1), dfncdy(0:1,0:1,0:1), dfncdz(0:1,0:1,0:1)
  real(long) :: d2fdxdx, d2fdxdy, d2fdxdz, &
              & d2fdydx, d2fdydy, d2fdydz, &
              & d2fdzdx, d2fdzdy, d2fdzdz
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii
  integer :: ip, ir, it, irg8, itg8, ipg8
  real(long), external :: dfdx_3rd
!
! --- Compute the DaDb fnc of a function in 3rd+2nd order.
! --- The DaDb fnc is evaluated on the mid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irg-2,nrg-3)
  it0 = itg-2
  ip0 = ipg-2
!
  do ip = 0, 1
    ipg8 = ipg - 1 + ip
    do it = 0, 1
      itg8 = itg - 1 + it
      do ir = 0, 1
        irg8 = irg - 1 + ir
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
    if (irg8.eq.0) then
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
      if (itg8.eq.0.or.itg8.eq.ntg) then
        irgex = irg8
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(0,itg0)
        fr4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg8
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(npgyzp,itg0)
        ft4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irgex_r(irg0)
        itgex = itgex_r(itg,irg0)
        ipgex = 0
        fp4(ii) = fnc(irgex,itgex,ipgex)
      else
        irgex = irgex_r(irg0)
        itgex = itgex_r(itg8,irg0)
        ipgex = ipgex_r(ipg8,irg0)
        fr4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg8
        itgex = itgex_th(itg0)
        ipgex = ipgex_th(ipg8,itg0)
        ft4(ii) = fnc(irgex,itgex,ipgex)
        irgex = irg8
        itgex = itg8
        ipgex = ipgex_phi(ipg0)
        fp4(ii) = fnc(irgex,itgex,ipgex)
      end if
    end if
  end do
!
! --- To cartesian component.
!
        rv = rg(irg8)
        tv = thg(itg8)
        pv = phig(ipg8)
!
        if (irg8.eq.0) then
          dfdx = dfdx_3rd(r4,fr4,rv)
          dfdy = dfdx_3rd(r4,ft4,rv)
          dfdz = dfdx_3rd(r4,fp4,rv)
        else
          if (itg8.eq.0) then
            dfdx = dfdx_3rd(th4,fr4,tv)*rginv(irg8)
            dfdy = dfdx_3rd(th4,ft4,tv)*rginv(irg8)
            dfdz = dfdx_3rd(r4,fp4,rv)
          else if (itg8.eq.ntg) then
            dfdx = - dfdx_3rd(th4,fr4,tv)*rginv(irg8)
            dfdy = - dfdx_3rd(th4,ft4,tv)*rginv(irg8)
            dfdz = - dfdx_3rd(r4,fp4,rv)
          else
            gr1  = dfdx_3rd(r4,fr4,rv)
            gr2  = dfdx_3rd(th4,ft4,tv)*rginv(irg8)
            gr3  = dfdx_3rd(phi4,fp4,pv)*rginv(irg8)*cosecthg(itg8)
            dfdx = gr1*sinthg(itg8)*cosphig(ipg8) &
          &      + gr2*costhg(itg8)*cosphig(ipg8) &
          &      - gr3*sinphig(ipg8)
            dfdy = gr1*sinthg(itg8)*sinphig(ipg8) &
          &      + gr2*costhg(itg8)*sinphig(ipg8) &
          &      + gr3*cosphig(ipg8)
            dfdz = gr1*costhg(itg8)  &
          &      - gr2*sinthg(itg8)
          end if
        end if
!
        dfncdx(ir,it,ip) = dfdx
        dfncdy(ir,it,ip) = dfdy
        dfncdz(ir,it,ip) = dfdz
      end do
    end do
  end do
!
  call grgrad_type0(dfncdx,d2fdxdx,d2fdxdy,d2fdxdz,irg,itg,ipg)
  call grgrad_type0(dfncdy,d2fdydx,d2fdydy,d2fdydz,irg,itg,ipg)
  call grgrad_type0(dfncdz,d2fdzdx,d2fdzdy,d2fdzdz,irg,itg,ipg)
  d2f(1,1) = d2fdxdx
  d2f(1,2) =(d2fdxdy + d2fdydx)*0.5d0
  d2f(1,3) =(d2fdxdz + d2fdzdx)*0.5d0
  d2f(2,1) =(d2fdxdy + d2fdydx)*0.5d0
  d2f(2,2) = d2fdydy
  d2f(2,3) =(d2fdydz + d2fdzdy)*0.5d0
  d2f(3,1) =(d2fdxdz + d2fdzdx)*0.5d0
  d2f(3,2) =(d2fdydz + d2fdzdy)*0.5d0
  d2f(3,3) = d2fdzdz
!
end subroutine dadbscalar_type0
