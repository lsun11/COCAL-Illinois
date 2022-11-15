subroutine interpo_fl2gr(flv,grv)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use coordinate_grav_extended
  use def_matter, only : rs
  implicit none
  real(long), external :: lagint_4th, lagint_2nd
  real(long), pointer :: grv(:,:,:), flv(:,:,:)
  real(long) :: x(4), f(4), x2(2), f2(2)
  real(long) :: rrff, rrgg, small = 1.0d-14
  integer :: irf, irg, itg, ipg, ir0, irf0, ii
  integer :: irgex, itgex, ipgex
!
  grv(0:nrg,0:ntg,0:npg) = 0.0d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        rrff = rs(itg,ipg)
        if (rg(irg).gt.rrff*rg(nrf)) exit
        do irf = 0, nrf
          if (rg(irg).le.rrff*rg(irf)) then 
            irf0= irf-1
            ir0 = min0(irf-2,nrf-3)
            exit
          end if
        end do
        if (irf0.eq.nrf-1) then
          x2(1:2) = rrff*rg(irf0:irf0+1)
          f2(1:2) = flv(irf0:irf0+1,itg,ipg)
          rrgg = rg(irg)
          grv(irg,itg,ipg) = lagint_2nd(x2,f2,rrgg)
        else 
          do ii = 1, 4
            irf0 = ir0 + ii - 1
            irgex = irgex_r(irf0)
            itgex = itgex_r(itg,irf0)
            ipgex = ipgex_r(ipg,irf0)
            x(ii) =  rs(itgex,ipgex)*rgex(irf0)
            f(ii) = flv(irgex,itgex,ipgex)
            rrgg = rg(irg)
            grv(irg,itg,ipg) = lagint_4th(x,f,rrgg)
          end do
        end if
      end do
    end do
  end do
!
end subroutine interpo_fl2gr
