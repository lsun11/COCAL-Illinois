subroutine interpo_flsph2flsfc(flsphv,flv)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use def_matter, only : rs
  use coordinate_grav_extended, only : irgex_r, itgex_r, ipgex_r, rgex
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer :: flsphv(:,:,:), flv(:,:,:)
  real(long) :: x(4), f(4)
  real(long) :: rrff,   small = 1.0d-14
  integer :: irg, irf, itf, ipf, ir0, irf0, irg0, itg0, ipg0, ii
!
  flv(0:nrf,0:ntf,0:npf) = 0.0d0
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rrff = rs(itf,ipf)*rg(irf)
        ir0 = nrf-3
        do irg = 1, nrf
          if (rrff.le.rg(irg)) then 
            ir0 = min0(irg-2,nrf-3)
            exit
          end if
        end do
        do ii = 1, 4
          irf0 = ir0 + ii - 1
          irg0 = irgex_r(irf0)
          itg0 = itgex_r(itf,irf0)
          ipg0 = ipgex_r(ipf,irf0)
          x(ii) = rgex(irf0)
          f(ii) = flsphv(irg0,itg0,ipg0)
        end do
        flv(irf,itf,ipf) = lagint_4th(x,f,rrff)
      end do
    end do
  end do
!
end subroutine interpo_flsph2flsfc
