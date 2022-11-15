subroutine interpo_flsfc2flsph_midpoint(flv,flsphv)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_extended, only : irgex_hr, itgex_hr, ipgex_hr, hrgex
  use def_matter, only : rs
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), external :: lagint_4th
  real(long), external :: lagint_2nd
  real(long), pointer :: flsphv(:,:,:), flv(:,:,:)
  real(long) :: x(4), f(4), x2(2), f2(2)
  real(long) :: rrff, rrgg, small = 1.0d-14
  integer :: irf, irg, itg, ipg, ir0, irf0, irg0, itg0, ipg0, ii
!
! Interpolete values on midpoints of surface fitted coordiantes
!          to values on midpoints of spherical coordinates 
!
  flsphv(1:nrf,1:ntf,1:npf) = 0.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(rrff,rs,itg,ipg)
      do irg = 1, nrf
        if (hrg(irg).gt.rrff*rg(nrf)) exit    ! rrff*rg is correct
        if (hrg(irg).gt.rrff*hrg(nrf-1)) then  ! rrff*hrg is correct
          ir0 = nrf-1
          do ii = 1, 2
            irf0 = ir0 + ii - 1
            irg0 = irgex_hr(irf0)
            itg0 = itgex_hr(itg,irf0)
            ipg0 = ipgex_hr(ipg,irf0)
            call interpo_linear_type0_2Dsurf(rrff,rs,itg0,ipg0)
            x2(ii) = rrff*hrgex(irf0)
            f2(ii) = flv(irg0,itg0,ipg0)
          end do
          rrgg = hrg(irg)
          flsphv(irg,itg,ipg) = lagint_2nd(x2,f2,rrgg)
          exit
        end if
!
        do irf = 1, nrf
          if (hrg(irg).le.rrff*hrg(irf)) then ! rrff*hrg is correct
            ir0 = min0(irf-2,nrf-3)
            exit
          end if
        end do
        do ii = 1, 4
          irf0 = ir0 + ii - 1
          irg0 = irgex_hr(irf0)
          itg0 = itgex_hr(itg,irf0)
          ipg0 = ipgex_hr(ipg,irf0)
          call interpo_linear_type0_2Dsurf(rrff,rs,itg0,ipg0)
          x(ii) = rrff*hrgex(irf0)
          f(ii) = flv(irg0,itg0,ipg0)
        end do
        rrgg = hrg(irg)
        flsphv(irg,itg,ipg) = lagint_4th(x,f,rrgg)
      end do
    end do
  end do
!
end subroutine interpo_flsfc2flsph_midpoint
