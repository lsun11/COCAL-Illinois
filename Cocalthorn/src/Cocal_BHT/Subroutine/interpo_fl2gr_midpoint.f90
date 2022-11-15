subroutine interpo_fl2gr_midpoint(flv,grv)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_theta, only : hthg
  use coordinate_grav_phi, only : hphig
  use coordinate_grav_extended
  use def_matter, only : rs
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), external :: lagint_4th, lagint_2nd
  real(long), pointer :: grv(:,:,:), flv(:,:,:)
  real(long) :: r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  real(long) :: r2(2), fr2(2)
  real(long) :: rrff, rrgg, small = 1.0d-14
  real(long) :: rc, thc, phic
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii, jj, kk
  integer :: irf, irf0
!
  grv(0:nrg,0:ntg,0:npg) = 0.0d0
!
  do ipg = 1, npg
    ip0 = ipg-2
    do itg = 1, ntg
      it0 = itg-2
      call interpo_linear_type0_2Dsurf(rrff,rs,itg,ipg)
      do irg = 1, nrg
        rc   = hrg(irg)/rrff
        thc  = hthg(itg)
        phic = hphig(ipg)
!
        if (hrg(irg).gt.rrff*rg(nrf)) exit    ! rrff*rg is correct
!
        do irf = 0, nrf
          if (hrg(irg).le.rrff*rg(irf)) then
            irf0= irf-1
            ir0 = min0(irf-2,nrf-3)
            exit
          end if
        end do
!
        if (irf0.eq.nrf-1) then
!
          r2(1:2) = rgex(irf0:irf0+1)
          th4(1:4) = thgex(it0:it0+3)
          phi4(1:4) = phigex(ip0:ip0+3)
!
          do kk = 1, 4
            ipg0 = ip0 + kk - 1
            do jj = 1, 4
              itg0 = it0 + jj - 1
              do ii = 1, 2
                irg0 = irf0 + ii - 1
                irgex = irgex_r(irg0)
                itgex = itgex_r(itgex_th(itg0),irg0)
                ipgex = ipgex_r(ipgex_th(ipgex_phi(ipg0),itg0),irg0)
                fr2(ii) = flv(irgex,itgex,ipgex)
              end do
              ft4(jj) = lagint_2nd(r2,fr2,rc)
            end do
            fp4(kk) = lagint_4th(th4,ft4,thc)
          end do
          grv(irg,itg,ipg) = lagint_4th(phi4,fp4,phic)
!
        else
!
          r4(1:4) = rgex(ir0:ir0+3)
          th4(1:4) = thgex(it0:it0+3)
          phi4(1:4) = phigex(ip0:ip0+3)
!
          do kk = 1, 4
            ipg0 = ip0 + kk - 1
            do jj = 1, 4
              itg0 = it0 + jj - 1
              do ii = 1, 4
                irg0 = ir0 + ii - 1
                irgex = irgex_r(irg0)
                itgex = itgex_r(itgex_th(itg0),irg0)
                ipgex = ipgex_r(ipgex_th(ipgex_phi(ipg0),itg0),irg0)
                fr4(ii) = flv(irgex,itgex,ipgex)
              end do
              ft4(jj) = lagint_4th(r4,fr4,rc)
            end do
            fp4(kk) = lagint_4th(th4,ft4,thc)
          end do
          grv(irg,itg,ipg) = lagint_4th(phi4,fp4,phic)
        end if
!             
      end do
    end do
  end do
!
end subroutine interpo_fl2gr_midpoint
