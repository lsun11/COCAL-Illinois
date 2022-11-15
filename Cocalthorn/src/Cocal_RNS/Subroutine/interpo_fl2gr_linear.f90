subroutine interpo_fl2gr_linear(flv,grv)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg
  use def_matter, only : rs
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer :: grv(:,:,:), flv(:,:,:)
  real(long) :: x(2), f(2)
  real(long) :: rrff, rrgg, small = 1.0d-14
  integer :: irf, irg, itg, ipg, ir0
!
  grv(0:nrg,0:ntg,0:npg) = 0.0d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      rrff = rs(itg,ipg)
      do irg = 0, nrg
        if (rg(irg).gt.rrff*rg(nrf)) exit
        do irf = 0, nrf
          if (rg(irg).le.rrff*rg(irf)) then 
            ir0 = min0(max0(0,irf-1),nrf-1)
            exit
          end if
        end do
        x(1:2) = rrff*rg(ir0:ir0+1)
        f(1:2) = flv(ir0:ir0+1,itg,ipg)
        rrgg = rg(irg)
        grv(irg,itg,ipg) = ((x(2)-rrgg)*f(1)+(rrgg-x(1))*f(2))/(x(2)-x(1))
      end do
    end do
  end do
!
end subroutine interpo_fl2gr_linear
