subroutine interpo_gr2fl_midpoint_type0(flv,grv,irf,itf,ipf)
  use phys_constant, only : long
  use grid_parameter, only : nrg
  use coordinate_grav_r, only : hrg
  use coordinate_grav_extended, only : hrgex, irgex_hr, itgex_hr, ipgex_hr
  use def_matter, only : rs
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer :: grv(:,:,:)
  real(long), intent(out) :: flv
  real(long) :: x(4), f(4)
  real(long) :: rrff, rrgg, small = 1.0d-14
  integer :: irf, itf, ipf, irg, ir0, irf0, irg0, itg0, ipg0
!
! Interpolete values on midpoints of spherical coordinates 
!          to values on midpoints of surface fitted coordiantes
!
  call interpo_linear_type0_2Dsurf(rrff,rs,itf,ipf)
  do irg = 1, nrg
    if (rrff*hrg(irf).le.hrg(irg)) then
      ir0 = min0(irg-2,nrg-3)
      exit
    end if
  end do
  do ii = 1, 4
    irf0 = ir0 + ii - 1
    irg0 = irgex_hr(irf0)
    itg0 = itgex_hr(itg,irf0)
    ipg0 = ipgex_hr(ipg,irf0)
    x(ii) = hrgex(irf0)
    f(ii) = grv(irg0,itg0,ipg0)
  end do
  rrgg = rrff*hrg(irf)
  flv = lagint_4th(x,f,rrgg)
!
end subroutine interpo_gr2fl_midpoint_type0
