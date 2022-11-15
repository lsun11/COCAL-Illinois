subroutine interpo_gr2fl_type0(flv,grv,irf,itf,ipf)
  use phys_constant, only : long
  use grid_parameter, only : nrg
  use coordinate_grav_r, only : rg
  use coordinate_grav_extended
  use def_matter, only : rs
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer :: grv(:,:,:)
  real(long), intent(out) :: flv
  real(long) :: x(4), f(4)
  real(long) :: rrff,   small = 1.0d-14
  integer :: irg, irf, itf, ipf, ir0, irf0, irg0, itg0, ipg0, ii
!
  rrff = rs(itf,ipf)*rg(irf)
  do irg = 0, nrg-1
    if (rrff.le.rg(irg)) then 
      ir0 = min0(irg-2,nrg-3)
      exit
    end if
  end do
  do ii = 1, 4
    irf0 = ir0 + ii - 1
    irg0 = irgex_r(irf0)
    itg0 = itgex_r(itf,irf0)
    ipg0 = ipgex_r(ipf,irf0)
    x(ii) = rgex(irf0)
    f(ii) = grv(irg0,itg0,ipg0)
  end do
  flv = lagint_4th(x,f,rrff)
!
end subroutine interpo_gr2fl_type0
