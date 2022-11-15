subroutine grdr_gridpoint_type0_3rd_nosym(fnc,deriv,irg,itg,ipg)
  use phys_constant, only : long
  use grid_parameter, only : nrg
  use coordinate_grav_extended, only : rgex, irgex_r, itgex_r, ipgex_r
  implicit none
  real(long), external :: dfdx_3rd
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: deriv
  real(long) :: rv
  real(long) :: r4(4), fr4(4)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, irg0, ii
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = max0(0,min0(irg-1,nrg-3))
!
  do ii = 1, 4
    irg0 = ir0 + ii - 1
    r4(ii) = rgex(irg0)
!
    irgex = irgex_r(irg0)
    itgex = itgex_r(itg,irg0)
    ipgex = ipgex_r(ipg,irg0)
    fr4(ii) = fnc(irgex,itgex,ipgex)
  end do
!
! --- To cartesian component.
!
  rv = rgex(irg)
  deriv = dfdx_3rd(r4,fr4,rv)
!
end subroutine grdr_gridpoint_type0_3rd_nosym
