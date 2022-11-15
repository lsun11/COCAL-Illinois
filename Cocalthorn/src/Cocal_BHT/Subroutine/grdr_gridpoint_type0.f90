subroutine grdr_gridpoint_type0(fnc,deriv,irg,itg,ipg)
  use phys_constant, only : long
  use grid_parameter, only : nrg
  use coordinate_grav_extended, only : rgex, irgex_r, itgex_r, ipgex_r
  implicit none
  real(long), external :: dfdx_4th
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: deriv
  real(long) :: rv
  real(long) :: r5(5), fr5(5)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, irg0, ii
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ir0 = min0(irg-2,nrg-4)
!
  do ii = 1, 5
    irg0 = ir0 + ii - 1
    r5(ii) = rgex(irg0)
!
    irgex = irgex_r(irg0)
    itgex = itgex_r(itg,irg0)
    ipgex = ipgex_r(ipg,irg0)
    fr5(ii) = fnc(irgex,itgex,ipgex)
  end do
!
! --- To cartesian component.
!
  rv = rgex(irg)
  deriv = dfdx_4th(r5,fr5,rv)
!
end subroutine grdr_gridpoint_type0
