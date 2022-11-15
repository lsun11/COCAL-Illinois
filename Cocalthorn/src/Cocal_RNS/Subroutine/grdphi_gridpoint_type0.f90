subroutine grdphi_gridpoint_type0(fnc,deriv,irg,itg,ipg)
  use phys_constant, only : long
  use coordinate_grav_extended
  implicit none
  real(long), external :: dfdx_4th
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: deriv
  real(long) :: pv
  real(long) :: phi5(5), fp5(5)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ip0, ipg0, ii
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ip0 = ipg-2
!
  do ii = 1, 5
    ipg0 = ip0 + ii - 1
    phi5(ii) = phigex(ipg0)
!
    irgex = irg
    itgex = itg
    ipgex = ipgex_phi(ipg0)
    fp5(ii) = fnc(irgex,itgex,ipgex)
  end do
!
! --- To cartesian component.
!
  pv = phig(ipg)
  deriv = dfdx_4th(phi5,fp5,pv)
!
end subroutine grdphi_gridpoint_type0
