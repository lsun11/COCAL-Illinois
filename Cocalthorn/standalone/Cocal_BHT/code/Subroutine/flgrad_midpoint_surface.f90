subroutine flgrad_midpoint_surface(fnc,dfdx,dfdy,dfdz)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: dfdx(:,:), dfdy(:,:), dfdz(:,:)
  real(long) :: dfdx0, dfdy0, dfdz0
  integer :: irf, itf, ipf
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives.
!
  do ipf = 1, npf
    do itf = 1, ntf
      irf = nrf
      call flgrad_midpoint_type0(fnc,dfdx0,dfdy0,dfdz0,irf,itf,ipf)
      dfdx(itf,ipf) = dfdx0
      dfdy(itf,ipf) = dfdy0
      dfdz(itf,ipf) = dfdz0
    end do
  end do
!
end subroutine flgrad_midpoint_surface
