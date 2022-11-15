subroutine flgrad(fnc,grad_fnc,order)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use interface_flgrad_2nd_gridpoint
  use interface_flgrad_4nd_gridpoint
  implicit none
  real(long) pointer :: fnc(:,:,:)
  real(long) pointer :: grad_fnc(:,:,:,:)
  real(long) :: dfdx, dfdy, dfdz
  integer :: irf, itf, ipf, order
!
  if (order.ne.2.or.order.ne.4) stop 'flgrad order wrong'
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
!
        if (order.eq.2) &
        & call flgrad_2nd_gridpoint(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
        if (order.eq.4) &
        & call flgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irf,itf,ipf)
        grad_fnc(irf,itf,ipf,1) = dfdx
        grad_fnc(irf,itf,ipf,2) = dfdy
        grad_fnc(irf,itf,ipf,3) = dfdz
!
      end do
    end do
  end do
end subroutine flgrad
