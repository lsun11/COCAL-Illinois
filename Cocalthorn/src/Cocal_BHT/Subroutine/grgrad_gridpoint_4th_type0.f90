subroutine grgrad_gridpoint_4th_type0(fnc,dfdx,dfdy,dfdz,irg,itg,ipg,cobj)
  use phys_constant, only : long
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_4th_gridpoint_bhex
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long) :: dfdx, dfdy, dfdz
  integer    :: irg, itg, ipg
  character(len=2), intent(in) :: cobj
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
! --- r, theta, phi derivatives.
!                                                                              
  if (cobj.eq.'bh') &
  &  call grgrad_4th_gridpoint_bhex(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
  if (cobj.eq.'ns') &
  &  call grgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
!
end subroutine grgrad_gridpoint_4th_type0
