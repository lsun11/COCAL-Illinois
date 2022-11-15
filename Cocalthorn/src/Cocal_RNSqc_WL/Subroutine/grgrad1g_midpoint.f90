subroutine grgrad1g_midpoint(fnc,grad1,irg,itg,ipg)
  use interface_grgrad_midpoint_type0
  implicit none
!
  real(8), pointer :: fnc(:,:,:)
  real(8) :: grad1(3), dfdx, dfdy, dfdz
  integer :: irg, itg, ipg
!
  call grgrad_midpoint_type0(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
  grad1(1) = dfdx
  grad1(2) = dfdy
  grad1(3) = dfdz
!
end subroutine grgrad1g_midpoint
