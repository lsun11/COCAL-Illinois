subroutine grgrad_gridpoint(fnc,dfdx,dfdy,dfdz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_4th_gridpoint_bhex
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: dfdx(:,:,:)
  real(long), pointer :: dfdy(:,:,:)
  real(long), pointer :: dfdz(:,:,:)
  real(long) :: dfncdx, dfncdy, dfncdz
  integer :: irg, itg, ipg
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
!
! --- r, theta, phi derivatives. 
!
  do irg = 0, nrg
    do itg = 0, ntg
      do ipg = 0, npg
        call grgrad_4th_gridpoint_bhex(fnc,dfncdx,dfncdy,dfncdz,irg,itg,ipg)
        dfdx(irg,itg,ipg) = dfncdx
        dfdy(irg,itg,ipg) = dfncdy
        dfdz(irg,itg,ipg) = dfncdz
      end do
    end do
  end do
end subroutine grgrad_gridpoint
