subroutine grgrad_midpoint_r3rd(fnc,dfdx,dfdy,dfdz,cobj)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_midpoint_r3rd_type0
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(long) :: dfncdx, dfncdy, dfncdz
  integer    :: irg, itg, ipg
  character(len=2), intent(in) :: cobj
!
! --- Compute the gradient of a function.
! --- The gradient is evaluated at mid points.
! --- r, theta, phi derivatives.
!
  do irg = 1, nrg
    do itg = 1, ntg
      do ipg = 1, npg
!
        call grgrad_midpoint_r3rd_type0(fnc,dfncdx,dfncdy,dfncdz, &
        &                               irg,itg,ipg,cobj)
!
        dfdx(irg,itg,ipg) = dfncdx
        dfdy(irg,itg,ipg) = dfncdy
        dfdz(irg,itg,ipg) = dfncdz
!
      end do 
    end do 
  end do
!
end subroutine grgrad_midpoint_r3rd
