module interface_grgrad_midpoint_type0
  implicit none
  interface 
    subroutine grgrad_midpoint_type0(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_midpoint_type0
  end interface
end module interface_grgrad_midpoint_type0
