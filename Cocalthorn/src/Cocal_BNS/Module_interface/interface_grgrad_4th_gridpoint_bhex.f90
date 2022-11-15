module interface_grgrad_4th_gridpoint_bhex
  implicit none
  interface
    subroutine grgrad_4th_gridpoint_bhex(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_4th_gridpoint_bhex
  end interface
end module interface_grgrad_4th_gridpoint_bhex
