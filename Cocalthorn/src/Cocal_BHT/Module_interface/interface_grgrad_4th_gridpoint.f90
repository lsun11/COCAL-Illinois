module interface_grgrad_4th_gridpoint
  implicit none
  interface
    subroutine grgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_4th_gridpoint
  end interface
end module interface_grgrad_4th_gridpoint
