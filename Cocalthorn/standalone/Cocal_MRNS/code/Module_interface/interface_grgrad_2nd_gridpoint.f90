module interface_grgrad_2nd_gridpoint
  implicit none
  interface 
    subroutine grgrad_2nd_gridpoint(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_2nd_gridpoint
  end interface
end module interface_grgrad_2nd_gridpoint
