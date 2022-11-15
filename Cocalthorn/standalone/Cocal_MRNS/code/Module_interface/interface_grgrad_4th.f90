module interface_grgrad_4th
  implicit none
  interface 
    subroutine grgrad_4th(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8), intent(out) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_4th
  end interface
end module interface_grgrad_4th
