module interface_grgrad_3rd
  implicit none
  interface 
    subroutine grgrad_3rd(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_3rd
  end interface
end module interface_grgrad_3rd
