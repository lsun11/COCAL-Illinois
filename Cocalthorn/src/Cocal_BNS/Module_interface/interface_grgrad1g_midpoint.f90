module interface_grgrad1g_midpoint
  implicit none
  interface 
    subroutine grgrad1g_midpoint(fnc,grad1,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: grad1(1:3)
      integer :: irg, itg, ipg
    end subroutine grgrad1g_midpoint
  end interface
end module interface_grgrad1g_midpoint
