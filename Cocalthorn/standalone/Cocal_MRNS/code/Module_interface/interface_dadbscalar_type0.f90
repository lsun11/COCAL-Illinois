module interface_dadbscalar_type0
  implicit none
  interface 
    subroutine dadbscalar_type0(fnc,d2f,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: d2f(1:3,1:3)
      integer :: irg, itg, ipg
    end subroutine dadbscalar_type0
  end interface
end module interface_dadbscalar_type0
