module interface_dadbscalar_3rd2nd_type0
  implicit none
  interface 
    subroutine dadbscalar_3rd2nd_type0(fnc,d2f,irg,itg,ipg,cobj)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: d2f(1:3,1:3)
      integer :: irg, itg, ipg
      character(len=2), intent(in) :: cobj
    end subroutine dadbscalar_3rd2nd_type0
  end interface
end module interface_dadbscalar_3rd2nd_type0
