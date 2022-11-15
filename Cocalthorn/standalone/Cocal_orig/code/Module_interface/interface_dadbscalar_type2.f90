module interface_dadbscalar_type2
  implicit none
  interface 
    subroutine dadbscalar_type2(fnc,dadbfnc,cobj)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: dadbfnc(:,:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine dadbscalar_type2
  end interface
end module interface_dadbscalar_type2
