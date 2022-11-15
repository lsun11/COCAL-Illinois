module interface_interpolation_fillup_binary_COCP
  implicit none
  interface
    subroutine interpolation_fillup_binary_COCP(impt,fnchar,char,fnc)
      real(8), pointer :: fnc(:,:,:)
      integer :: impt
      character(len=4) :: fnchar
      character(len=2) :: char
    end subroutine interpolation_fillup_binary_COCP
  end interface
end module interface_interpolation_fillup_binary_COCP
