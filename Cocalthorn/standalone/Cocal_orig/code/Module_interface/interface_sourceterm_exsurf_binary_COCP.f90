module interface_sourceterm_exsurf_binary_COCP
  implicit none
  interface 
    subroutine sourceterm_exsurf_binary_COCP(impt,fnchar,char,sou_ex,dsou_ex)
      real(8), pointer :: sou_ex(:,:), dsou_ex(:,:)
      integer :: impt
      character(len=4) :: fnchar
      character(len=2) :: char
    end subroutine sourceterm_exsurf_binary_COCP
  end interface
end module interface_sourceterm_exsurf_binary_COCP
