module interface_sourceterm_insurf_ARCP_from_COCP
  implicit none
  interface 
    subroutine sourceterm_insurf_ARCP_from_COCP(impt,fnchar,char,sou_in, &
    &                                                           dsou_in)
      real(8), pointer :: sou_in(:,:), dsou_in(:,:)
      integer :: impt
      character(len=4) :: fnchar
      character(len=2) :: char
    end subroutine sourceterm_insurf_ARCP_from_COCP
  end interface
end module interface_sourceterm_insurf_ARCP_from_COCP
