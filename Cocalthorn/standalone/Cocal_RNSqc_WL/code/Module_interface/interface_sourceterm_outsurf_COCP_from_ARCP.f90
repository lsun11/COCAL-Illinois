module interface_sourceterm_outsurf_COCP_from_ARCP
  implicit none
  interface 
    subroutine sourceterm_outsurf_COCP_from_ARCP(impt,fnchar,char,sou_out, &
    &                                                            dsou_out)
      real(8), pointer :: sou_out(:,:), dsou_out(:,:)
      integer :: impt
      character(len=4) :: fnchar
      character(len=2) :: char
    end subroutine sourceterm_outsurf_COCP_from_ARCP
  end interface
end module interface_sourceterm_outsurf_COCP_from_ARCP
