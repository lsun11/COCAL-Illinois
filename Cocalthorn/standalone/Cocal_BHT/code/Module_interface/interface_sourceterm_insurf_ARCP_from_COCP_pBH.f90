module interface_sourceterm_insurf_ARCP_from_COCP_pBH
  implicit none
  interface 
    subroutine sourceterm_insurf_ARCP_from_COCP_pBH(impt,fnchar,char,sou_in, &
    &                                                               dsou_in)
      real(8), pointer :: sou_in(:,:), dsou_in(:,:)
      integer :: impt
      character(len=4) :: fnchar
      character(len=2) :: char
    end subroutine sourceterm_insurf_ARCP_from_COCP_pBH
  end interface
end module interface_sourceterm_insurf_ARCP_from_COCP_pBH
