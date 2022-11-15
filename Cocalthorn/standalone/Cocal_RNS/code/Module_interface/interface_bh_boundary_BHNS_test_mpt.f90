module interface_bh_boundary_BHNS_test_mpt
  implicit none
  interface 
    subroutine bh_boundary_BHNS_test_mpt(char_bc,sou_surf,dsou_surf)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=2), intent(in) :: char_bc
    end subroutine bh_boundary_BHNS_test_mpt
  end interface
end module interface_bh_boundary_BHNS_test_mpt
