module interface_bh_boundary_test
  implicit none
  interface 
    subroutine bh_boundary_test(char_bc,sou_surf,dsou_surf)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=2), intent(in) :: char_bc
    end subroutine bh_boundary_test
  end interface
end module interface_bh_boundary_test
