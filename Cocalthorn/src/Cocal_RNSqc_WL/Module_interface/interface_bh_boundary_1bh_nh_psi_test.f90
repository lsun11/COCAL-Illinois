module interface_bh_boundary_1bh_nh_psi_test
  implicit none
  interface
    subroutine bh_boundary_1bh_nh_psi_test(dsou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: dsou_surf(:,:)
    end subroutine bh_boundary_1bh_nh_psi_test
  end interface
end module interface_bh_boundary_1bh_nh_psi_test
