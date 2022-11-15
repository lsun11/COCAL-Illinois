module interface_bh_boundary_AH_n_alps
  implicit none
  interface 
    subroutine bh_boundary_AH_n_alps(dsou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: dsou_surf(:,:)
    end subroutine bh_boundary_AH_n_alps
  end interface
end module interface_bh_boundary_AH_n_alps
