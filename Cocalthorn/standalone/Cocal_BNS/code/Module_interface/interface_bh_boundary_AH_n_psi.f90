module interface_bh_boundary_AH_n_psi
  implicit none
  interface 
    subroutine bh_boundary_AH_n_psi(dsou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: dsou_surf(:,:)
    end subroutine bh_boundary_AH_n_psi
  end interface
end module interface_bh_boundary_AH_n_psi
