module interface_bh_boundary_AH_d_psi
  implicit none
  interface 
    subroutine bh_boundary_AH_d_psi(sou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:)
    end subroutine bh_boundary_AH_d_psi
  end interface
end module interface_bh_boundary_AH_d_psi
