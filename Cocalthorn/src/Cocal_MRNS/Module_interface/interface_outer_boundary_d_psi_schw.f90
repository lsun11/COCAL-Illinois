module interface_outer_boundary_d_psi_schw
  implicit none
  interface 
    subroutine outer_boundary_d_psi_schw(sou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:)
    end subroutine outer_boundary_d_psi_schw
  end interface
end module interface_outer_boundary_d_psi_schw
