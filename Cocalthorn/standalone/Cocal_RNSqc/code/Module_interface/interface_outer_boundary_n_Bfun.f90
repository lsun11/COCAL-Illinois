module interface_outer_boundary_n_Bfun
  implicit none
  interface 
    subroutine outer_boundary_n_Bfun(dsou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: dsou_surf(:,:)
    end subroutine outer_boundary_n_Bfun
  end interface
end module interface_outer_boundary_n_Bfun
