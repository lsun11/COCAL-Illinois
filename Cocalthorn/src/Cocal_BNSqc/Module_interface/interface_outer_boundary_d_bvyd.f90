module interface_outer_boundary_d_bvyd
  implicit none
  interface 
    subroutine outer_boundary_d_bvyd(sou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:)
    end subroutine outer_boundary_d_bvyd
  end interface
end module interface_outer_boundary_d_bvyd
