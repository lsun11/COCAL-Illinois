module interface_outer_boundary_d_bvxd
  implicit none
  interface 
    subroutine outer_boundary_d_bvxd(sou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:)
    end subroutine outer_boundary_d_bvxd
  end interface
end module interface_outer_boundary_d_bvxd
