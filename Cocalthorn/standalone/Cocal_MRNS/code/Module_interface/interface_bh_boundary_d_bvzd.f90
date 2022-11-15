module interface_bh_boundary_d_bvzd
  implicit none
  interface 
    subroutine bh_boundary_d_bvzd(sou_surf)
!      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
!      character(len=2), intent(in) :: char_bc
      real(8), pointer :: sou_surf(:,:)
    end subroutine bh_boundary_d_bvzd
  end interface
end module interface_bh_boundary_d_bvzd
