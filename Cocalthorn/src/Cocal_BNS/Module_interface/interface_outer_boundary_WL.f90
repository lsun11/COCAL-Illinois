module interface_outer_boundary_WL
  implicit none
  interface 
    subroutine outer_boundary_WL(sou_surf,dsou_surf,char_mp)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=4), intent(in) :: char_mp
    end subroutine outer_boundary_WL
  end interface
end module interface_outer_boundary_WL
