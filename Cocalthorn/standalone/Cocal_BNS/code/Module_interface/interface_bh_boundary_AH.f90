module interface_bh_boundary_AH
  implicit none
  interface 
    subroutine bh_boundary_AH(sou_surf,dsou_surf,char_mp)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=4), intent(in) :: char_mp
    end subroutine bh_boundary_AH
  end interface
end module interface_bh_boundary_AH
