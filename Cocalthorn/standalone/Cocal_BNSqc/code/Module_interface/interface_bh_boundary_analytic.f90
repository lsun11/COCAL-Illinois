module interface_bh_boundary_analytic
  implicit none
  interface 
    subroutine bh_boundary_analytic(char_mp,irg,sou_surf,dsou_surf)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
      character(len=4), intent(in) :: char_mp
      integer :: irg
    end subroutine bh_boundary_analytic
  end interface
end module interface_bh_boundary_analytic
