module interface_plane_surf_int_fluid
  implicit none
  interface 
    subroutine plane_surf_int_fluid(cline, ia, ib, souf, surf_int)
      character(len=2), intent(in) :: cline
      integer,          intent(in) :: ia, ib
      real(8), pointer             :: souf(:,:)
      real(8), intent(out)         :: surf_int
    end subroutine plane_surf_int_fluid
  end interface
end module interface_plane_surf_int_fluid
