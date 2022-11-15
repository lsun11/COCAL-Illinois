module interface_sourceterm_surface_int
  implicit none
  interface 
    subroutine sourceterm_surface_int(fnc,irg,sou_surf,dsou_surf)
      real(8), pointer :: fnc(:,:,:), sou_surf(:,:), dsou_surf(:,:)
      integer :: irg
    end subroutine sourceterm_surface_int
  end interface
end module interface_sourceterm_surface_int
