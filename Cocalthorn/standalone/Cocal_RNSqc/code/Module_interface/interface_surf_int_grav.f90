module interface_surf_int_grav
  implicit none
  interface 
    subroutine surf_int_grav(sousf,surf,irg)
      real(8), pointer     :: sousf(:,:)
      real(8), intent(out) :: surf
      integer, intent(in)  :: irg
    end subroutine surf_int_grav
  end interface
end module interface_surf_int_grav
