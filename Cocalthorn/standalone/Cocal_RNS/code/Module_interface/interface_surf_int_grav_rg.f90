module interface_surf_int_grav_rg
  implicit none
  interface 
    subroutine surf_int_grav_rg(sousf,surf,irg)
      real(8), pointer     :: sousf(:,:)
      real(8), intent(out) :: surf
      integer, intent(in)  :: irg
    end subroutine surf_int_grav_rg
  end interface
end module interface_surf_int_grav_rg
