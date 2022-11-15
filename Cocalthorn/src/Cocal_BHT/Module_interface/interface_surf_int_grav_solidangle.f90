module interface_surf_int_grav_solidangle
  implicit none
  interface 
    subroutine surf_int_grav_solidangle(sousf,surf)
      real(8), pointer     :: sousf(:,:)
      real(8), intent(out) :: surf
    end subroutine surf_int_grav_solidangle
  end interface
end module interface_surf_int_grav_solidangle
