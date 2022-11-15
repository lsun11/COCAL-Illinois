module interface_surf_int_grav_rg_noweights
  implicit none
  interface 
    subroutine surf_int_grav_rg_noweights(sousg,surf)
      real(8), pointer     :: sousg(:,:)
      real(8), intent(out) :: surf
    end subroutine surf_int_grav_rg_noweights
  end interface
end module interface_surf_int_grav_rg_noweights
