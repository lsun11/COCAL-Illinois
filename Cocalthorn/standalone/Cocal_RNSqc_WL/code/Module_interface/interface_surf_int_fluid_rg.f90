module interface_surf_int_fluid_rg
  implicit none
  interface 
    subroutine surf_int_fluid_rg(souf,surf,irf)
      real(8), pointer     :: souf(:,:)
      real(8), intent(out) :: surf
      integer, intent(in)  :: irf
    end subroutine surf_int_fluid_rg
  end interface
end module interface_surf_int_fluid_rg
