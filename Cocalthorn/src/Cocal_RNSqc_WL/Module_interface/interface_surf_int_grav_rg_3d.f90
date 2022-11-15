module interface_surf_int_grav_rg_3d
  implicit none
  interface 
    subroutine surf_int_grav_rg_3d(sousf,surfx,surfy,surfz,irg)
      real(8), pointer     :: sousf(:,:,:)
      real(8), intent(out) :: surfx,surfy,surfz
      integer, intent(in)  :: irg
    end subroutine surf_int_grav_rg_3d
  end interface
end module interface_surf_int_grav_rg_3d
