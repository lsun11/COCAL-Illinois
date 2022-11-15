module interface_interpo_radial1p_grav
  implicit none
  interface 
    subroutine interpo_radial1p_grav(grv,val,rv,it,ip)
      real(8), pointer :: grv(:,:,:)
      real(8), intent(out) :: val
      real(8), intent(in)  :: rv
      integer, intent(in)  :: it, ip
    end subroutine interpo_radial1p_grav
  end interface
end module interface_interpo_radial1p_grav
