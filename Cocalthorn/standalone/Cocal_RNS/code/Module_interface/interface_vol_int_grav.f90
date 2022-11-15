module interface_vol_int_grav
  implicit none
  interface 
    subroutine vol_int_grav(soug,vol)
      real(8), pointer     :: soug(:,:,:)
      real(8), intent(out) :: vol
    end subroutine vol_int_grav
  end interface
end module interface_vol_int_grav
