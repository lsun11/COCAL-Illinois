module interface_vol_int_grav_bhex
  implicit none
  interface 
    subroutine vol_int_grav_bhex(sou,vol,irg_vol)
      real(8), pointer     :: sou(:,:,:)
      real(8), intent(out) :: vol
      integer, intent(in)  :: irg_vol
    end subroutine vol_int_grav_bhex
  end interface
end module interface_vol_int_grav_bhex
