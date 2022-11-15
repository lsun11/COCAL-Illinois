module interface_vol_int_fluid
  implicit none
  interface 
    subroutine vol_int_fluid(souf,vol)
      real(8), pointer     :: souf(:,:,:)
      real(8), intent(out) :: vol
    end subroutine vol_int_fluid
  end interface
end module interface_vol_int_fluid
