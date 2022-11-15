module interface_error_velocity_potential
  implicit none
  interface 
    subroutine error_velocity_potential(pot,pot_bak,error,flag) 
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_velocity_potential
  end interface
end module interface_error_velocity_potential
