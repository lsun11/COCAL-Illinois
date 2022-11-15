module interface_error_surface
  implicit none
  interface 
    subroutine error_surface(pot,pot_bak,error,flag)
      real(8), pointer     :: pot(:,:), pot_bak(:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_surface
  end interface
end module interface_error_surface
