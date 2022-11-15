module interface_update_parameter
  implicit none
  interface 
    subroutine update_parameter(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter
  end interface
end module interface_update_parameter
