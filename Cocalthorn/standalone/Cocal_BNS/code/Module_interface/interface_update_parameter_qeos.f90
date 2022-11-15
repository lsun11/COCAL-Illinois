module interface_update_parameter_qeos
  implicit none
  interface 
    subroutine update_parameter_qeos(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_qeos
  end interface
end module interface_update_parameter_qeos
