module interface_update_parameter_axisym_qeos
  implicit none
  interface 
    subroutine update_parameter_axisym_qeos(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_axisym_qeos
  end interface
end module interface_update_parameter_axisym_qeos
