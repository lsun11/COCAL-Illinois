module interface_update_parameter_axisym_qeos_drot
  implicit none
  interface 
    subroutine update_parameter_axisym_qeos_drot(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_axisym_qeos_drot
  end interface
end module interface_update_parameter_axisym_qeos_drot
