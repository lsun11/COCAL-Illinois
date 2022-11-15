module interface_update_parameter_axisym_peos
  implicit none
  interface 
    subroutine update_parameter_axisym_peos(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_axisym_peos
  end interface
end module interface_update_parameter_axisym_peos
