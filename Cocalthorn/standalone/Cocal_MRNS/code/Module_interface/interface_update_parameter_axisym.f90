module interface_update_parameter_axisym
  implicit none
  interface 
    subroutine update_parameter_axisym(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_axisym
  end interface
end module interface_update_parameter_axisym
