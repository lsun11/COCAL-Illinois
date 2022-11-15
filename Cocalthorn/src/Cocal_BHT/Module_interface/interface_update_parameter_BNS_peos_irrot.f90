module interface_update_parameter_BNS_peos_irrot
  implicit none
  interface 
    subroutine update_parameter_BNS_peos_irrot(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS_peos_irrot
  end interface
end module interface_update_parameter_BNS_peos_irrot
