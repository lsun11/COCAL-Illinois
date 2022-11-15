module interface_update_parameter_BNS_peos_lecc_spin
  implicit none
  interface 
    subroutine update_parameter_BNS_peos_lecc_spin(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS_peos_lecc_spin
  end interface
end module interface_update_parameter_BNS_peos_lecc_spin
