module interface_update_parameter_BNS_peos_lecc
  implicit none
  interface 
    subroutine update_parameter_BNS_peos_lecc(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS_peos_lecc
  end interface
end module interface_update_parameter_BNS_peos_lecc
