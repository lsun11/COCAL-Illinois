module interface_update_parameter_BNS_peos
  implicit none
  interface 
    subroutine update_parameter_BNS_peos(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS_peos
!
    subroutine update_parameter_BNS_peos_radi(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS_peos_radi
  end interface
end module interface_update_parameter_BNS_peos
