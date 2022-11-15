module interface_update_parameter_BNS
  implicit none
  interface 
    subroutine update_parameter_BNS(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_BNS
  end interface
end module interface_update_parameter_BNS
