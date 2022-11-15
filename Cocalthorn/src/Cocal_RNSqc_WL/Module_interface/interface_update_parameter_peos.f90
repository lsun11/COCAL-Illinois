module interface_update_parameter_peos
  implicit none
  interface 
    subroutine update_parameter_peos(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_peos
  end interface
end module interface_update_parameter_peos
