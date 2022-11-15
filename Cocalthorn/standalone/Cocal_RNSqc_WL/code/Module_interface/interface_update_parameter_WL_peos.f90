module interface_update_parameter_WL_peos
  implicit none
  interface 
    subroutine update_parameter_WL_peos(conv_den)
      real(8), intent(in) :: conv_den
    end subroutine update_parameter_WL_peos
  end interface
end module interface_update_parameter_WL_peos
