module interface_update_parameter_triaxial
  implicit none
  interface 
    subroutine update_parameter_triaxial(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_triaxial
  end interface
end module interface_update_parameter_triaxial
