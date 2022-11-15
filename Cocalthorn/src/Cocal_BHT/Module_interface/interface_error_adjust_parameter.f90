module interface_error_adjust_parameter
  implicit none
  interface 
    subroutine error_adjust_parameter(niq,pot,pot_bak,error,flag)
      integer              :: niq
      real(8), intent(in)  :: pot(niq), pot_bak(niq)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_adjust_parameter
  end interface
end module interface_error_adjust_parameter
