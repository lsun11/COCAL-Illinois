module interface_sourceterm_MoC_CF_with_2nd_derivatives
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_2nd_derivatives(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_2nd_derivatives
  end interface
end module interface_sourceterm_MoC_CF_with_2nd_derivatives
