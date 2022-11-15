module interface_sourceterm_MoC_CF_with_divshift_peos_spin
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_divshift_peos_spin(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_divshift_peos_spin
  end interface
end module interface_sourceterm_MoC_CF_with_divshift_peos_spin
