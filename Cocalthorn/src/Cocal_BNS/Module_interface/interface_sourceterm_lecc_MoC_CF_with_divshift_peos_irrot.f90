module interface_sourceterm_lecc_MoC_CF_with_divshift_peos_irrot
  implicit none
  interface 
    subroutine sourceterm_lecc_MoC_CF_with_divshift_peos_irrot(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_lecc_MoC_CF_with_divshift_peos_irrot
  end interface
end module interface_sourceterm_lecc_MoC_CF_with_divshift_peos_irrot
