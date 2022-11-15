module interface_sourceterm_MoC_CF_with_divshift_r3rd
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_divshift_r3rd(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_divshift_r3rd
  end interface
end module interface_sourceterm_MoC_CF_with_divshift_r3rd
