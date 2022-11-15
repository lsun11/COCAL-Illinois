module interface_sourceterm_MoC_CF_with_divshift_peos
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_divshift_peos(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_divshift_peos
  end interface
end module interface_sourceterm_MoC_CF_with_divshift_peos
