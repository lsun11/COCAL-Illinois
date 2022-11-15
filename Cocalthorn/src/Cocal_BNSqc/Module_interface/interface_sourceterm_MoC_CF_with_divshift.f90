module interface_sourceterm_MoC_CF_with_divshift
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_divshift(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_divshift
  end interface
end module interface_sourceterm_MoC_CF_with_divshift
