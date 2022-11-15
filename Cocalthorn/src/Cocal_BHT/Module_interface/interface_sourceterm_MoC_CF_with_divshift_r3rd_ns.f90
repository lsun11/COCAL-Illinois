module interface_sourceterm_MoC_CF_with_divshift_r3rd_ns
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_with_divshift_r3rd_ns(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_with_divshift_r3rd_ns
  end interface
end module interface_sourceterm_MoC_CF_with_divshift_r3rd_ns
