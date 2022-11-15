module interface_sourceterm_MoC_CF_type1_bhex
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_type1_bhex(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_type1_bhex
  end interface
end module interface_sourceterm_MoC_CF_type1_bhex
