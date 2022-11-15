module interface_sourceterm_MoC_WL_CTT_all_bhex
  implicit none
  interface 
    subroutine sourceterm_MoC_WL_CTT_all_bhex(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL_CTT_all_bhex
  end interface
end module interface_sourceterm_MoC_WL_CTT_all_bhex
