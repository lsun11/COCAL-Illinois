module interface_sourceterm_lecc_MoC_CF_drot_SFC
  implicit none
  interface 
    subroutine sourceterm_lecc_MoC_CF_drot_SFC(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_lecc_MoC_CF_drot_SFC
  end interface
end module interface_sourceterm_lecc_MoC_CF_drot_SFC
