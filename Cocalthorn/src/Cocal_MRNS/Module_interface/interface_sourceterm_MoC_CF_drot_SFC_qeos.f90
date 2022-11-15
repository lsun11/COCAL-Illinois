module interface_sourceterm_MoC_CF_drot_SFC_qeos
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_drot_SFC_qeos(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_drot_SFC_qeos
  end interface
end module interface_sourceterm_MoC_CF_drot_SFC_qeos
