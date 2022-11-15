module interface_sourceterm_MoC_WL_drot_SFC
  implicit none
  interface 
    subroutine sourceterm_MoC_WL_drot_SFC(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL_drot_SFC
  end interface
end module interface_sourceterm_MoC_WL_drot_SFC
