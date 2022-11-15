module interface_sourceterm_MoC_WL_drot_all
  implicit none
  interface 
    subroutine sourceterm_MoC_WL_drot_all(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL_drot_all
  end interface
end module interface_sourceterm_MoC_WL_drot_all
