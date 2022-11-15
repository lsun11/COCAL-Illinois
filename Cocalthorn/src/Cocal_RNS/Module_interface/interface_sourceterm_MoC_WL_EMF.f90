module interface_sourceterm_MoC_WL_EMF
  implicit none
  interface 
    subroutine sourceterm_MoC_WL_EMF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL_EMF
  end interface
end module interface_sourceterm_MoC_WL_EMF
