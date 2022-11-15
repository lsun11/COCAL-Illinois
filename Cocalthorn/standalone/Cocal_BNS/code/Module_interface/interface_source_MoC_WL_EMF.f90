module interface_source_MoC_WL_EMF
  implicit none
  interface 
    subroutine source_MoC_WL_EMF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_WL_EMF
  end interface
end module interface_source_MoC_WL_EMF
