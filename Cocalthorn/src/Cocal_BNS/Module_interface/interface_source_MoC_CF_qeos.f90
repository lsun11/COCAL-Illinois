module interface_source_MoC_CF_qeos
  implicit none
  interface 
    subroutine source_MoC_CF_qeos(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_CF_qeos
  end interface
end module interface_source_MoC_CF_qeos
