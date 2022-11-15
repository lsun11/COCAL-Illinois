module interface_source_MoC_CF_gp4th
  implicit none
  interface 
    subroutine source_MoC_CF_gp4th(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_CF_gp4th
  end interface
end module interface_source_MoC_CF_gp4th
