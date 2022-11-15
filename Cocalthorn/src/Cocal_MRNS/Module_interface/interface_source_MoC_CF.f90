module interface_source_MoC_CF
  implicit none
  interface 
    subroutine source_MoC_CF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_CF
  end interface
end module interface_source_MoC_CF
