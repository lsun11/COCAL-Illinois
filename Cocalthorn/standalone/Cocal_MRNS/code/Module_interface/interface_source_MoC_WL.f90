module interface_source_MoC_WL
  implicit none
  interface 
    subroutine source_MoC_WL(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_WL
  end interface
end module interface_source_MoC_WL
