module interface_source_MoC_WL_BH
  implicit none
  interface 
    subroutine source_MoC_WL_BH(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine source_MoC_WL_BH
  end interface
end module interface_source_MoC_WL_BH
