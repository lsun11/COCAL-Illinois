module interface_sourceterm_MoC_WL
  implicit none
  interface 
    subroutine sourceterm_MoC_WL(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL
  end interface
end module interface_sourceterm_MoC_WL
