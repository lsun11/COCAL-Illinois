module interface_sourceterm_MoC_CF
  implicit none
  interface 
    subroutine sourceterm_MoC_CF(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF
  end interface
end module interface_sourceterm_MoC_CF
