module interface_sourceterm_MoC_CF_corot
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_corot(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_corot
  end interface
end module interface_sourceterm_MoC_CF_corot
