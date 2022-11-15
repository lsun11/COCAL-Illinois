module interface_sourceterm_MoC_CF_drot
  implicit none
  interface 
    subroutine sourceterm_MoC_CF_drot(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_CF_drot
  end interface
end module interface_sourceterm_MoC_CF_drot
