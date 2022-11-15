module interface_sourceterm_MoC_WL_SEM
  implicit none
  interface 
    subroutine sourceterm_MoC_WL_SEM(souvec)
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC_WL_SEM
  end interface
end module interface_sourceterm_MoC_WL_SEM
