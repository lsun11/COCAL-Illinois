module interface_violation_midpoint_MoC_CF_peos_irrot
  implicit none
  interface 
    subroutine violation_midpoint_MoC_CF_peos_irrot(MoC_vio,cobj)
      real(8), pointer :: MoC_vio(:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine violation_midpoint_MoC_CF_peos_irrot
  end interface
end module interface_violation_midpoint_MoC_CF_peos_irrot
