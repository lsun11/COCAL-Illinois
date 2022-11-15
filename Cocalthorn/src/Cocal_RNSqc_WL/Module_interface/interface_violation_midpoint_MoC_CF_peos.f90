module interface_violation_midpoint_MoC_CF_peos
  implicit none
  interface 
    subroutine violation_midpoint_MoC_CF_peos(MoC_vio)
      real(8), pointer :: MoC_vio(:,:,:,:)
    end subroutine violation_midpoint_MoC_CF_peos
  end interface
end module interface_violation_midpoint_MoC_CF_peos
