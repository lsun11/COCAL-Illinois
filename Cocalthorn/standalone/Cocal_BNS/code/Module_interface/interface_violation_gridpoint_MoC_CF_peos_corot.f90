module interface_violation_gridpoint_MoC_CF_peos_corot
  implicit none
  interface 
    subroutine violation_gridpoint_MoC_CF_peos_corot(MoC_vio,cobj)
      real(8), pointer :: MoC_vio(:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine violation_gridpoint_MoC_CF_peos_corot
  end interface
end module interface_violation_gridpoint_MoC_CF_peos_corot
