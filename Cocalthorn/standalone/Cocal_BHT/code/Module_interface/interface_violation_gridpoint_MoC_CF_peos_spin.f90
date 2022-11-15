module interface_violation_gridpoint_MoC_CF_peos_spin
  implicit none
  interface 
    subroutine violation_gridpoint_MoC_CF_peos_spin(MoC_vio,cobj)
      real(8), pointer :: MoC_vio(:,:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine violation_gridpoint_MoC_CF_peos_spin
  end interface
end module interface_violation_gridpoint_MoC_CF_peos_spin
