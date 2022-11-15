module interface_violation_gridpoint_HaC_CF_peos_irrot
  implicit none
  interface 
    subroutine violation_gridpoint_HaC_CF_peos_irrot(HaC_vio,cobj)
      real(8), pointer :: HaC_vio(:,:,:)
      character(len=2), intent(in) :: cobj
    end subroutine violation_gridpoint_HaC_CF_peos_irrot
  end interface
end module interface_violation_gridpoint_HaC_CF_peos_irrot
