module interface_violation_gridpoint_HaC_CF_peos
  implicit none
  interface 
    subroutine violation_gridpoint_HaC_CF_peos(HaC_vio)
      real(8), pointer :: HaC_vio(:,:,:)
    end subroutine violation_gridpoint_HaC_CF_peos
  end interface
end module interface_violation_gridpoint_HaC_CF_peos
