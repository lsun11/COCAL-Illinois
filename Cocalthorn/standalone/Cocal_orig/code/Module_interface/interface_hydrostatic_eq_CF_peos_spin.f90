module interface_hydrostatic_eq_CF_peos_spin
  implicit none
  interface 
    subroutine hydrostatic_eq_CF_peos_spin(emd,utf,vepxf,vepyf,vepzf)
      real(8), pointer :: emd(:,:,:),  utf(:,:,:), vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
    end subroutine hydrostatic_eq_CF_peos_spin
  end interface
end module interface_hydrostatic_eq_CF_peos_spin
