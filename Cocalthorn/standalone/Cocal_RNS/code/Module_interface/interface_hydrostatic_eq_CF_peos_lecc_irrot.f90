module interface_hydrostatic_eq_CF_peos_lecc_irrot
  implicit none
  interface 
    subroutine hydrostatic_eq_CF_peos_lecc_irrot(emd,utf,vepxf,vepyf,vepzf)
      real(8), pointer :: emd(:,:,:),  utf(:,:,:), vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
    end subroutine hydrostatic_eq_CF_peos_lecc_irrot
  end interface
end module interface_hydrostatic_eq_CF_peos_lecc_irrot
