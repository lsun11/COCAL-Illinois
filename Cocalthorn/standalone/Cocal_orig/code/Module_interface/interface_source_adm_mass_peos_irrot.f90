module interface_source_adm_mass_peos_irrot
  implicit none
  interface 
    subroutine source_adm_mass_peos_irrot(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass_peos_irrot
  end interface
end module interface_source_adm_mass_peos_irrot
