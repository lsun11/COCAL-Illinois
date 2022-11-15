module interface_source_proper_mass_peos_irrot
  implicit none
  interface 
    subroutine source_proper_mass_peos_irrot(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_proper_mass_peos_irrot
  end interface
end module interface_source_proper_mass_peos_irrot
