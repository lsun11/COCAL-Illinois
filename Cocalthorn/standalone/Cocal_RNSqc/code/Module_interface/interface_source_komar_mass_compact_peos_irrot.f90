module interface_source_komar_mass_compact_peos_irrot
  implicit none
  interface 
    subroutine source_komar_mass_compact_peos_irrot(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact_peos_irrot
  end interface
end module interface_source_komar_mass_compact_peos_irrot
