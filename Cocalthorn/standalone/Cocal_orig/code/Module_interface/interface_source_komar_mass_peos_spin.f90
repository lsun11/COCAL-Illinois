module interface_source_komar_mass_peos_spin
  implicit none
  interface 
    subroutine source_komar_mass_peos_spin(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_komar_mass_peos_spin
  end interface
end module interface_source_komar_mass_peos_spin
