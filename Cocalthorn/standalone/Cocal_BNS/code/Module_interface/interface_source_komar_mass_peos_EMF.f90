module interface_source_komar_mass_peos_EMF
  implicit none
  interface 
    subroutine source_komar_mass_peos_EMF(soug)
      real(8), pointer     :: soug(:,:,:)
    end subroutine source_komar_mass_peos_EMF
  end interface
end module interface_source_komar_mass_peos_EMF
