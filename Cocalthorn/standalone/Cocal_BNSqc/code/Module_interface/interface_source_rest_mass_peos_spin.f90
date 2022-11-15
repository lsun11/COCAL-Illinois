module interface_source_rest_mass_peos_spin
  implicit none
  interface 
    subroutine source_rest_mass_peos_spin(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_rest_mass_peos_spin
  end interface
end module interface_source_rest_mass_peos_spin
