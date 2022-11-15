module interface_source_adm_mass_peos
  implicit none
  interface 
    subroutine source_adm_mass_peos(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass_peos
  end interface
end module interface_source_adm_mass_peos
