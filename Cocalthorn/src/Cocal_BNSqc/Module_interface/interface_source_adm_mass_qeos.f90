module interface_source_adm_mass_qeos
  implicit none
  interface 
    subroutine source_adm_mass_qeos(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass_qeos
  end interface
end module interface_source_adm_mass_qeos
