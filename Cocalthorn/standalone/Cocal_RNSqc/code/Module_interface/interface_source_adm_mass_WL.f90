module interface_source_adm_mass_WL
  implicit none
  interface 
    subroutine source_adm_mass_WL(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass_WL
  end interface
end module interface_source_adm_mass_WL
