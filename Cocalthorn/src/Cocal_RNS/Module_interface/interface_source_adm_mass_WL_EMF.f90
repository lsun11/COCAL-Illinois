module interface_source_adm_mass_WL_EMF
  implicit none
  interface 
    subroutine source_adm_mass_WL_EMF(soug)
      real(8), pointer     :: soug(:,:,:)
    end subroutine source_adm_mass_WL_EMF
  end interface
end module interface_source_adm_mass_WL_EMF
