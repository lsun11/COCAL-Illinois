module interface_source_adm_mass
  implicit none
  interface 
    subroutine source_adm_mass(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass
  end interface
end module interface_source_adm_mass
