module interface_source_rest_mass_peos
  implicit none
  interface 
    subroutine source_rest_mass_peos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_rest_mass_peos
  end interface
end module interface_source_rest_mass_peos
