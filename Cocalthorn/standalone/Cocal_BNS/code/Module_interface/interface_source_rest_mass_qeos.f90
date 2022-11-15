module interface_source_rest_mass_qeos
  implicit none
  interface 
    subroutine source_rest_mass_qeos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_rest_mass_qeos
  end interface
end module interface_source_rest_mass_qeos
