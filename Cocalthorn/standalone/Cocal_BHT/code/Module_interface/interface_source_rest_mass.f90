module interface_source_rest_mass
  implicit none
  interface 
    subroutine source_rest_mass(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_rest_mass
  end interface
end module interface_source_rest_mass
