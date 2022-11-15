module interface_source_proper_mass
  implicit none
  interface 
    subroutine source_proper_mass(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_proper_mass
  end interface
end module interface_source_proper_mass
