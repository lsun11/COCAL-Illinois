module interface_source_komar_mass
  implicit none
  interface 
    subroutine source_komar_mass(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_komar_mass
  end interface
end module interface_source_komar_mass
