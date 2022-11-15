module interface_source_komar_mass_compact
  implicit none
  interface 
    subroutine source_komar_mass_compact(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact
  end interface
end module interface_source_komar_mass_compact
