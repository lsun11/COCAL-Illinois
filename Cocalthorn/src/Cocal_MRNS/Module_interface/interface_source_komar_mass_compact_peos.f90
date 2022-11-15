module interface_source_komar_mass_compact_peos
  implicit none
  interface 
    subroutine source_komar_mass_compact_peos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact_peos
  end interface
end module interface_source_komar_mass_compact_peos
