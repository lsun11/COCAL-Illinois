module interface_source_komar_mass_compact_WL
  implicit none
  interface 
    subroutine source_komar_mass_compact_WL(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact_WL
  end interface
end module interface_source_komar_mass_compact_WL
