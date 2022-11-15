module interface_source_komar_mass_compact_qeos
  implicit none
  interface 
    subroutine source_komar_mass_compact_qeos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact_qeos
  end interface
end module interface_source_komar_mass_compact_qeos
