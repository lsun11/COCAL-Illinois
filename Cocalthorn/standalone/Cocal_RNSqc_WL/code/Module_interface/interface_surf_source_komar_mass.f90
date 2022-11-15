module interface_surf_source_komar_mass
  implicit none
  interface 
    subroutine surf_source_komar_mass(surf_sou_komar, irg)
      real(8), pointer  :: surf_sou_komar(:,:)
      integer ::  irg
    end subroutine surf_source_komar_mass
  end interface
end module interface_surf_source_komar_mass
