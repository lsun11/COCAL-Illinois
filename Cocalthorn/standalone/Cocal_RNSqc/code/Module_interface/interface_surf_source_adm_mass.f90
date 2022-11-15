module interface_surf_source_adm_mass
  implicit none
  interface 
    subroutine surf_source_adm_mass(surf_sou_adm, irg)
      real(8), pointer  :: surf_sou_adm(:,:)
      integer :: irg
    end subroutine surf_source_adm_mass
  end interface
end module interface_surf_source_adm_mass
