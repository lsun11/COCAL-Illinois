module interface_source_mass_asympto
  implicit none
  interface 
    subroutine source_mass_asympto(fnc,sousf,irg)
      real(8), pointer :: fnc(:,:,:), sousf(:,:)
      integer          :: irg
    end subroutine source_mass_asympto
  end interface
end module interface_source_mass_asympto
