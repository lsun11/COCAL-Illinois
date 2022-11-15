module interface_helmholtz_solver_asymptotic_patch_homosol_outgoing
  implicit none
  interface 
    subroutine helmholtz_solver_asymptotic_patch_homosol_outgoing &
    &                             (sou,sou_insurf,dsou_insurf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_insurf(:,:), dsou_insurf(:,:)
    end subroutine helmholtz_solver_asymptotic_patch_homosol_outgoing
  end interface
end module interface_helmholtz_solver_asymptotic_patch_homosol_outgoing
