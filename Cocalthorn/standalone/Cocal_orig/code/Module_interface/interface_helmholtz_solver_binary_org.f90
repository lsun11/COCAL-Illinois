module interface_helmholtz_solver_binary
  implicit none
  interface 
    subroutine helmholtz_solver_binary(sou,sou_surf,dsou_surf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
      real(8), pointer :: sou_surf(:,:), dsou_surf(:,:)
    end subroutine helmholtz_solver_binary
  end interface
end module interface_helmholtz_solver_binary
