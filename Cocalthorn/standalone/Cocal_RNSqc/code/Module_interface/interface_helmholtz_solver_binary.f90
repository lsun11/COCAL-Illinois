module interface_helmholtz_solver_binary
  implicit none
  interface 
    subroutine helmholtz_solver_binary(fnc,sou,sou_exsurf,dsou_exsurf,sou_outsurf,dsou_outsurf,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:), fnc(:,:,:)
      real(8), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
      real(8), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
    end subroutine helmholtz_solver_binary
  end interface
end module interface_helmholtz_solver_binary
