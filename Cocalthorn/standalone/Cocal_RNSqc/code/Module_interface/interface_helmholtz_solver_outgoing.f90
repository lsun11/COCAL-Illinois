module interface_helmholtz_solver_outgoing
  implicit none
  interface 
    subroutine helmholtz_solver_outgoing(sou,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
    end subroutine helmholtz_solver_outgoing
  end interface
end module interface_helmholtz_solver_outgoing
