module interface_helmholtz_solver
  implicit none
  interface 
    subroutine helmholtz_solver(sou,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
    end subroutine helmholtz_solver
  end interface
end module interface_helmholtz_solver
