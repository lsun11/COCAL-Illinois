module interface_helmholtz_solver_binary_bhex_homosol
  implicit none
  interface 
    subroutine helmholtz_solver_binary_bhex_homosol(fnc,sou,pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:), fnc(:,:,:)
    end subroutine helmholtz_solver_binary_bhex_homosol
  end interface
end module interface_helmholtz_solver_binary_bhex_homosol
