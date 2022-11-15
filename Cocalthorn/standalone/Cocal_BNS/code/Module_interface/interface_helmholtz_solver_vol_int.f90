module interface_helmholtz_solver_vol_int
  implicit none
  interface 
    subroutine helmholtz_solver_vol_int(sou,pot)
      real(8), pointer :: sou(:,:,:),pot(:,:,:)
    end subroutine helmholtz_solver_vol_int
  end interface
end module interface_helmholtz_solver_vol_int
