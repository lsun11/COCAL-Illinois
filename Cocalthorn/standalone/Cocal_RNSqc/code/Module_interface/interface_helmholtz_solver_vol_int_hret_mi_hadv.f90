module interface_helmholtz_solver_vol_int_hret_mi_hadv
  implicit none
  interface 
    subroutine helmholtz_solver_vol_int_hret_mi_hadv(sou,pot)
      real(8), pointer :: sou(:,:,:),pot(:,:,:)
    end subroutine helmholtz_solver_vol_int_hret_mi_hadv
  end interface
end module interface_helmholtz_solver_vol_int_hret_mi_hadv
