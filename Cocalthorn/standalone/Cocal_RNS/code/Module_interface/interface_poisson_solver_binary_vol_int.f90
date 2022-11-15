module interface_poisson_solver_binary_vol_int
  implicit none
  interface 
    subroutine poisson_solver_binary_vol_int(sou,pot)
      real(8), pointer :: sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_binary_vol_int
  end interface
end module interface_poisson_solver_binary_vol_int
