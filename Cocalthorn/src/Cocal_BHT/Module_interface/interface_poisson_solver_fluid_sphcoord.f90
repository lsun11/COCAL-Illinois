module interface_poisson_solver_fluid_sphcoord
  implicit none
  interface 
    subroutine poisson_solver_fluid_sphcoord(sou, pot)
      real(8), pointer :: pot(:,:,:), sou(:,:,:)
    end subroutine poisson_solver_fluid_sphcoord
  end interface
end module interface_poisson_solver_fluid_sphcoord
