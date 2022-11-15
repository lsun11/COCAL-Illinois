module interface_poisson_solver_no_dipole
  implicit none
  interface 
    subroutine poisson_solver_no_dipole(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_no_dipole
  end interface
end module interface_poisson_solver_no_dipole
