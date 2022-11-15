module interface_poisson_solver
  implicit none
  interface 
    subroutine poisson_solver(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver
  end interface
end module interface_poisson_solver
