module interface_poisson_solver_Az
  implicit none
  interface 
    subroutine poisson_solver_Az(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_Az
  end interface
end module interface_poisson_solver_Az
