module interface_poisson_solver_no_lm1
  implicit none
  interface 
    subroutine poisson_solver_no_lm1(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_no_lm1
  end interface
end module interface_poisson_solver_no_lm1
