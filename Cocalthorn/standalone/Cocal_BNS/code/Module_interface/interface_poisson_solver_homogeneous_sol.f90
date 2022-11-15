module interface_poisson_solver_homogeneous_sol
  implicit none
  interface 
    subroutine poisson_solver_homogeneous_sol(surp, vpot_b)
      real(8),pointer  ::  surp(:,:), vpot_b(:,:,:)
    end subroutine  poisson_solver_homogeneous_sol
  end interface
end module interface_poisson_solver_homogeneous_sol
